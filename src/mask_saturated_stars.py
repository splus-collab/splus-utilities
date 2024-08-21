#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script is part of the set of tools used to prepare the S-PLUS data
# Current version: prepared by Herpich F. R., March, 2024 - fabio.herpich@ast.cam.ac.uk
# Original version: prepared by Luisa Buzzo, 2021, which was used to produce
# the masks for the S-PLUS DR3

import logging
import colorlog
import os
import sys
import argparse
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits
from astropy.wcs import WCS
from regions import CirclePixelRegion, PixCoord
from astropy.wcs.utils import skycoord_to_pixel as sky2pix
from astropy.stats import sigma_clipped_stats
from astroquery.vizier import Vizier
import pandas as pd
import matplotlib.pyplot as plt
import glob
from itertools import repeat
from multiprocessing import Pool, freeze_support
import gc

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    'font.size': 20,
})


def get_args():
    parser = argparse.ArgumentParser(
        description='Apply masks to the SPLUS catalogues')
    parser.add_argument('--workdir', type=str, default=os.getcwd(),
                        help='Working directory')
    parser.add_argument('--datadir', type=str, default=None,
                        help='Directory with the SPLUS catalogues')
    parser.add_argument('--imgdir', type=str, default=None,
                        help='Directory with the SPLUS images')
    parser.add_argument('--outdir', type=str, default=None,
                        help='Output directory')
    parser.add_argument('--splusfootprint', type=str,
                        help='S-PLUS footprint')
    parser.add_argument('--cat_stars', type=str,
                        help='Catalogue of stars. If None, will query specified by ref2use')
    parser.add_argument('--ref2use', type=str, default='gsc',
                        help='Reference catalogue to use')
    parser.add_argument('--field', type=str,
                        help='Field to process')
    parser.add_argument('--listfields', type=str,
                        help='List of fields to process')
    parser.add_argument('--prefix', type=str, default='',
                        help='Prefix for the S-PLUS catalogue')
    parser.add_argument('--postfix', type=str, default='',
                        help='Postfix for the S-PLUS catalogue')
    parser.add_argument('--nprocs', type=int, default=1,
                        help='Number of processes to use')
    parser.add_argument('--plotstars', action='store_true',
                        help='Plot stars')
    parser.add_argument('--showfig', action='store_true',
                        help='Show figure')
    parser.add_argument('--savefig', action='store_true',
                        help='Save figure')
    parser.add_argument('--istest', action='store_true',
                        help='Test mode. Runs only one field')

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def call_logger():
    """Configure the logger."""
    logging.shutdown()
    logging.root.handlers.clear()

    # configure the module with colorlog
    logger = colorlog.getLogger()
    logger.setLevel(logging.INFO)

    # create a formatter with green color for INFO level
    formatter = colorlog.ColoredFormatter(
        '%(log_color)s%(levelname)s:%(name)s:%(message)s',
        log_colors={
            'DEBUG':    'cyan',
            'INFO':     'blue',
            'WARNING':  'yellow',
            'ERROR':    'red',
            'CRITICAL': 'red,bg_white',
        })

    # create handler and set the formatter
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)

    # add the handler to the logger
    logger.addHandler(ch)


def get_splusfootprint(
    args: argparse.Namespace
):
    """
    Get the S-PLUS footprint

    Parameters
    ----------
    args: argparse.Namespace
        Arguments

    Returns
    -------
    footprint: pd.DataFrame
        S-PLUS footprint
    """
    footppath = os.path.join(args.splusfootprint)
    try:
        footprint = pd.read_csv(footppath, delimiter=',')
    except FileNotFoundError:
        raise FileNotFoundError('S-PLUS footprint file not found...')

    gc.collect()

    return footprint


def query_gsc(
    ra: float,
    dec: float,
    radius: float = 1.0
):
    """
    Query the GSC1.2 catalog using the Vizier service

    Parameters
    ----------
    ra: float
        Right ascension in degrees
    dec: float
        Declination in degrees
    radius: float
        Radius in degrees

    Returns
    -------
    result: pd.DataFrame
        Result of the query
    """

    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'Pmag', 'e_Pmag', 'PosErr',
               'n_Pmag', 'Epoch', 'Class'],
               column_filters={"Pmag": "<14"},
               catalog="I/254")
    v.ROW_LIMIT = 10000
    result = v.query_region(SkyCoord(ra=ra, dec=dec, unit=(
        u.hour, u.deg), frame='icrs'),
        radius=Angle(radius, "deg"))[0]

    return result[result['Class'] == 0]


def get_stars(
    gsccat: pd.DataFrame,
    image: str = 'None',
    scatname: str = 'None',
    fieldname: str = 'None',
):
    """
    Get the stars from the GSC1.2 catalog and the S-PLUS catalogue

    Parameters
    ----------

    gsccat: pd.DataFrame
        GSC1.2 catalog
    image: str
        Image to plot
    spluscat: pd.DataFrame
        S-PLUS catalogue
    fieldname: str
        Field name

    Returns
    -------
    objects2plot: dict
        Dictionary with the objects to plot
    """
    if image == 'None':
        logging.warning('No image provided. Skipping...')
        f = None
    else:
        logging.info('Processing image %s' % image)
        try:
            f = fits.open(image)
        except FileNotFoundError:
            logging.error('Image not found. Skipping...')
            image = 'None'
            f = None
    spluscat = fits.open(scatname)[1].data
    gsc_coords = SkyCoord(ra=gsccat['RAJ2000'].value.data,
                          dec=gsccat['DEJ2000'].value.data,
                          unit=(u.deg, u.deg), frame='icrs')
    mask = np.ones(len(spluscat), dtype=bool)
    if ('ALPHA_J2000' in spluscat.columns.names) and ('DELTA_J2000' in spluscat.columns.names):
        sra = spluscat['ALPHA_J2000']
        sdec = spluscat['DELTA_J2000']
    elif ('RA' in spluscat.columns.names) and ('DEC' in spluscat.columns.names):
        sra = spluscat['RA']
        sdec = spluscat['DEC']
    else:
        raise ValueError(
            'No coordinates columns found in S-PLUS catalogue')
    splus_coords = SkyCoord(ra=sra,
                            dec=sdec,
                            unit=(u.deg, u.deg), frame='icrs')

    if f is not None:
        mean, median, std = sigma_clipped_stats(f[1].data, sigma=3.0)
        print('mean, median, std:', mean, median, std)
        wcs = WCS(f[1].header)
        pixcoords = np.transpose(sky2pix(gsc_coords, wcs))
        spixcoords = sky2pix(splus_coords, wcs)
    else:
        mean, median, std = None, None, None
        wcs = None
        pixcoords = None
        spixcoords = None

    gscrad = 19223 * \
        np.e ** (-(gsccat['Pmag'] - (-11.1)) ** 2/(2 * 6.2 ** 2))
    gscrad[gscrad < 30] = 30

    objects2plot = {'fieldname': fieldname,
                    'catname': scatname,
                    'image': f,
                    'wcs': wcs,
                    'std': std,
                    'gsc': {'coords': gsc_coords,
                            'rad': gscrad,
                            'pixcoords': pixcoords,
                            'colour': 'red',
                            'mag': gsccat['Pmag'].value.data},
                    }
    objects2plot['splus'] = {'ID': spluscat['ID'],
                             'coords': splus_coords,
                             'rad': 5,
                             'pixcoords': spixcoords,
                             'colour': 'blue',
                             'mask': mask}

    objects2plot['masksat'] = None
    gc.collect()

    return objects2plot


def plot_stars(
    args: argparse.Namespace,
    objects2plot: dict,
):
    """
    Plot the stars on the image

    Parameters
    ----------
    args: argparse.Namespace
        Arguments
    objects2plot: dict
        Dictionary with the objects to plot
    """

    image2plot = os.path.join(
        args.imgdir,  objects2plot['fieldname'] + '_trilogy.png')
    if not os.path.exists(image2plot):
        logging.info('No colour image found. Trying to find fits image...')
        image2plot = objects2plot['fieldname'] + '_R_swp.fz'
        if not os.path.exists(image2plot):
            logging.warning('No image found. Skipping...')
            return
        else:
            imgdata = objects2plot['image'][1].data
            imgdata[imgdata < objects2plot['std']] = 0
    else:
        imgdata = np.flip(plt.imread(image2plot), axis=0)

    outputname = os.path.join(
        objects2plot['catname'].strip('.fits') + '_masked.png')

    if not os.path.exists(outputname):
        wcs = objects2plot['wcs']
        pixcoords = objects2plot['gsc']['pixcoords']
        gscrad = objects2plot['gsc']['rad']
        gscregions = [CirclePixelRegion(center=PixCoord(x, y), radius=z)
                      for (x, y), z in zip(pixcoords, gscrad)]

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection=wcs)
        ax.imshow(imgdata, origin='lower', vmin=0, vmax=3.5)
        logging.info('Plotting stars...')
        for reg, mag in zip(gscregions, objects2plot['gsc']['mag']):
            reg.plot(ax=ax, color='c', lw=3)
            ax.text(reg.center.x, reg.center.y,
                    f'{mag:.2f}', color='white', fontsize=8)
        if 'splus' in objects2plot and objects2plot['masksat'] is None:
            spixcoords = objects2plot['splus']['pixcoords']
            mask = objects2plot['splus']['mask']
            splusregions = [CirclePixelRegion(center=PixCoord(x, y), radius=3)
                            for x, y in zip(spixcoords[0][mask], spixcoords[1][mask])]
            for reg in splusregions:
                reg.plot(ax=ax, color='blue', lw=3)
        elif 'splus' in objects2plot and objects2plot['masksat'] is not None:
            spixcoords = objects2plot['splus']['pixcoords']
            mask = objects2plot['masksat']
            splusregions = [CirclePixelRegion(center=PixCoord(x, y), radius=3)
                            for x, y in zip(spixcoords[0], spixcoords[1])]
            for l, reg in zip(mask, splusregions):
                if l == 1:
                    reg.plot(ax=ax, color='r', lw=3)
                elif l == 2:
                    reg.plot(ax=ax, color='y', lw=3)
                elif l == 3:
                    reg.plot(ax=ax, color='orange', lw=3)
                else:
                    reg.plot(ax=ax, color='blue', lw=3)
        elif 'splus' not in objects2plot:
            logging.info('No S-PLUS catalogue provided')
        else:
            logging.info('No mask provided')

        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')

        if args.savefig and args.showfig:
            logging.info(f'Saving figure to {outputname}')
            plt.savefig(outputname, dpi=300)
            plt.show()
        elif args.savefig and not args.showfig:
            logging.info(f'Saving figure to {outputname}')
            plt.savefig(outputname, dpi=300)
            plt.close()
        elif not args.savefig and args.showfig:
            plt.show()
        else:
            logging.info(f'Closing figure...')
            plt.close()
    else:
        print('File already exists. Skipping...')

    gc.collect()


def make_masks(
    args: argparse.Namespace,
    objects2plot: dict,
):
    """
    Make the masks for the saturated stars

    Parameters
    ----------
    args: argparse.Namespace
        Arguments
    objects2plot: dict

    Returns
    -------
    objects2plot: dict
        Dictionary with the objects to plot
    """

    if 'splus' in objects2plot:
        scoords = objects2plot['splus']['coords']
    else:
        raise ValueError('No S-PLUS catalogue provided')
    mask = np.zeros(len(scoords), dtype=int)
    for coord, radius in zip(objects2plot['gsc']['coords'], objects2plot['gsc']['rad']):
        separation = coord.separation(scoords)
        mask[separation < radius * 0.55 * u.arcsec] = 1

    objects2plot['masksat'] = mask

    # if args.showfig or args.savefig:
    #     plot_stars(args, objects2plot)

    gc.collect()

    return objects2plot


def check_distance_to_border(
    objects2plot: dict,
    sfoot: pd.DataFrame,
    fieldname: str,
):
    """
    Check the distance of the objects to the border of the field

    Parameters
    ----------
    objects2plot: dict
        Dictionary with the objects to plot
    sfoot: pd.DataFrame
        S-PLUS footprint
    fieldname: str
        Field name

    Returns
    -------
    objects2plot: dict
        Dictionary with the objects to plot
    """
    fields = np.array([f.replace('_', '-') for f in sfoot['NAME']])
    field_coords = sfoot[fields == fieldname]
    field_coords = SkyCoord(ra=field_coords['RA'],
                            dec=field_coords['DEC'],
                            unit=(u.deg, u.deg))
    scoords = objects2plot['splus']['coords']
    if (len(scoords) == 0) or (len(field_coords) == 0):
        objects2plot['masksat'] = None
    else:
        sepxy = field_coords.spherical_offsets_to(scoords)
        mask = abs(sepxy[0]) > 0.7 * u.deg - 30 * u.arcsec
        mask |= abs(sepxy[1]) > 0.7 * u.deg - 30 * u.arcsec
        objects2plot['masksat'][mask] += 2

    gc.collect()

    return objects2plot


def process_field(
    args: argparse.Namespace,
    sfoot: pd.DataFrame,
    fieldname: str,
):
    """
    Process the field

    Parameters
    ----------
    args: argparse.Namespace
        Arguments
    sfoot: pd.DataFrame
        S-PLUS footprint
    fieldname: str
        Field name
    """
    try:
        cats2proc = glob.glob(os.path.join(
            args.datadir, f'{args.prefix}{fieldname}{args.postfix}.fits'))
    except FileNotFoundError:
        f = open(os.path.join(args.datadir, f'{fieldname}_failed.txt'), 'w')
        f.write('No S-PLUS catalog found for field: %s' % fieldname)
        f.close()
        return
    if len(cats2proc) == 0:
        logging.warning('No S-PLUS catalog found for field: %s' % fieldname)
        return
    imagename = os.path.join(args.imgdir, f'{fieldname}_R_swp.fz')
    field_coords = sfoot[sfoot['NAME'] == fieldname.replace('-', '_')]
    gsccat = query_gsc(field_coords['RA'], field_coords['DEC'])
    for catname in cats2proc:
        outcat = os.path.join(
            args.outdir, f'{catname.split("/")[-1].replace(".fits", "_mask.csv")}')
        failed_cat = outcat.replace('_mask.csv', '_failed.txt')
        if not os.path.exists(outcat) and not os.path.exists(failed_cat):
            objects2plot = get_stars(
                gsccat, image=imagename, scatname=catname, fieldname=fieldname)
            objects2plot = make_masks(args, objects2plot)
            objects2plot = check_distance_to_border(
                objects2plot, sfoot, fieldname)
            if objects2plot['masksat'] is not None:
                newdf = pd.DataFrame()
                newdf['ID'] = objects2plot['splus']['ID']
                newdf['RA'] = objects2plot['splus']['coords'].ra.value
                newdf['Dec'] = objects2plot['splus']['coords'].dec.value
                newdf['MASK'] = objects2plot['masksat']
                logging.info(f'Writing mask catalogue to disk: {outcat}')
                newdf.to_csv(outcat, index=False)
                if args.plotstars:
                    logging.info('Plotting stars...')
                    plot_stars(args, objects2plot)
                else:
                    logging.info('No plot requested. Skipping...')
            else:
                with open(outcat.replace('_mask.csv', '_failed.txt'), 'w') as f:
                    f.write('Failed to calculate masks for file %s' % catname)
                    f.close()
                logging.error(
                    'Failed to calculate masks for file %s' % catname)
        else:
            if args.plotstars:
                objects2plot = get_stars(
                    gsccat, image=imagename, scatname=catname, fieldname=fieldname)
                objects2plot = make_masks(args, objects2plot)
                logging.info('Plotting stars...')
                plot_stars(args, objects2plot)
            else:
                logging.info('No plot requested. Skipping...')
            logging.info('File already exists. Skipping...')

    gc.collect()

    return


def main():

    args = get_args()
    if args.datadir is None:
        logging.info('No data directory provided. Using workdir instead...')
        args.datadir = args.workdir
    if args.imgdir is None:
        logging.info('No image directory provided. Using workdir instead...')
        args.imgdir = args.workdir
    if args.outdir is None:
        logging.info('No output directory provided. Using workdir instead...')
        args.outdir = args.workdir
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if (args.showfig or args.savefig) and not args.plotstars:
        logging.warning('No plot requested. Will not show nor save figures')
    sfoot = get_splusfootprint(args)
    if args.field is not None:
        fieldname = args.field.replace(
            '_', '-')
        process_field(args, sfoot, fieldname)
    elif args.listfields is not None:
        f = pd.read_csv(args.listfields)
        try:
            fields = f['NAME']
        except KeyError:
            fields = f['Field']
        fields = [field.replace('_', '-') for field in fields]
        if args.istest:
            fields = fields[:6]
            nprocs = 3
        else:
            nprocs = args.nprocs
        pool = Pool(nprocs)
        try:
            pool.starmap(process_field, zip(
                repeat(args), repeat(sfoot), fields))
            pool.close()
            pool.join()
        except Exception as e:
            print(e)
            pool.close()
            pool.join()
            log_file = 'pool_exception_fail.log'
            if not os.path.exists(log_file):
                with open(log_file, 'w') as f:
                    f.write(str(e))
                    f.close()
            else:
                with open(log_file, 'a') as f:
                    f.write(str(e))
                    f.close()

    else:
        raise ValueError('No field or list of fields provided')

    gc.collect()

    return


if __name__ == '__main__':
    call_logger()
    freeze_support()
    main()
