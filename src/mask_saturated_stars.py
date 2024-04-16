#!/usr/bin/env python
# -*- coding: utf-8 -*-
# This script is part of the set of tools used to prepare the S-PLUS data
# Current version: prepared by Herpich F. R., March, 2024 - fabio.herpich@ast.cam.ac.uk
# Original version: prepared by Luisa Buzzo, 2021, which was used to produce
# the masks for the S-PLUS DR3

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
import warnings

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
    args, _ = parser.parse_known_args()
    parser.add_argument('--datadir', type=str, default=args.workdir,
                        help='Directory with the SPLUS catalogues')
    parser.add_argument('--imgdir', type=str, default=args.workdir,
                        help='Directory with the SPLUS images')
    parser.add_argument('--outdir', type=str, default=args.workdir,
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
    parser.add_argument('--spluscat', type=str,
                        help='S-PLUS catalogue prefix to serach for')
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

    Returns
    -------
    objects2plot: dict
        Dictionary with the objects to plot
    """
    if image == 'None':
        raise ValueError('No image provided')
    try:
        f = fits.open(image)
    except FileNotFoundError:
        raise FileNotFoundError('Image not found')
    spluscat = fits.open(scatname)[1].data
    mean, median, std = sigma_clipped_stats(f[1].data, sigma=3.0)
    print('mean, median, std:', mean, median, std)
    wcs = WCS(f[1].header)
    gsc_coords = SkyCoord(ra=gsccat['RAJ2000'].value.data,
                          dec=gsccat['DEJ2000'].value.data,
                          unit=(u.deg, u.deg), frame='icrs')
    pixcoords = np.transpose(sky2pix(gsc_coords, wcs))
    gscrad = 19223 * \
        np.e ** (-(gsccat['Pmag'] - (-11.1)) ** 2/(2 * 6.2 ** 2))
    gscrad[gscrad < 30] = 30

    objects2plot = {'catname': scatname,
                    'image': f,
                    'wcs': wcs,
                    'std': std,
                    'gsc': {'coords': gsc_coords,
                            'rad': gscrad,
                            'pixcoords': pixcoords,
                            'colour': 'red',
                            'mag': gsccat['Pmag'].value.data},
                    }
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
    spixcoords = sky2pix(splus_coords, wcs)
    objects2plot['splus'] = {'coords': splus_coords,
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
    Plot the stars in the image

    Parameters
    ----------
    args: argparse.Namespace
        Arguments
    objects2plot: dict
        Dictionary with the objects to plot
    """

    if 'splus' in objects2plot and objects2plot['masksat'] is None:
        outputname = os.path.join(
            args.workdir, objects2plot['catname'].strip('.fits') + '_splus.png')
    elif 'splus' in objects2plot and objects2plot['masksat'] is not None:
        outputname = os.path.join(
            args.workdir, objects2plot['catname'].strip('.fits') + '_masked.png')
    elif 'splus' not in objects2plot:
        outputname = os.path.join(
            args.workdir, objects2plot['catname'].strip('.fits') + '_gsc.png')
    else:
        outputname = os.path.join(
            args.workdir, objects2plot['catname'].strip('.fits') + '_nomask.png')
    if not os.path.exists(outputname):
        data = objects2plot['image'][1].data
        wcs = objects2plot['wcs']
        pixcoords = objects2plot['gsc']['pixcoords']
        gscrad = objects2plot['gsc']['rad']
        gscregions = [CirclePixelRegion(center=PixCoord(x, y), radius=z)
                      for (x, y), z in zip(pixcoords, gscrad)]

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection=wcs)
        data[data < objects2plot['std']] = 0
        ax.imshow(data, origin='lower', vmin=0, vmax=3.5)
        print('Plotting stars')
        for reg, mag in zip(gscregions, objects2plot['gsc']['mag']):
            reg.plot(ax=ax, color='c', lw=3)
            ax.text(reg.center.x, reg.center.y, f'{mag:.2f}', color='white')
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
            masked_regs = [CirclePixelRegion(center=PixCoord(x, y), radius=3)
                           for x, y in zip(spixcoords[0][mask == 1], spixcoords[1][mask == 1])]
            for reg in masked_regs:
                reg.plot(ax=ax, color='r', lw=3)
            unmasked_regs = [CirclePixelRegion(center=PixCoord(x, y), radius=3)
                             for x, y in zip(spixcoords[0][mask == 0], spixcoords[1][mask == 0])]
            for reg in unmasked_regs:
                reg.plot(ax=ax, color='forestgreen', lw=3)
        elif 'splus' not in objects2plot:
            print('No S-PLUS catalog provided')
        else:
            print('No mask provided')

        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')

        if args.savefig:
            plt.savefig(outputname, dpi=300)
            if args.showfig:
                plt.show()
            else:
                plt.close()
        if args.showfig:
            plt.show()
        else:
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
        raise ValueError('No S-PLUS catalog provided')
    mask = np.zeros(len(scoords), dtype=int)
    for coord, radius in zip(objects2plot['gsc']['coords'], objects2plot['gsc']['rad']):
        separation = coord.separation(scoords)
        mask[separation < radius * 0.55 * u.arcsec] = 1

    objects2plot['masksat'] = mask

    if args.showfig or args.savefig:
        plot_stars(args, objects2plot)

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
    objects2plot['masksat'][(abs(scoords.ra - field_coords.ra + 0.7 * u.deg) < 30 * u.arcsec) |
                            (abs(scoords.dec - field_coords.dec + 0.7 * u.deg) < 30 * u.arcsec)] += 2

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
    cats2proc = glob.glob(os.path.join(args.datadir, f'{fieldname}_*.fits'))
    if len(cats2proc) == 0:
        warnings.warn('No S-PLUS catalog found for field: %s' % fieldname)
        return
    imagename = os.path.join(args.imgdir, fieldname, f'{fieldname}_R_swp.fz')
    field_coords = sfoot[sfoot['NAME'] == fieldname]
    gsccat = query_gsc(field_coords['RA'], field_coords['DEC'])
    for catname in cats2proc:
        outcat = os.path.join(
            args.outdir, f'{catname.split("/")[-1].replace(".fits", "_mask.csv")}')
        if not os.path.exists(outcat):
            objects2plot = get_stars(gsccat, image=imagename, scatname=catname)
            if args.plotstars:
                plot_stars(args, objects2plot)
            objects2plot = make_masks(args, objects2plot)
            objects2plot = check_distance_to_border(
                objects2plot, sfoot, fieldname)
            newdf = pd.DataFrame()
            newdf['RA'] = objects2plot['splus']['coords'].ra.value
            newdf['Dec'] = objects2plot['splus']['coords'].dec.value
            newdf['MASK'] = objects2plot['masksat']
            print('Writing mask catalogue to disk:', outcat)
            newdf.to_csv(outcat)
        else:
            print('File already exists. Skipping...')

    gc.collect()

    return


def main():

    args = get_args()
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if (args.showfig or args.savefig) and not args.plotstars:
        warnings.warn('No plot requested. Will not show nor save figures')
    sfoot = get_splusfootprint(args)
    if args.field is not None:
        fieldname = args.field.replace(
            '_', '-')
        process_field(args, sfoot, [fieldname])
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
        pool.starmap(process_field, zip(repeat(args), repeat(sfoot), fields))
        pool.close()
        pool.join()
    else:
        raise ValueError('No field or list of fields provided')

    gc.collect()

    return


if __name__ == '__main__':
    freeze_support()
    main()
