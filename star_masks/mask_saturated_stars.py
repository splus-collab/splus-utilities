#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
import astropy.units as u
# from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import SkyCoord, Angle
# from astropy.coordinates import search_around_sky
# from astropy.table import Table, vstack, hstack, unique
from astropy.io import fits
from astropy.wcs import WCS
from regions import CirclePixelRegion, PixCoord
from astropy.wcs.utils import skycoord_to_pixel as sky2pix
from photutils import DAOStarFinder
from astropy.stats import sigma_clipped_stats
# from skyplot import *
# import astropy.coordinates as coord
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from astroquery.vizier import Vizier
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
    'font.size': 20,
})


def get_args():
    parser = argparse.ArgumentParser(
        description='Apply masks to the SPLUS catalogs')
    parser.add_argument('--workdir', type=str, default=os.getcwd(),
                        help='Working directory')
    parser.add_argument('--datadir', type=str, default=os.getcwd(),
                        help='Directory with the SPLUS catalogs')
    parser.add_argument('--splusfootprint', type=str,
                        help='S-PLUS footprint')
    parser.add_argument('--cat_stars', type=str,
                        help='Catalog of stars')
    parser.add_argument('--catfile', type=str,
                        help='Catalog of objects')
    parser.add_argument('--photoflag', type=str,
                        help='Photometric flag')
    parser.add_argument('--field', type=str,
                        help='Field name')
    parser.add_argument('--ref2use', type=str, default='gsc',
                        help='Reference catalog to use')
    parser.add_argument('--spluscat', type=str,
                        help='S-PLUS catalog')
    args, _ = parser.parse_known_args()
    parser.add_argument('--outdir', type=str, default=args.workdir,
                        help='Output directory')

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def get_splusfootprint(args):
    """
    Get the S-PLUS footprint
    """
    footppath = os.path.join(args.splusfootprint)
    footprint = pd.read_csv(footppath, delimiter=',')

    return footprint


def query_gsc(ra, dec, radius=1.0):
    """
    Query the GSC1.2 catalog using the Vizier service
    """

    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'Pmag', 'e_Pmag', 'PosErr',
               'n_Pmag', 'Epoch', 'Class'],
               column_filters={"Pmag": "<12"},
               catalog="I/254")
    v.ROW_LIMIT = 10000
    result = v.query_region(SkyCoord(ra=ra, dec=dec, unit=(
        u.hour, u.deg), frame='icrs'),
        radius=Angle(radius, "deg"))[0]

    return result[result['Class'] == 0]


def plot_stars(gsccat, image, spluscat=None, showfig=False, savefig=False):
    f = fits.open(image)
    mean, median, std = sigma_clipped_stats(f[1].data, sigma=3.0)
    print('mean, median, std:', mean, median, std)
    wcs = WCS(f[1].header)
    gsc_coords = SkyCoord(ra=gsccat['RAJ2000'].value.data,
                          dec=gsccat['DEJ2000'].value.data,
                          unit=(u.deg, u.deg), frame='icrs')
    pixcoords = np.transpose(sky2pix(gsc_coords, wcs))
    gscrad = 23979.31 * \
        np.e ** (-(gsccat['Pmag'] - (-10.58363)) ** 2/(2 * 5.5 ** 2))
    gscrad[gscrad < 50] = 50

    if showfig:
        gscregions = [CirclePixelRegion(center=PixCoord(x, y), radius=z)
                      for (x, y), z in zip(pixcoords, gscrad)]

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection=wcs)
        f[1].data[f[1].data < std] = 0
        ax.imshow(f[1].data, origin='lower', vmin=0, vmax=3.5)
        print('Plotting stars')
        for reg, mag in zip(gscregions, gsccat['Pmag']):
            reg.plot(ax=ax, color='red', lw=3)
            ax.text(reg.center.x, reg.center.y, f'{mag:.2f}', color='white')
        if spluscat is not None:
            mask = spluscat['MAG_AUTO'] < 16
            # mask &= spluscat['FLAGS'] != 0
            # mask &= spluscat['CLASS_STAR'] > 0.7
            splus_coords = SkyCoord(ra=spluscat['ALPHA_J2000'][mask],
                                    dec=spluscat['DELTA_J2000'][mask],
                                    unit=(u.deg, u.deg), frame='icrs')
            spixcoords = sky2pix(splus_coords, wcs)
            splusregions = [CirclePixelRegion(center=PixCoord(x, y), radius=3)
                            for x, y in zip(spixcoords[0], spixcoords[1])]
            for reg in splusregions:
                reg.plot(ax=ax, color='blue', lw=3)

        ax.set_xlabel('RA')
        ax.set_ylabel('Dec')
        plt.tight_layout()

        if savefig:
            outputname = os.path.join(args.workdir,
                                      image.split('/')[-1].replace('swp.fz', '_gsc.png'))
            plt.savefig(outputname, dpi=300)
        plt.show()

    return gsc_coords, gscrad


def run_daofinder(image, fwhm=50.0, threshold=1.0):
    """
    Run the DAOStarFinder algorithm to find stars in the image
    """
    f = fits.open(image)
    data = f[1].data
    mean, median, std = sigma_clipped_stats(data, sigma=5.0)
    # median, std = 0.0, 0.5
    # daofind = DAOStarFinder(fwhm=fwhm, sharplo=0.2, sharphi=0.9,
    #                         roundlo=-0.5, roundhi=0.5, threshold=threshold*std)
    daofind = DAOStarFinder(
        fwhm=fwhm, threshold=threshold*std, ratio=0.5, theta=90.0)
    sources = daofind(data - median)
    norm_flux = sources['flux'] / np.max(sources['flux'])
    mask = norm_flux > 0.5
    import pdb
    pdb.set_trace()
    return sources[mask]


def plot_daostars(image, sources):
    f = fits.open(image)
    wcs = WCS(f[1].header)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=wcs)
    ax.imshow(f[1].data, origin='lower', vmin=-0.1, vmax=3.5)
    daosources = np.transpose([sources['xcentroid'], sources['ycentroid']])
    daorad = 4 * (abs(sources['sharpness']) +
                  abs(sources['roundness1']) + abs(sources['roundness2']))
    daoregions = [CirclePixelRegion(center=PixCoord(x, y), radius=r)
                  for (x, y), r in zip(daosources, daorad)]
    for reg in daoregions:
        reg.plot(ax=ax, color='red', lw=3)
    plt.show()


def make_masks(coords, rad, image, spluscat, showfig=False, savefig=False):
    f = fits.open(image)
    wcs = WCS(f[1].header)
    scoords = SkyCoord(ra=spluscat['ALPHA_J2000'], dec=spluscat['DELTA_J2000'],
                       unit=(u.deg, u.deg), frame='icrs')
    mask = np.zeros(len(spluscat))
    for coord, radius in zip(coords, rad):
        separation = coord.separation(scoords)
        mask[separation < radius * 0.55 * u.arcsec] = 1

    if showfig:
        gscregions = [CirclePixelRegion(center=PixCoord(x, y), radius=z)
                      for (x, y), z in zip(np.transpose(sky2pix(coords, wcs)), rad)]
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection=wcs)
        ax.imshow(f[1].data, origin='lower', vmin=0, vmax=3.5)
        for reg in gscregions:
            reg.plot(ax=ax, color='red', lw=3)

        masked_regs = [CirclePixelRegion(center=PixCoord(x, y), radius=3)
                       for x, y in np.transpose(sky2pix(scoords[mask == 1], wcs))]
        for reg in masked_regs:
            reg.plot(ax=ax, color='blue', lw=5)
        unmasked_regs = [CirclePixelRegion(center=PixCoord(x, y), radius=3)
                         for x, y in np.transpose(sky2pix(scoords[mask == 0], wcs))]
        for reg in unmasked_regs:
            reg.plot(ax=ax, color='green', lw=5)
        if savefig:
            outputname = os.path.join(args.workdir,
                                      image.split('/')[-1].replace('swp.fz', '_masked.png'))
            plt.savefig(outputname, dpi=300)
        plt.show()

    return mask


def check_distance_to_border(mask, spluscat, sfoot, fieldname):
    """
    Check the distance of the objects to the border of the field
    """
    field_coords = sfoot[sfoot['NAME'] == fieldname]
    field_coords = SkyCoord(ra=field_coords['RA'], dec=field_coords['DEC'], unit=(
        u.deg, u.deg))
    scoords = SkyCoord(ra=spluscat['ALPHA_J2000'], dec=spluscat['DELTA_J2000'],
                       unit=(u.deg, u.deg), frame='icrs')
    mask[(abs(scoords.ra - field_coords.ra + 0.7 * u.deg) < 10 * u.arcsec) |
         (abs(scoords.dec - field_coords.dec + 0.7 * u.deg) < 10 * u.arcsec)] += 2
    return mask


def main():

    args = get_args()
    sfoot = get_splusfootprint(args)
    fieldname = args.field.replace(
        '-', '_') if 'STRIPE82' in args.field else args.field
    imagename = os.path.join(args.datadir, f'{args.field}_R_swp.fz')
    if args.ref2use in ['gsc', 'GSC']:
        splus_cat = fits.open(os.path.join(
            args.datadir, args.spluscat))[2].data
        field_coords = sfoot[sfoot['NAME'] == fieldname]
        gsccat = query_gsc(field_coords['RA'], field_coords['DEC'])
        gsc_coords, gsc_rad = plot_stars(gsccat, imagename, showfig=True)
        masked_bright_stars = make_masks(
            gsc_coords, gsc_rad, imagename, splus_cat)

        mask = check_distance_to_border(
            masked_bright_stars, splus_cat, sfoot, fieldname)
        newdf = pd.DataFrame()
        newdf['RA'] = splus_cat['ALPHA_J2000']
        newdf['Dec'] = splus_cat['DELTA_J2000']
        newdf['MASK'] = mask
        print('Writing masked catalog', os.path.join(
            args.outdir, f'{args.field}_mask.csv'))
        newdf.to_csv(os.path.join(args.outdir, f'{args.field}_mask.csv'))
    elif args.ref2use in ['dao', 'DAO', 'daofinder']:
        sources = run_daofinder(imagename)
        plot_daostars(imagename, sources)
    return


if __name__ == '__main__':
    main()
