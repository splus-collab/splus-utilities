#!/usr/bin/env python3
# Herpich F.R. -- 2019-12-13
# Plots the S-PLUS footprint on the sky
# using the tiles_nc.csv file
#
# Usage: python3 plot_footprint.py --splusfoot tiles_nc.csv
#
# Output: splus_footprint.png
#
# Notes: The 2MASS catalog is available at
# http://cdsarc.u-strasbg.fr/viz-bin/Cat?II/281
#
# Last modified: 2023-07-19
# Modified by: Herpich F.R. CASU/IoA Cambridge
# email: fabio.herpich@ast.cam.ac.uk
# --------------------------------------------------

import os
import sys
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
import numpy as np
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description='Plots the S-PLUS footprint on the sky')
    parser.add_argument('--splusfoot', type=str, required=True,
                        help='S-PLUS footprint file')
    parser.add_argument('--workdir', type=str, default=os.getcwd(),
                        help='Working directory')
    parser.add_argument('--is2mass', action='store_true',
                        help='If 2mass allssky catalogue is available')
    parser.add_argument('--twomasscat', type=str, default=None,
                        help='2MASS allssky catalogue name')
    # print help if no arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def plot_foot(args):
    plt.figure(figsize=(16, 8.4))
    ax = plt.subplot(111, projection="aitoff")
    plt.grid(True)

    if args.is2mass:
        if args.twomasscat is None:
            print('If 2mass allssky catalogue is available, please provide the catalogue name')
            raise ValueError('Please provide the 2mass allssky catalogue name')
        print('Using 2MASS catalogue %s' % args.twomasscat)
        try:
            cat = fits.open(args.twomasscat)[1].data
        except (UnicodeDecodeError, OSError):
            cat = ascii.read(args.twomasscat, format='csv')
        else:
            raise ValueError('Please provide a valid 2mass allssky catalogue on FITS or CSV format')

        c0 = SkyCoord(ra=cat['RAJ2000'], dec=cat['DEJ2000'],
                      unit=(u.deg, u.deg), frame='icrs')
        ra_rad0 = c0.ra.wrap_at(180 * u.deg).radian
        dec_rad0 = c0.dec.radian
        h, ex, ey = np.histogram2d(ra_rad0, dec_rad0, bins=(500, 500))
        xpos, ypos = np.meshgrid(ex[:-1], ey[:-1])
        h[h < 1] = 1.

        print('Ploting 2MASS...')
        ax.scatter(xpos, ypos, c=np.log10(h.T), s=10,
                   marker='H', edgecolor='None', cmap='Greys')

    print('Reading splus table...')
    t = ascii.read(args.splusfoot, format='csv', fast_reader=False)

    c = SkyCoord(ra=t['RA'], dec=t['DEC'], unit=(u.hour, u.deg), frame='icrs')
    ra_rad = c.ra.wrap_at(180 * u.deg).radian
    dec_rad = c.dec.radian

    mask_obs = (t['STATUS'] == 1) | (t['STATUS'] == 2) | (
        t['STATUS'] == 4) | (t['STATUS'] == 5) | (t['STATUS'] == 6)
    lbl = r'$\mathrm{Observed:\ %i}$' % mask_obs.sum()
    ax.scatter(ra_rad[mask_obs], dec_rad[mask_obs],
               marker='H', s=8, color='limegreen', label=lbl)
    mask_foo = (t['STATUS'] == -1) | (t['STATUS'] == -
                                      2) | (t['STATUS'] == 0) | (t['STATUS'] == 3)
    lbl = r'$\mathrm{To\ be\ observed:\ %i}$' % mask_foo.sum()
    ax.scatter(ra_rad[mask_foo], dec_rad[mask_foo],
               marker='H', s=8, color='r', label=lbl)

    plt.setp(ax.get_xticklabels(), fontsize=18)
    plt.setp(ax.get_yticklabels(), fontsize=18)
    plt.legend(loc='upper right', scatterpoints=1, markerscale=3, shadow=True,
               bbox_to_anchor=[1.02, 1.07], fancybox=True, fontsize=20)
    plt.subplots_adjust(top=0.95, bottom=0.05, right=0.9, left=0.05)
    plt.grid(True)

    pathtosave = './splus_footprint.png'
    print('Saving fig %s...' % pathtosave)
    plt.savefig(pathtosave, format='png', dpi=300)
    plt.show()

if __name__ == '__main__':
    args = parse_args()
    plot_foot(args)
