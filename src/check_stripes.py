#!/urs/bin/env python3

import os
import sys
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from astroquery.vizier import Vizier


def get_args():
    parser = argparse.ArgumentParser(
        description='Check for stripes in the image')
    parser.add_argument('--workdir', type=str,
                        default=os.getcwd(), help='The working directory')
    parser.add_argument('--filename', type=str,
                        help='The filename of the image to check')
    parser.add_argument('--radius', type=float, default=1.0,
                        help='The radius to search for Gaia stars')
    parser.add_argument('--splusfoot', type=str,
                        help='The filename of the SPLUS footprint')
    parser.add_argument('--sfield', type=str,
                        help='The SPLUS field to search for')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def get_splus_foot(args):
    df = pd.read_csv(args.splusfoot, delimiter=',')
    return df


def get_gaia_stars(
    coords: pd.DataFrame,
    radius: float = 1.0
):

    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'Gmag'])
    v.ROW_LIMIT = 100000
    result = v.query_region(coords, radius=f'{radius}d', catalog='I/355')

    return result[0]


def match_gaia_splus(
    args: argparse.Namespace,
    gaia_stars: pd.DataFrame,
    sfield: str = 'SPLUS-s28s25'
):
    with fits.open(os.path.join(args.workdir, sfield + '_dual_VAC_features.fits')) as hdul:
        data = hdul[1].data
        mask = data['SEX_FLAGS_DET'] >= 0
        coords = SkyCoord(ra=data['RA'][mask],
                          dec=data['DEC'][mask], unit=(u.deg, u.deg))
        splusdata = data[mask]

    gaia_coords = SkyCoord(ra=gaia_stars['RAJ2000'],
                           dec=gaia_stars['DEJ2000'], unit=(u.deg, u.deg))
    idx, d2d, d3d = coords.match_to_catalog_sky(gaia_coords)
    mask = d2d < 1 * u.arcsec

    gaia_matches = gaia_stars[idx[mask]]
    splus_matches = splusdata[mask]

    return data, gaia_stars, gaia_matches, splus_matches


def plot_gaia_splus(
    args: argparse.Namespace,
    splusdetect: pd.DataFrame,
    gaia_stars: pd.DataFrame,
    gaia_matches: pd.DataFrame,
    splus_matches: pd.DataFrame,
    fieldname: str = 'SPLUS-s28s25'
):
    mask = splus_matches['SEX_FLAGS_DET'] < 8
    # create a multiplot with a larger main panel and three smaller ones on a 3x1 scale
    left, width = 0.1, 0.5
    bottom, height = 0.05, 0.5
    spacing = 0.005
    small_width, small_height = 0.28, 0.3

    rect_scatter = [left, bottom + 0.35 + spacing, width, height]
    bottom_left_pad = [left, bottom, small_width, small_height]
    bottom_middle_pad = [left + small_width + spacing, bottom,
                         small_width, small_height]
    bottom_right_pad = [left + 2 * small_width + 2 * spacing, bottom,
                        small_width, small_height]

    fig = plt.figure(figsize=(10, 10))

    ax_scatter = fig.add_axes(rect_scatter)
    ax_scatter.scatter(splusdetect['RA'], splusdetect['DEC'],
                       s=10, c='gray', label='SPLUS all', alpha=0.5)
    ax_scatter.scatter(gaia_stars['RAJ2000'], gaia_stars['DEJ2000'],
                       s=5, c='g', label='Gaia all', alpha=0.5)
    ax_scatter.scatter(gaia_matches['RAJ2000'], gaia_matches['DEJ2000'],
                       s=20, c='r', label='Gaia mathces')
    ax_scatter.scatter(splus_matches['RA'][mask], splus_matches['DEC'][mask],
                       s=10, c='b', label=r"$\mathrm{SPLUS\ SF < 8}$")
    mask = splus_matches['SEX_FLAGS_DET'] == 0
    ax_scatter.scatter(splus_matches['RA'][mask], splus_matches['DEC'][mask],
                       s=5, c='y', label=r"$\mathrm{SPLUS\ DET = 0}$")
    ax_scatter.legend(loc='upper right', fontsize=16,
                      bbox_to_anchor=(1.7, 1.0), ncol=1)
    ax_scatter.set_title(fieldname)

    ax_bottom_left_pad = fig.add_axes(bottom_left_pad)
    ax_bottom_left_pad.scatter(splusdetect['r_auto'], splusdetect['e_r_auto'],
                               s=25, c='gray', label='SPLUS all', alpha=0.5)
    ax_bottom_left_pad.set_xlim(10, 25)
    ax_bottom_left_pad.set_ylim(0, 1)
    ax_bottom_left_pad.legend(loc='upper left', fontsize=8)

    mask = splusdetect['SEX_FLAGS_DET'] < 8
    ax_bottom_middle_pad = fig.add_axes(
        bottom_middle_pad, sharex=ax_bottom_left_pad, sharey=ax_bottom_left_pad)
    ax_bottom_middle_pad.scatter(splusdetect['r_auto'][mask], splusdetect['e_r_auto'][mask],
                                 s=25, c='b', label=r"$\mathrm{SPLUS\ SF < 8}$", alpha=0.5)
    ax_bottom_middle_pad.set_xlim(10, 25)
    ax_bottom_middle_pad.set_ylim(0, 1)
    plt.setp(ax_bottom_middle_pad.get_yticklabels(), visible=False)
    ax_bottom_middle_pad.legend(loc='upper left', fontsize=8)

    mask = splusdetect['SEX_FLAGS_DET'] == 0
    ax_bottom_right_pad = fig.add_axes(
        bottom_right_pad, sharex=ax_bottom_left_pad, sharey=ax_bottom_left_pad)
    ax_bottom_right_pad.scatter(splusdetect['r_auto'][mask], splusdetect['e_r_auto'][mask],
                                s=25, c='b', label=r"$\mathrm{SPLUS\ DET = 0}$", alpha=0.5)
    ax_bottom_right_pad.set_xlim(10, 25)
    ax_bottom_right_pad.set_ylim(0, 1)
    plt.setp(ax_bottom_right_pad.get_yticklabels(), visible=False)
    ax_bottom_right_pad.legend(loc='upper left', fontsize=8)

    plt.show()
    return {'spluscat': splusdetect,
            'gaiacat': gaia_stars,
            'gaiamatches': gaia_matches,
            'splusmatches': splus_matches}


def main(args):
    fields = ['SPLUS-s28s25']
    df = get_splus_foot(args)
    for field in fields:
        fielddf = df[df['NAME'].replace('_', '-') == field]
        coords = SkyCoord(ra=fielddf['RA'],
                          dec=fielddf['DEC'], unit=(u.hour, u.deg))
        gaia_stars = get_gaia_stars(coords, radius=args.radius)
        splusall, gaia_stars, gaia_matches, splus_matches = match_gaia_splus(
            args, gaia_stars, sfield=field)
        output = plot_gaia_splus(args, splusall, gaia_stars,
                                 gaia_matches, splus_matches, fieldname=field)
        return output


if __name__ == '__main__':
    args = get_args()
    output = main(args)
