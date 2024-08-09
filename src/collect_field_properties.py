#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Collect depth, reddening, FWHM, and sky brightness for each field and filter

import os
import argparse
from astropy.io import fits
import logging
import pandas as pd
import numpy as np
import sfdmap
from astropy.coordinates import SkyCoord


def parseargs():
    parser = argparse.ArgumentParser(
        description="Collect depth, reddening, FWHM, and sky brightness for each field and filter")
    parser.add_argument(
        '--workdir', default='./',
        help='Working directory. Default is current directory')
    parser.add_argument(
        "--fields", help="List of fields to process")
    parser.add_argument(
        "--data_dir", help="Directory containing data files")
    parser.add_argument(
        "--output_file", help="Output file")
    parser.add_argument(
        '--footprint', help='Footprint file')
    parser.add_argument(
        '--test', action='store_true', help='Test mode')

    return parser.parse_args()


def call_logger(logfile=None, loglevel=logging.INFO):

    logger = logging.getLogger(__name__)

    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] @%(module)s.%(funcName)s() %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")

    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(loglevel)
    if logfile is not None:
        file_handler = logging.FileHandler(logfile)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    return logger


class FieldProperties:
    def __init__(self, args, logger, field):
        self.logger = logger
        self.workdir = args.workdir
        self.data_dir = args.data_dir
        self.field = field
        self.footprint = None
        self.filters = {'U': 'u',
                        'F378': 'J0378',
                        'F395': 'J0395',
                        'F410': 'J0410',
                        'F430': 'J0430',
                        'G': 'g',
                        'F515': 'J0515',
                        'R': 'r',
                        'I': 'i',
                        'F660': 'J0660',
                        'F861': 'J0861',
                        'Z': 'z'}
        self.catalogues = None

    def main(self):
        self.catalogues = self.get_catalogues()
        depths = self.get_depth()
        fwhms = self.get_fwhm()
        reddening = self.get_reddening()
        sky_brightness = self.get_sky_brightness()

        return depths, fwhms, reddening, sky_brightness

    def get_catalogues(self):
        catalogues = {}
        for f in self.filters.keys():
            cat_file = os.path.join(
                self.data_dir, f'{self.field}_{f}_dual.fits')
            if not os.path.exists(cat_file):
                self.logger.error(f'Catalogue {cat_file} not found')
                catalogues[f] = None
            else:
                catalogues[f] = fits.open(cat_file)[1].data
        return catalogues

    def get_depth(self):
        depths = {}
        mags_type = ['auto', 'petro', 'iso']
        s2ns = [3, 5, 10, 50]
        for f in self.filters.keys():
            for m in mags_type:
                for sn in s2ns:
                    tag = f'depth_{self.filters[f]}_{m}_{repr(sn)}'
                    if self.catalogues[f] is None:
                        depths[tag] = None
                    else:
                        scat = self.catalogues[f]
                        mask = scat[f'SEX_FLAGS_{self.filters[f]}'] == 0
                        mask &= scat[f'{self.filters[f]}_{m}'] > 10
                        mask &= scat[f'{self.filters[f]}_{m}'] < 30
                        mask &= scat[f's2n_{self.filters[f]}_{m}'] >= sn
                        depths[tag] = np.median(
                            scat[f'{self.filters[f]}_{m}'][mask])

        return depths

    def get_fwhm(self):
        fwhms = {}
        mags_type = ['auto', 'petro', 'iso']
        for f in self.filters.keys():
            for m in mags_type:
                tag = f'FWHM_{self.filters[f]}'
                if self.catalogues[f] is None:
                    fwhms[tag] = None
                else:
                    scat = self.catalogues[f]
                    mask = scat[f'SEX_FLAGS_{self.filters[f]}'] == 0
                    mask &= scat[f'{self.filters[f]}_{m}'] > 10
                    mask &= scat[f'{self.filters[f]}_{m}'] < 30
                    mask &= scat[f'FWHM_{self.filters[f]}'] > 0
                    fwhms[tag] = np.median(
                        scat[f'FWHM_{self.filters[f]}'][mask] * 3600)

        return fwhms

    def get_reddening(self):
        # Initialize SFD map
        sfd = sfdmap.SFDMap(os.path.join(
            os.getenv('HOME'), 'Documents/sfddata-master/'))
        # get field coordinates
        _field = self.footprint[self.footprint['NAME'] == self.field]
        central_coords = SkyCoord(
            ra=_field['RA'], dec=_field['DEC'], unit=('hour', 'deg'))
        _ra = [central_coords.ra.deg]
        _dec = [central_coords.dec.deg]
        # additional positions around central coordinates
        deltas = [-.35, 0, .35]
        # make pairs with all possible combinations of deltas
        _ = [(_ra.append(_ra[0] + i), _dec.append(_dec[0] + j))
             for i in deltas for j in deltas]
        scoords = SkyCoord(ra=_ra, dec=_dec, unit=('deg', 'deg'))
        reddening = sfd.ebv(scoords).mean()

        return reddening

    def get_sky_brightness(self):
        pass


if __name__ == '__main__':
    args = parseargs()
    logger = call_logger()
    idr5_fields = pd.read_csv(args.fields)
    footprint = pd.read_csv(args.footprint)
    footprint['NAME'] = [n.replace('_', '-') for n in footprint['NAME']]
    if args.test:
        field = idr5_fields['NAME'][0]
        fprops = FieldProperties(args, logger, field)
        fprops.footprint = footprint
        depths, fwhms, reddening, sky_brightness = fprops.main()
