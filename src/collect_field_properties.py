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
from astropy.coordinates import SkyCoord, get_body, EarthLocation, AltAz
from astropy.time import Time
import multiprocessing as mp


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
        "--output_file", default='fields_properties.csv', help="Output file")
    parser.add_argument(
        '--footprint', help='Footprint file')
    parser.add_argument(
        '--istest', action='store_true', help='Test mode')
    parser.add_argument(
        '--obstimefile', default=None, help='Observation time file')
    parser.add_argument(
        '--nprocesses', default=1, type=int, help='Number of processes to use')
    parser.add_argument(
        '--loglevel', default=logging.INFO, help='Log level')
    parser.add_argument(
        '--logfile', default='collect_field_properties.log', help='Log file')

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
        self.obstimefile = args.obstimefile

    def main(self):
        logger.info(f'Processing field {self.field}')
        self.catalogues = self.get_catalogues()
        depths = self.get_depth()
        fwhms = self.get_fwhm()
        reddening = self.get_reddening()
        sky_brightness = self.get_sky_brightness()

        # build a Dataframe with the results where each field is a row
        df = pd.DataFrame()
        df['Field'] = [self.field]
        df_depths = pd.DataFrame([depths])
        df = pd.concat([df, df_depths], axis=1)
        df_fwhm = pd.DataFrame([fwhms])
        df = pd.concat([df, df_fwhm], axis=1)
        df['reddening'] = reddening
        df_skybrightness = pd.DataFrame([sky_brightness])
        df = pd.concat([df, df_skybrightness], axis=1)

        return df

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
        if self.obstimefile is None:
            return None
        else:
            times4field = self.obstimefile[self.obstimefile['Field']
                                           == self.field]
            targcoords = SkyCoord(
                ra=self.footprint[self.footprint['NAME'] == self.field]['RA'],
                dec=self.footprint[self.footprint['NAME']
                                   == self.field]['DEC'],
                unit=('hour', 'deg'))
            sky_brightness = {}
            for f in self.filters.keys():
                filterobstime = Time(times4field[f'MJD_{f}'], format='mjd')
                sun_position = get_body(
                    'sun', filterobstime, location=EarthLocation.of_site('ctio'))
                moon_position = get_body(
                    'moon', filterobstime, location=EarthLocation.of_site('ctio'))
                elongation = sun_position.separation(moon_position)
                moon_phase = np.arctan2(
                    sun_position.distance * np.sin(elongation),
                    moon_position.distance -
                    sun_position.distance * np.cos(elongation)
                )
                moon_illumination = (1 + np.cos(moon_phase))/2.0
                sky_brightness[f'skybr_{self.filters[f]}'] = moon_illumination.value[0]
                targaltaz = targcoords.transform_to(
                    AltAz(obstime=filterobstime, location=EarthLocation.of_site('ctio')))
                sky_brightness[f'targalt_{self.filters[f]}'] = targaltaz.alt.deg[0]

                moon_separation = targcoords.separation(moon_position)
                sky_brightness[f'moonsep_{self.filters[f]}'] = moon_separation.deg[0]

            return sky_brightness


if __name__ == '__main__':
    args = parseargs()
    logger = call_logger()
    idr5_fields = pd.read_csv(args.fields)
    footprint = pd.read_csv(args.footprint)
    footprint['NAME'] = [n.replace('_', '-') for n in footprint['NAME']]
    obstimetab = pd.read_csv(os.path.join(args.workdir, args.obstimefile))
    if args.istest:
        field = idr5_fields['NAME'][0]
        fprops = FieldProperties(args, logger, field)
        fprops.footprint = footprint
        fprops.obstimefile = obstimetab
        df = fprops.main()
    else:
        results = []
        with mp.Pool(args.nprocesses) as pool:
            for field in idr5_fields['NAME']:
                fprops = FieldProperties(args, logger, field)
                fprops.footprint = footprint
                fprops.obstimefile = obstimetab
                results.append(pool.apply_async(fprops.main))
            pool.close()
            pool.join()

        df = pd.concat([r.get() for r in results], ignore_index=True)

    output_file = os.path.join(args.workdir, args.output_file)
    logger.info(f'Saving results to {output_file}')
    df.to_csv(output_file, index=False)

# End of file
