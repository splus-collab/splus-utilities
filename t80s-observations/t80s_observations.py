#!/bin/env python

import os
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import sfdmap
import gc
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description='Plot T80S observations since starting operations')
    parser.add_argument('--workdir', type=str, default=os.getcwd(),
                        help='Working directory')
    parser.add_argument('--sfddata', type=str, default=os.getcwd(),
                        help='Path to sfddata')
    parser.add_argument('--savefig', action='store_true',
                        help='Save figures to disk')
    parser.add_argument('--all', action='store_true',
                        help='Plot all data at once')
    return parser.parse_args()


def prepare_dateframe(args):
    workdir = args.workdir
    df = pd.read_csv(os.path.join(
        workdir, 't80s_observations-20160801-20240201.txt'))
    df = df.applymap(lambda x: x.replace(' ', '') if isinstance(x, str) else x)
    df.columns = df.columns.str.strip()

    dates = [date[-15:] for date in df['Name']]

    pids = np.array([pid.split('-')[0] if pid[:9] !=
                     'SPLUS-GAL' else 'SPLUS-GAL' for pid in df['Name']])
    datetimes = [(datetime.datetime.strptime(
        date, '%Y%m%d-%H%M%S').replace(tzinfo=datetime.timezone.utc) - datetime.timedelta(hours=12)).strftime('%Y-%m-%d') for date in dates]
    is_valid = np.int_(['0' if a == 'None' else a for a in df['is_valid']])

    date_col = []
    pid_col = []
    object_col = []
    ra_col = []
    dec_col = []
    is_valid_col = []
    numobs_col = []
    for date in np.unique(datetimes):
        mask = np.array(datetimes) == date
        for obj in np.unique(df['Object'][mask]):
            sumisvalid = sum(is_valid[mask][df['Object'][mask] == obj] > 0)
            date_col += [date]
            pid_col += [pids[mask][df['Object'][mask] == obj][0]]
            object_col += [obj]
            ra_col += [df['RA'][mask][df['Object'][mask] == obj].iloc[0]]
            dec_col += [df['DEC'][mask][df['Object'][mask] == obj].iloc[0]]
            is_valid_col += [sumisvalid]
            numobs_col += [sum(mask & (df['Object'] == obj))]

    print('Saving cleaned dataframe to t80s_observations_clean.csv')
    new_df = pd.DataFrame({'date': date_col, 'pid': pid_col, 'object': object_col,
                           'ra': ra_col, 'dec': dec_col, 'is_valid': is_valid_col, 'numobs': numobs_col})
    new_df.to_csv(os.path.join(
        workdir, 't80s_observations_clean.csv'), index=False)


def get_ebv(args, res=200, save=True):
    ra = np.linspace(-180, 180, res)
    dec = np.linspace(-90, 90, res)
    mesra, mesdec = np.meshgrid(ra, dec)
    coords = SkyCoord(ra=mesra.flatten(), dec=mesdec.flatten(),
                      unit='degree', frame='icrs')
    m = sfdmap.SFDMap(args.sfddata)
    mesebv = m.ebv(coords, unit='degree')
    if save:
        print('Saving ebv to ebv_%i.csv' % res)
        df = pd.DataFrame(
            {'ra': mesra.flatten(), 'dec': mesdec.flatten(), 'ebv': mesebv})
        df.to_csv(os.path.join(args.workdir, 'ebv_%i.csv' % res), index=False)
    return coords, mesebv


def plot_data(args):
    df = pd.read_csv(os.path.join(args.workdir, 't80s_observations_clean.csv'))
    c = SkyCoord(ra=df['ra'], dec=df['dec'],
                 unit=(u.deg, u.deg), frame='icrs')
    ra_rad = c.ra.wrap_at(180 * u.deg).radian
    dec_rad = c.dec.radian

    res = 500
    if os.path.exists(os.path.join(args.workdir, 'ebv_%i.csv' % res)):
        ebvdf = pd.read_csv(os.path.join(args.workdir, 'ebv_%i.csv' % res))
        coords = SkyCoord(ra=ebvdf['ra'], dec=ebvdf['dec'],
                          unit=(u.deg, u.deg), frame='icrs')
        ebv = ebvdf['ebv']
    else:
        coords, ebv = get_ebv(res=res, save=True)

    if not os.path.isdir(os.path.join(args.workdir, 'plots')):
        os.mkdir(os.path.join(args.workdir, 'plots'))

    pidcolours = {'SPLUS': 'b', 'SPLUS-GAL': 'g', 'STRIPE': 'r', 'HYDRA': 'c',
                  'MC': 'm', 'SHORTS': 'y', 'CN': 'orange', 'EXTMONI': 'limegreen',
                  'M2CAL': 'pink', 'SN': 'purple', 'other': 'k'}
    print('Assigning colours to PIDs')
    for i, pid in enumerate(df['pid']):
        found = False
        for key in pidcolours.keys():
            if key in pid:
                found = True
                break
        if not found:
            df['pid'][i] = 'other'
        if 'MC' in pid:
            df['pid'][i] = 'MC'
        if 'SHORTS' in pid:
            df['pid'][i] = 'SHORTS'
        if 'HYDRA' in pid:
            df['pid'][i] = 'HYDRA'
    pids = df['pid']
    print('Plotting data')

    for date in np.unique(df['date']):
        figname = os.path.join(
            args.workdir, 'plots/t80s_observations_' + date + '.png')
        if os.path.exists(figname):
            print('Plot for date', date, 'already exists. Skipping...')
            continue
        else:
            mask = df['date'] <= date
            plt.figure(figsize=(16, 8.4))
            ax = plt.subplot(111, projection="aitoff")
            plt.grid(True)
            ax.scatter(coords.ra.wrap_at(180 * u.deg).radian,
                       coords.dec.radian,
                       c=ebv,
                       cmap='gray_r', marker='H', s=5, vmin=-.5, vmax=10)
            for pid in np.unique(pids[mask]):
                if pid == 'SPLUS-GAL':
                    colour = pidcolours['SPLUS-GAL']
                elif pid == 'SPLUS':
                    colour = pidcolours['SPLUS']
                elif 'STRIPE' in pid:
                    colour = pidcolours['STRIPE']
                elif 'HYDRA' in pid:
                    colour = pidcolours['HYDRA']
                elif 'MC' in pid:
                    colour = pidcolours['MC']
                elif 'SHORTS' in pid:
                    colour = pidcolours['SHORTS']
                elif 'CN' in pid:
                    colour = pidcolours['CN']
                elif 'EXTMONI' in pid:
                    colour = pidcolours['EXTMONI']
                elif 'M2CAL' in pid:
                    colour = pidcolours['M2CAL']
                elif 'SN' in pid:
                    colour = pidcolours['SN']
                else:
                    colour = pidcolours['other']

                alpha = 1
                maskpid = mask & (df['pid'] == pid)
                ax.scatter(ra_rad[maskpid],
                           dec_rad[maskpid],
                           label=pid,
                           marker='o',
                           color=colour,
                           s=5,
                           alpha=alpha)
            plt.legend(loc='upper left', bbox_to_anchor=(-.165, 1.1), ncol=1)
            plt.title('T80S observations since starting operations', y=1.1)
            textstr = '\n'.join(('Date: %s' % date,
                                 'N observations: %10.0f' % sum(df['numobs'][mask])))
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.8, 1.1, textstr, transform=ax.transAxes,
                    fontsize=14, verticalalignment='top', bbox=props)
            if args.savefig:
                print('Saving figure', figname)
                plt.savefig(figname, dpi=150, bbox_inches='tight')
                plt.close()
            else:
                plt.show()
        gc.collect()


def plot_data_all(args):
    df = pd.read_csv(os.path.join(args.workdir, 't80s_observations_clean.csv'))
    c = SkyCoord(ra=df['ra'], dec=df['dec'],
                 unit=(u.deg, u.deg), frame='icrs')
    ra_rad = c.ra.wrap_at(180 * u.deg).radian
    dec_rad = c.dec.radian

    if os.path.exists(os.path.join(args.workdir, 'ebv.csv')):
        ebvdf = pd.read_csv(os.path.join(args.workdir, 'ebv.csv'))
        coords = SkyCoord(ra=ebvdf['ra'], dec=ebvdf['dec'],
                          unit=(u.deg, u.deg), frame='icrs')
        ebv = ebvdf['ebv']
    else:
        coords, ebv = get_ebv(save=True)

    if not os.path.isdir('plots'):
        os.mkdir('plots')

    for date in np.unique(df['date']):
        figname = os.path.join(
            args.workdir, 'plots/t80s_observations_' + date + '.png')
        if os.path.exists(figname):
            print('Plot for date', date, 'already exists. Skipping...')
            continue
        else:
            mask = df['date'] <= date
            plt.figure(figsize=(16, 8.4))
            ax = plt.subplot(111, projection="aitoff")
            plt.grid(True)
            ax.scatter(coords.ra.wrap_at(180 * u.deg).radian,
                       coords.dec.radian,
                       c=ebv,
                       cmap='gray_r', marker='H', s=5, vmin=-.5, vmax=10)
            c_bar = ax.scatter(ra_rad[mask],
                               dec_rad[mask],
                               c=df['is_valid'][mask],
                               label=date,
                               cmap='plasma',
                               marker='s',
                               vmin=0,
                               vmax=max(df['is_valid']))
            cb = plt.colorbar(c_bar, ax=ax)
            cb.set_label('N invalid')
            plt.legend(loc='upper left', bbox_to_anchor=(0, 1.2), ncol=1)
            plt.title('T80S observations since starting operations', y=1.1)
            if args.savefig:
                print('Saving figure', figname)
                plt.savefig(figname, dpi=150, bbox_inches='tight')
                plt.close()
            else:
                plt.show()
                break


def main(args):
    if not os.path.exists('t80s_observations_clean.csv'):
        prepare_dateframe(args)
    plot_data(args)


if __name__ == '__main__':
    args = parse_args()
    main(args)
