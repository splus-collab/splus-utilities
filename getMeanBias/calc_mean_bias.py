#!/bin/python3
"""
This script calculates the mean bias of the S-PLUS data against time.
Author: Fabio R Herpich CASU/IOA/UK
Date: 19/07/2023
email: fabio.herpich@ast.cam.ac.uk
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.io import fits
import argparse
import multiprocessing as mp
import os

# create a parser function to receive rawdir, ncores and workdir


def parser():
    parser = argparse.ArgumentParser(
        description='This script calculates the mean bias of the S-PLUS data against time.')
    parser.add_argument('-r', '--rawdir', type=str, help='raw data directory')
    parser.add_argument('-w', '--workdir', type=str, help='working directory')
    parser.add_argument('-n', '--ncores', type=int, help='number of cores')
    # save image if requested
    parser.add_argument('-s', '--save', type=bool,
                        help='save image', default=False)
    args = parser.parse_args()

    return args

# this function receives a bias image and calculates the mean value


def calc_mean_bias(bias):
    mean = np.mean(fits.getdata(bias, ext=1))
    return mean


# calculate the mean of a list using multiprocessing
def mean_over_list(bias_list, ncores):
    print('Creating pool...')
    pool = mp.Pool(processes=ncores)
    print('Starting calculations...')
    mean_list = pool.map(calc_mean_bias, bias_list)
    print('Closing pool...')
    pool.close()
    print('Joining pool...')
    pool.join()
    print(mean_list)
    return mean_list

# define main function


def main():
    # get arguments
    args = parser()
    rawdir = args.rawdir
    workdir = args.workdir
    ncores = args.ncores
    save = args.save

    # get bias list
    bias_list = glob.glob(os.path.join(rawdir, '*/bias*.fz'))
    list_of_dates = [bias.split('/')[-2] for bias in bias_list]

    # calculate mean bias
    print('Calculating mean bias...')
    mean_list = mean_over_list(bias_list, ncores)

    # plot mean bias against time
    plt.figure(figsize=(10, 8))
    plt.plot(list_of_dates, mean_list, 'o')
    plt.xlabel('Date')
    plt.ylabel('Mean Bias value')
    if save:
        print('Saving image...')
        plt.savefig(os.path.join(workdir, 'mean_bias.png'))
    else:
        plt.show()


# run main function
if __name__ == '__main__':
    main()
