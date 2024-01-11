# 2021-09-17 herpich@usp.br

import pandas as pd
import numpy as np
import argparse
import os


def get_parser():
    parser = argparse.ArgumentParser(
        description='Get arguments to initialize the module.')
    parser.add_argument('--workdir', type=str,
                        default=os.path.dirname(os.path.realpath(__file__)),
                        help='Working directory. Default: script location')
    parser.add_argument('--filters', type=str,
                        default=os.path.join(os.path.dirname(
                            os.path.realpath(__file__)), 'extras/list_filters.csv'),
                        help=" ".join(['File containing list of filters and',
                                       'information regarding observations.',
                                       'Default is {workdir}/list_filters.csv']))
    parser.add_argument('--outfilename', type=str,
                        default=os.path.join(os.path.dirname(
                            os.path.realpath(__file__)), 'extras/block_list_filters.yaml'),
                        help="Output block file"
                        )

    return parser.parse_args()


def build_blockfile(args):
    df = pd.read_csv(args.filters)

    print('block size is', np.sum(df['EXPTIME'] *
          df['N']) + np.sum(df['N'])*30, 'seconds')

    ini_text = "pre-actions:\n"
    ini_text += "  - action: autofocus\n"
    ini_text += "    step: -1\n"
    ini_text += "\n"
    ini_text += "pos-actions:\n"
    ini_text += "  - action: autofocus\n"
    ini_text += "    step: 0\n"

    side = ['east', 'west']

    blockname = 'myScripts/block_'
    k = 0
    for i in range(df['FILTER'].size):
        for j in range(df['N'][i]):
            ini_text += "  - action: expose\n"
            ini_text += "    filter: %s\n" % df['FILTER'][i]
            ini_text += "    frames: 1\n"
            ini_text += "    exptime: %i\n" % np.int(df['EXPTIME'][i])
            ini_text += "    imageType: OBJECT\n"
            ini_text += "    compress_format: fits_rice\n"
            ini_text += '    objectName: "{name}"\n'
            ini_text += '    filename: "{pid}-$DATE-$TIME"\n'
            if k < np.sum(df['N']) - 1:
                ini_text += "  - action: point\n"
                ini_text += "    offset:\n"
                ini_text += "      %s: 10\n" % side[k % 2]
            k += 1
        blockname += df['FILTER'][i] + 'x' + \
            repr(df['N'][i]) + 'x' + repr(df['EXPTIME'][i])
        if i < df['FILTER'].size - 1:
            blockname += '-'

    blockname += '.yaml'

    block = open(blockname, 'w')
    block.write(ini_text)
    block.close()

    print('block file', blockname, 'is ready')
    return


if __name__ == '__main__':
    args = get_parser()
