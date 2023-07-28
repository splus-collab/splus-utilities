# calc_mean_bias.py

This script calculates the mean bias of the S-PLUS data against time.

## Usage

```
python3 calc_mean_bias.py -r <raw_data_directory> -w <working_directory> -n <number_of_cores> [-s]
```

* `-r` is the path to the raw data directory.
* `-w` is the path to the working directory.
* `-n` is the number of cores to use for the calculations.
* `-s` is a flag to save the plot.

## Example

```
python3 calc_mean_bias.py -r /data/splus/raw -w /data/splus/work -n 4
```

This will calculate the mean bias of the S-PLUS data in the `/data/splus/raw` directory using 4 cores. The plot will be saved in the `/data/splus/work` directory.

## Output

The script will output a plot of the mean bias against time. The plot will also be saved in the working directory.

## Author

Fabio R Herpich CASU/IoA Cambridge
