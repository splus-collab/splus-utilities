# Plot S-PLUS Footprint on the Sky

This script plots the S-PLUS footprint on the sky using the tiles_nc.csv file and the 2MASS all sky catalog.

## Usage

```
python3 plot_footprint.py --splusfoot tiles_nc.csv
```

## Output

The above example will output a PNG file called `splus_footprint.png`.

![Example Output Plot](splus_footprint.png)

## Notes

The 2MASS catalog is available at http://cdsarc.u-strasbg.fr/viz-bin/Cat?II/281

The S-PLUS footprint contain the column STATUS, which can be used to filter the tiles. The STATUS column can have the following values:
- 0: Not observed
- 1: Observed OK
- 2: Observed with marginal seeing
- 3: Observed incomplete
- 4: Observed with potential contamination by moonlight
- 5: Observed with potential low transparency
- 6: Observed but some filters have less than 3 exposures
- -10: Removed from the footprint
- -5: New tiles added to the footprint
- -2,-1: Tiles with priority status at the moment of file creation

## License

This code is distributed under the [S-PLUS Member Exclusive License Agreement](../LICENSE). Please refer to the `LICENSE` file in the repository for more details.

## Author

The script is authored by Fabio R Herpich and can be reached at fabio.herpich@ast.cam.ac.uk.
