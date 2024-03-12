## Mask Regions Around Saturated Stars

This Python script is designed to prepare S-PLUS (Southern Photometric Local Universe Survey) data by applying masks to the S-PLUS catalogues. Here's a breakdown of what the script does:

## Function Definitions:

The script contains several function definitions:
```python
get_splusfootprint():
```
Retrieves the S-PLUS footprint from a CSV file.
```python
query_gsc():
```
Queries the GSC1.2 catalog using the Vizier service based on RA and DEC coordinates.
```python
get_stars():
```
Processes FITS image files to extract information about stars and S-PLUS objects.
```python
plot_stars():
```
Plots stars and S-PLUS objects on FITS images.
```python
make_masks():
```
Generates masks based on S-PLUS and GSC catalogs.
```python
check_distance_to_border():
```
Checks the distance of objects to the field borders.
```python
main():
```
The main function that orchestrates the entire process by parsing arguments, retrieving necessary data, processing it, and generating output.

## Execution

Finally, the script executes the main() function if it's run as the main program.

This script is designed to be executed from the command line, taking various parameters to customize its behavior. It processes S-PLUS data, applies masks based on GSC catalogs, and outputs masked catalogs. Additionally, it can plot the processed data for visualization.
