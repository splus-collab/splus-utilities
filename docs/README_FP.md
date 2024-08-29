## Field properties VAC code generator


The VAC of field properties contains a series of parameters that may be useful for different science cases. The catalog includes depth per field per filter, median FWHM per field per filter, average reddening per field, and moon brightness. The parameters are estimated as follows:

- The depths are given for 4 SN ratios, 3, 5, 10, and 50 for the photometries auto, petro, and iso. To calculate the depth only objects with SExtractor photometric flags = 0, with the depth calculated as the median of the selected sample.

- The FWHM is also calculated for the magnitudes auto, petro and iso, photometric flags = 0, and only objects with FWHM in a given catalog bigger than 0 are acceptable. As per the depths, the value per field per filter is the result of the median of the selected sample.

- The reddening per field we use the package [SFDmap](https://github.com/kbarbary/sfdmap) to estimate the reddening per field. We estimate the reddening from the mean of the extinctions calculated for 10 positions across the field.

- The moon brightness is calculated for the moment of the observation of the first image used for the coadding of a given filter/field. We also provide the distance between the central coordinates of the field and the center of the Moon for the same point in time. This estimation is the best approximation possible given the complexity inherited from the observation strategy and the constraints imposed by the coadding process.

The code used to produce the VAC is part of the splus-utilities suite and is available in the [source directory](../src/).
