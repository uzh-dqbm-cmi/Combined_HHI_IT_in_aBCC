# Combined_HHI_IT_in_aBCC
Code used to process the data and generate the figures for the paper.

The R markdown files contain the code for processing the count matrices (downloadable from https://doi.org/10.5281/zenodo.17814765) and for producing the figures of the manuscript.
`PBMC_figures.Rmd` has the code used to process the PBMC data and plot Figure 3 and Suppl. Figures 2 and 3. `Spatial_figures.Rmd` has the code for processing the spatial data and to plot Figure 4 and Suppl. Figures 4 and 5 of the manuscript. Also note, that an external dataset is used to provide reference for the cell-typing. Access to this dataset can be requested here: https://ega-archive.org/studies/EGAS00001005891 . `mesoscale_analysis.R` was used to generate the boxplots of the mesoscale data. The underlying data is uploaded to this repository as `mesoscale_absolute_values.csv`.

The clinical figures were produced using the `clinical.R` script and tables with the tables.Rmd notebook. The corresponding data have been uploaded to this repository.

The code was jointly written by Zsolt Balázs and Xiaocheng Yang. The code used for the analysis of the Mesoscale data were written by Patrick Turko and edited by Zsolt Balázs.
