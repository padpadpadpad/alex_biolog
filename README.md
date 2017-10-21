## Cleaning, analysis and visualisation of sequencing and biolog data

This is the repository for the cleaning, analysis and visualisation of biolog plates and sequencing data, based on the data provided by Alex Vujakovic

### To run

- Download the folder using Git and GitHub or if you are not familiar with GitHub you can download the whole folder using the big green button on the top right of the screen (called "Clone or Download") and select "Download zip".
- In RStudio (I recommend using [RStudio](https://www.rstudio.com/products/rstudio/download/), open the Rstudio project file `alex.Rproj`. This will set all the paths relative to where this folder resides and will allow the analysis to be run across multiple machines.

#### Biolog analysis

- Open `package_install.R` and run the script. This will install all of the packages necessary to run the script.
- Open `analysis.R` and this will run the script for the cleaning, analysis and visualisation of the biolog plate data

#### Sequencing

- Open `package_install.R` and run the script to install all necessary packages for the analysis of sequence data
- The raw data files & filtered `.fastq` are not included in this online repository because they are too big. The raw files can be downloaded directly from the Liverpool CGS website or alternatively from Floh's iMac within the ESI.
- `raw_read_processing.R` inside `scripts/data_processing` takes a very long time to run and is probably best being run on the RStudio Server only.

Any problems please email d.padfield@exeter.ac.uk or post in the [Issues tab](https://github.com/padpadpadpad/alex_biolog/issues)