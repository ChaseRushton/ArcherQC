# ArcherQC

A comprehensive QC analysis tool for Archer NGS data, providing automated quality control processing and reporting.

## Overview

ArcherQC is an R-based tool that performs quality control analysis on Archer Next Generation Sequencing (NGS) data. It processes molecular barcode (molbar) data, applies various QC thresholds, and generates detailed reports and visualizations.

## Features

- Processes molbar API data and historical records
- Applies multiple QC thresholds and validations
- Generates detailed reports and visualizations
- Handles historical data comparison
- Performs primer-specific analysis

## Requirements

- R version >= 4.0.0
- Minimum 8GB RAM

### R Dependencies

- tidyverse (data manipulation)
- RCurl (API requests)
- jsonlite (JSON processing)

## Installation

```R
# Install required R packages
install.packages(c("tidyverse", "RCurl", "jsonlite"))
```

## Usage

```R
Rscript ArcherQC.R
```

### Important Notes

- Ensure the working directory is empty except for molbar files when rerunning
- Verify all input files exist in specified paths
- Check file permissions for output directories
- Confirm historic data file format matches expected structure

## Development

- Created: 05/06/2024
- Last Updated: 07/18/2024
- Developer: Chase Rushton
- Contact: chase.rushton@pennmedicine.upenn.edu

## License



## Contributing

For internal use only. Please contact the developer for any modifications or improvements.
