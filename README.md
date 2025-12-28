# AqEquil

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13530217.svg)](https://doi.org/10.5281/zenodo.13530217)

Boyer, G., Robare, J., Park, N., Ely, T., Shock, E.L.

## About

AqEquil is a Python 3 package that enables users to rapidly perform aqueous speciation calculations of water chemistry data for multiple samples by interfacing with [geochemical speciation software EQ3/6](https://github.com/LLNL/EQ3_6) (Wolery 2013, [Wolery 1979](https://inis.iaea.org/collection/NCLCollectionStore/_Public/10/474/10474294.pdf)). AqEquil uses [a modified version of EQ3/6](https://github.com/39alpha/eq3_6/tree/main) created by the 39Alpha research team for easy local installation.

Water sample data in CSV format is automatically converted to a format readable by EQ3 and then speciated. Distributions of aqueous species, mineral saturation indices, oxidation reduction potentials, and more are data-mined and returned as Pandas tables and interactive Plotly visualizations.

Speciated fluids can be further reacted with minerals or other fluids in mass transfer calculations to produce tables and interactive diagrams of reaction paths and composition changes as a function of reaction progress.

## Requirements

AqEquil works on Linux, macOS, and Windows.

**Python Requirements:**
- Python >= 3.10
- pandas, numpy, matplotlib, plotly, and other dependencies (automatically installed)

**Additional Requirements:**
- R installation with CHNOSZ package version 1.4.3 (see instructions below)

**Note:** As of version 0.41.0, EQ3/6 executables are bundled with aqequil and no longer need to be installed separately. The package includes pre-compiled binaries for Linux, macOS, and Windows.

## Installation

### Installing CHNOSZ version 1.4.3

AqEquil requires an R installation with the CHNOSZ package. See [these instructions](https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/) for installing R with Anaconda.

Open an R session and install CHNOSZ version 1.4.3 package with:

```r
install.packages('http://cran.r-project.org/src/contrib/Archive/CHNOSZ/CHNOSZ_1.4.3.tar.gz', repos=NULL, type='source')
```

Once CHNOSZ is installed you can quit the R session.

Compatibility with CHNOSZ v.2.0.0 is forthcoming.

### Installing AqEquil

Install AqEquil using pip:

```bash
pip install aqequil
```

The bundled EQ3/6 executables will be automatically installed with the package. No additional configuration or environment variables are needed.

## Usage

See this [demo notebook](https://nbviewer.jupyter.org/github/worm-portal/WORM-Library/blob/master/3-Aqueous-Speciation/1-Introduction-to-Aq-Speciation/2-Intro-to-Multi-Aq-Speciation.ipynb) for usage examples.

## Bundled Software

This package includes pre-compiled binaries from [EQ3/6 v8.0a](https://github.com/39alpha/eq3_6), a software package for geochemical modeling developed by Thomas Wolery at Lawrence Livermore National Laboratory.

**EQ3/6 License:** BSD 3-Clause License
**Copyright:** (c) 1987, 1990-1993, 1995, 1997, 2002, 2013 The Regents of the University of California, Lawrence Livermore National Laboratory.

See `THIRD_PARTY_LICENSES.txt` for the full EQ3/6 license text.

**References:**
- Wolery, T. J., and USDOE. EQ3/6 A Software Package for Geochemical Modeling. Computer software. December 13, 2010. https://www.osti.gov//servlets/purl/1231666. doi:https://doi.org/10.11578/dc.20210416.44.
- Wolery, T. J. and R. L. Jarek. Software User's Manual EQ36, Version 8.0. U.S. Tech. Rep. 2003. Department of Energy, Office of Civilian Radioactive Waste Management, Office of Repository Development. 10813-UM-8.0-00.
