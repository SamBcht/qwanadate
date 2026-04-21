
<!-- README.md is generated from README.Rmd. Please edit that file -->

# qwanadate (PROTOTYPE) <img src="man/figures/logo.png" align="right" height="136" alt="" />

<!-- badges: start -->
<!-- badges: end -->

Quantitative Wood Anatomy dating : Shiny application to interactively
date tree-rings in QWA images

qwanadate is an interactive app designed to help absolute dating of wood
anatomical images analysed with the qwanamiz Python package.

The app is currently a prototype.

## Installation

The app is not yet entirely packaged. You can download the following R
script :

[qwanadate_app](https://github.com/SamBcht/qwanadate/blob/ec0b8fdb6593d237f9bc30337a98b05385b59d75/R/qwanadate_app.R)

Then open and run it in R or RStudio

You can install the development version of qwanadate from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("SamBcht/qwanadate")
```

## Example

``` r
library(qwanadate)
## basic example code
```

Main crossdating window of qwanadate :

<img src="man/figures/README-qwanadate-1.PNG" alt="" width="800px" />

## Required inputs

The app requires two types of input data:

### 1. Anatomical ring data

- A **base directory** containing one or more QWAnamiz output folders  

- Each folder must follow the pattern:

- Inside these folders, the app expects files named:

These files are automatically detected after entering the **Tree ID**
and selecting the base directory.

------------------------------------------------------------------------

### 2. Reference dendrochronological series

A ring-width file in one of the following formats:

- `.csv` file containing:

- `year` → calendar year  

- `TRW` → ring width values  

- `Tree.ID` → tree identifier (optional but recommended)

- `.rwl` file (standard dendro format)

The user must also specify the **unit scaling** of the TRW values
(e.g. mm, µm).

------------------------------------------------------------------------

### 3. Tree identifier

- A **Tree ID** must be provided to:
- match anatomical files  
- select the corresponding reference series

If no match is found in the dendro file, the app allows manual
selection.

------------------------------------------------------------------------

## Minimal workflow

1.  Enter **Tree ID**  
2.  Select **base directory** (QWAnamiz outputs)  
3.  Load **reference TRW file**  
4.  Click **“Search files”**  
5.  Run crossdating
