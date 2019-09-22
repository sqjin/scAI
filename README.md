# scAI: a single cell Aggregation and Integration method for analyzing single cell multi-omics data

- scAI is an unsupervised approach for integrative analysis of gene expression and chromatin accessibility or DNA methylation proflies measured in the same individual cells.
- scAI infers a set of biologically relevant factors, which enable various downstream analyses, including the identification of cell clusters, cluster-specific markers and regulatory relationships. 
- scAI provides an intuitive way to visualize features (i.e., genes and loci) alongside the cells in two dimensions.
- scAI aggegrates chromatin profiles of similar cells in an unsupervised and iterative manner, which opens up new avenues for analyzing extremely sparse, binary scATAC-seq data. 

Once the single cell multi-omics data are decomposed into multiple biologically relevant factors, the package provides functionality for further data exploration, analysis, and visualization. Users can:

- Visualize the latent biological patterns of the multi-omics data
- Visualize both genes and loci alongside cells onto the same two-dimensional space
- Identify cell clusters from the inferred joint cell loading matrix and cluster-specific markers
- Visualize clusters and gene expression in the low-dimensional space such as VscAI, t-SNE and UMAP
- Infer regulatory relationships between cluster-specific chromatin region and marker genes


## Packages
scAI has been implemented as both **R package** and **MATLAB package**. In each package, we provide example workflows that outline the key steps and unique features of scAI. The MATLAB package and examples are available [here](https://github.com/amsszlh/scAI).


## Examples and Walkthroughs

All the R markdown used to generate the walkthroughs can be found under the /examples directory. 

- Simulated single cell RNA-seq and ATAC-seq data [(Walkthrough)](https://htmlpreview.github.io/?https://github.com/sqjin/scAI/blob/master/examples/walkthrough_simulation.html):  This simulated data were generated based on bulk RNA-seq and DNase-seq profiles from the same sample using MOSim package. 
- Paired single cell RNA-seq and ATAC-seq data of A549 cells [(Walkthrough)](https://htmlpreview.github.io/?https://github.com/sqjin/scAI/blob/master/examples/walkthrough_A549dataset.html): This data describes lung adenocarcinoma-derived A549 cells after 0, 1 and 3 hours of 100 nM dexamethasone treatment. 


## Installation 

Before installing scAI R package, please install [RcppEigen](https://github.com/RcppCore/RcppEigen) and [rfunctions](https://github.com/jaredhuling/rfunctions).
```
if(!require(devtools)){ install.packages("devtools")}
install.packages("RcppEigen")
devtools::install_github("jaredhuling/rfunctions")
```
**Troubleshooting**: Installing RcppEigen and rfunctions on R>=3.5 requires Clang >= 6 and gfortran-6.1. For MacOS, it's recommended to follow guidance on the official R page [here](https://cloud.r-project.org/bin/macosx/tools/) OR the [post](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos-before-r-3.6.0/).  For Windows, please ensure that [Rtools](https://cran.r-project.org/bin/windows/Rtools/) is installed. 


Then scAI R package can be installed by the following two ways. 

### Install from Github using devtools
```
devtools::install_github("sqjin/scAI")
```

### Install from R source codes
Download source codes [here](https://github.com/sqjin/scAI/blob/master/scAI_0.1.1.tar.gz) and type (in R)
```
install.packages(path_to_file, type = 'source', rep = NULL)
```
The path_to_file would represent the full path and file name:
- On macOS/UNIX it will look like this: "~/Downloads/scAI_0.1.0.tar.gz".
- On Windows it will look something like this: "C:\Downloads\scAI_0.1.0.tar.gz".


### Additional installation steps (possibly)
- Install other dependencies

scAI provides functionality for further data exploration, analysis, and visualization. A couple of excellent packages need to be installed. 
```
library(devtools)
install_github('linxihui/NNLM')
install_github("yanwu2014/swne")
install_github("jokergoo/ComplexHeatmap")
```
- Install Leiden algorithm for identifying cell clusters: pip install leidenalg. Please check [here](https://github.com/vtraag/leidenalg) if there is any trouble.

- Install UMAP and FIt-SNE for faster dimension reduction in `reducedDims`

Using UMAP and FIt-SNE is recommended for computational efficiency when using `reducedDims` on very large datasets.

-- install UMAP Python package: pip install umap-learn. Please check [here](https://github.com/lmcinnes/umap) if there is any trouble. 

-- install FIt-SNE R package:  Installing and compiling the necessary software requires the use of FIt-SNE and FFTW. For detailed instructions of installation, please visit this [page](https://github.com/KlugerLab/FIt-SNE).


### Troubleshooting on the R Compiler Tools for Rcpp on macOS
If you get the error "clang: error: unsupported option '-fopenmp'" when installing R package, please consider the configuration in ~/.R/Makevars and see this [post](https://thecoatlessprofessor.com/programming/r-compiler-tools-for-rcpp-on-macos-before-r-3.6.0/) for detailed configuration. In addition, you may can also reinstall your R because -fopenmp option is usually added by R automatically if openmp is available.

If you are using macOS Mojave Version (10.14) and you might get the error "/usr/local/clang6/bin/../include/c++/v1/math.h:301:15: fatal error: 'math.h' file not found", please check the [post](https://github.com/RcppCore/Rcpp/issues/922). This error can be solved if running the following on the terminal:

```
sudo installer -pkg \
/Library/Developer/CommandLineTools/Packages/macOS_SDK_headers_for_macOS_10.14.pkg \
-target /
```

## Help
If you have any problems, comments or suggestions, please contact us at Suoqin Jin (suoqin.jin@uci.edu) or Lihua Zhang (lihuaz1@uci.edu).




