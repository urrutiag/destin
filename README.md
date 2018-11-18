## Destin

This is the Jiang lab scATAC-seq processing and cell clustering pipeline.

## Manuscript

Urrutia, Eugene, et al. "Destin: toolkit for single-cell analysis of chromatin accessibility." bioRxiv (2018): 461905. [link](https://www.biorxiv.org/content/early/2018/11/05/461905?rss=1)

## Questions & issues
  If you have any questions or problems when using destin, please feel free to open a new issue [here](https://github.com/urrutiag/destin/issues). You can also email the maintainers of the corresponding packages -- the contact information is below.
  
## Installation

The bioinformatic pipeline requires cloning the git repostory from github, where yourPathToDestinRepo is your path to the local cloned repository

Running the vignettes also requires cloning the git repository

```bash
cd yourPathToDestinRepo
git clone https://github.com/urrutiag/destin.git
```

install dependencies
```r
installed <- rownames(installed.packages())
pkgs = c("ChIPpeakAnno", "cluster", "ClusterR", "data.table", 
         "GenomicAlignments", "ggplot2", "gridExtra", "irlba",  "Matrix", 
         "parallel", "rtracklayer", "Rtsne")
pkgs <- setdiff(pkgs, installed)
if (length(pkgs))
  install.packages(pkgs, dep=c("Depends", "Imports"))
```

Running the R package requires either installing from the above git repostory locally
```{r} 
install.packages("yourPathToDestinRepo/package", repos = NULL, type = "source")
library(destin)
```

or downloading from github directly (note that this will not allow for the bioinformatics pipeline or the vignette):
```r
install.packages("devtools")
devtools::install_github("urrutiag/destin/package")
library(destin)
```

## Dependencies

- software: SRAtoolkit, cutadapt, BOWTIE2, 
            samtools, picard, MACS2, bedtools, awk,
            R, python

- R packages:  
 ChIPpeakAnno, 
 cluster,
 ClusterR,
 data.table, 
 GenomicAlignments,
 ggplot2,
 gridExtra,
 irlba, 
 Matrix, 
 parallel, 
 rtracklayer,
 Rtsne
 
## Overview 

### Bioinformatics Pipeline

Input: fastq files of entire experiment or individual fastq files by cell, 
       set of 2 for each of paired reads

Output: bam files by cell, peaks file

---

- download fastq

---

- separate fastq by cell (if combinatorial indexed)

---

- cut adapters 
- align 
- sam to bam
- sort
- Add read group and index
- mark duplicates
- remove mitochondrial, unmapped and chr Y
- adjust for Tn5 insertion
- alignment quality >= 30
- index

---

- call peaks (p < 0.01)
- filter blacklist

---

### Destin

Input: bam files by cell, peaks file

Output: cluster membership, differential accessibility

- create ranged summarized experiment from bam files and peaks file
- append experimental information if available
- annotate regions
- quality control on cells and regions
- determine number of clusters
- cluster cells by destin which optimizes hyperparameters via multinomial likelihood
- calculate differential accessibility

### GWAS association

Determine whether GWAS results are associated with increased chromatin accessibility in a particular cell type cluster.  We utilize 2 methods originally developed for scRNA-seq expression: ECWE and MAGMA.


## Example Workflows

-Bioinformatics and Clustering: Buenrostro mouse cells, Fluidigm microfluidic technology 
[html](https://rawgit.com/urrutiag/destin/master/package/vignettes/destinBuenrostroMouse.html)
[markdown](https://github.com/urrutiag/destin/blob/master/package/vignettes/destinBuenrostroMouse.Rmd)

-Clustering: Preissl P56 forebrain mouse cells, combinatorial barcode technology
[html](https://rawgit.com/urrutiag/destin/master/package/vignettes/destinPreisslP56.html)
[markdown](https://github.com/urrutiag/destin/blob/master/package/vignettes/destinPreisslP56.Rmd)

-GWAS cell-type specific association: Preissl P56 forebrain mouse cells 
[html](https://rawgit.com/urrutiag/destin/master/package/vignettes/GWAS.html)
[markdown](https://github.com/urrutiag/destin/blob/master/package/vignettes/GWAS.Rmd)

## Citations

de Leeuw, C. A., et al. (2015). Magma: generalized gene-set analysis of gwas data.
PLoS comput. biol., 11 (4), e1004219.

Skene,  N. G. et al. (2016).   Identification of vulnerable cell types in major brain
disorders  using  single  cell  transcriptomes  and  expression  weighted  cell  type
enrichment. Front. neurosci-switz,10, 16.

## Developers & Maintainers

* Gene Urrutia (gene dot urrutia at gmail dot com)
  <br>
  Department of Biostatistics, UNC-Chapel Hill

* [Yuchao Jiang](http://sph.unc.edu/adv_profile/yuchao-jiang-phd/) (yuchaoj at email dot unc dot edu)
  <br>
  Department of Biostatistics & Department of Genetics, UNC-Chapel Hill


