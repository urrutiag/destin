## Destin

This is the Jiang lab scATAC-seq processing and cell clustering pipeline.

## Manuscript

## Questions & issues
  If you have any questions or problems when using destin, please feel free to open a new issue [here](https://github.com/urrutiag/destin/issues). You can also email the maintainers of the corresponding packages -- the contact information is below.
  
## Installation

```bash
git clone https://github.com/urrutiag/destin.git
```

```r
install.packages("devtools")
devtools::install_github("urrutiag/destin/destin")
```

## Dependencies

- software: SRAtoolkit, cutadapt, BOWTIE2, 
            samtools, picard, MACS2, bedtools, awk,
            R, python

- R packages:  ChIPpeakAnno, data.table, GenomicAlignments,
 ggplot2, irlba, Matrix, parallel, rtracklayer

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

## Example Workflows

-Buenrostro mouse cells, Fluidigm microfluidic technology 
[html](https://rawgit.com/urrutiag/destin/master/destin/vignettes/destin.html)
[markdown](https://github.com/urrutiag/destin/blob/master/destin/vignettes/destin.Rmd)

-Preissl P56 forebrain mouse cells, combinatorial barcode technology

## Citations

## Developers & Maintainers

* Gene Urrutia (urrutia at email dot unc dot edu)
  <br>
  Department of Biostatistics, UNC-Chapel Hill

* [Yuchao Jiang](http://sph.unc.edu/adv_profile/yuchao-jiang-phd/) (yuchaoj at email dot unc dot edu)
  <br>
  Department of Biostatistics & Department of Genetics, UNC-Chapel Hill


