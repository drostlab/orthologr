---
title: "Performing Divergence Stratigraphy"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Performing Divergence Stratigraphy}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
---

```{r, echo = FALSE, message = FALSE}
options(width = 750)
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = FALSE)
```

The `orthologr` package allows users to perform __Divergence Stratigraphy__ for any query and subject organisms of interest.

__Divergence Stratigraphy__ is the process of quantifying the selection pressure (in terms of protein evolutionary rate) acting on orthologous genes between closely related species. The resulting sequence divergence map (short divergence map), stores the divergence stratum in the first column and the `query_id` of inferred orthologous genes in the second column ( Quint et al., 2012 _Nature_; [Drost et al., 2015 _Mol. Biol. Evol._](http://mbe.oxfordjournals.org/content/32/5/1221.full); [Drost et al., 2016 _Mol. Biol. Evol._](http://mbe.oxfordjournals.org/content/early/2016/02/23/molbev.msw039.abstract); [Introduction to myTAI](https://github.com/HajkD/myTAI) ).

The following Algorithm implemented in `divergence_stratigraphy()` defines __Divergence Stratigraphy__ as method (see [Drost et al., 2015](http://mbe.oxfordjournals.org/content/32/5/1221.full)):

1) Orthology Inference using BLAST best reciprocal hit ("RBH") based on blastp

2) Pairwise global amino acid alignments of orthologous genes using the [Needleman-Wunsch algorithm](http://www.sciencedirect.com/science/article/pii/0022283670900574)

3) Codon alignments of orthologous genes using [PAL2NAL](http://www.bork.embl.de/pal2nal/)

4) dNdS estimation using [Comeron's method (1995)](http://link.springer.com/article/10.1007/BF00173196)

5) Categorize estimated dNdS values into `divergence strata` (= deciles of all dNdS values)


In `orthologr` the [Needleman-Wunsch algorithm](http://www.sciencedirect.com/science/article/pii/0022283670900574), [PAL2NAL](http://www.bork.embl.de/pal2nal/) and [Comeron's method (1995)](http://link.springer.com/article/10.1007/BF00173196) 
are already included in the `orthologr` package and do not have to be installed separately. Nevertheless, users need to make sure __they have BLAST installed on their machine before using the `divergence_stratigraphy()`function__.  

__Note__: The following examples assume that the __BLAST__ program is installed and stored in the default execution path `usr/local/bin`. In case users do not have __BLAST__ installed yet or the following command in R produces a different output, please consult the [Installation Vignette](https://github.com/HajkD/orthologr/blob/master/vignettes/Install.Rmd#orthology-inference-tools) to corretly set up the __BLAST__ program to perform __Divergence Stratigraphy__.  

```{r,eval=FALSE}
system("blastp -version")
```

```
blastp: 2.2.30+
Package: blast 2.2.30, build Oct 27 2014 17:10:51
```



## Divergence Map Computations

In [Drost et al., 2015 _Mol. Biol. Evol._](http://mbe.oxfordjournals.org/content/32/5/1221.full) we define a `Divergence Map` as table storing the degree of selection pressure (= `divergence strata`) for each protein coding gene of a given query organism. In this case selection pressure was quantified by dNdS estimation (ratio of synonymous versus non-synonymous `codon -> amino acid` sequence substitution rates).
The resulting dNdS values for all protein coding genes of the query organism are then categorized into deciles (10%-quantiles) allowing users to compare the results obtained from `Phylostratigraphy` with results obtained form `Divergence Stratigraphy`. 

To perform __Divergence Stratigraphy__ using `orthologr` users need to retrieve the following input files:

* a CDS file covering all protein coding genes of the query organism of interest
* a CDS file covering all protein coding genes of the subject organism of interest

### Sequence Data Retrieval

In the following example, we will use _Arabidopsis thaliana_ as query organism
and _Arabidopsis lyrata_ as subject organism.

First, we need to download the CDS sequences for all protein coding genes of _A. thaliana_
and _A. lyrata_.

#### Option 1:

The CDS retrieval can be done using a `Terminal` or by manual downloading the files

* `Arabidopsis_thaliana.TAIR10.23.cds.all.fa.gz`
* `Arabidopsis_lyrata.v.1.0.23.cds.all.fa.gz`

```shell

# download CDS file of A. thaliana
curl ftp://ftp.ensemblgenomes.org/pub/
plants/release-23/fasta/arabidopsis_thaliana/
cds/Arabidopsis_thaliana.TAIR10.23.cds.all.fa.gz 
-o Arabidopsis_thaliana.TAIR10.23.cds.all.fa.gz

# unzip the fasta file
gunzip -d Arabidopsis_thaliana.TAIR10.23.cds.all.fa.gz

# download CDS file of A. lyrata

curl ftp://ftp.ensemblgenomes.org/pub/plants/
release-23/fasta/arabidopsis_lyrata/cds/
Arabidopsis_lyrata.v.1.0.23.cds.all.fa.gz 
-o Arabidopsis_lyrata.v.1.0.23.cds.all.fa.gz

# unzip the fasta file
gunzip -d Arabidopsis_lyrata.v.1.0.23.cds.all.fa.gz

```
When the download is finished you need to unzip the files.


#### Option 2:

We implemented the [biomartr](https://github.com/HajkD/biomartr) package to automate the process of performing biological data retrieval. The [Sequence Retrieval Vignette](https://github.com/HajkD/biomartr/blob/master/vignettes/Sequence_Retrieval.Rmd) stored in `biomartr` provides detailed use cases for the automation of biological sequence retrieval.

__Note__: Users need to make sure they have [biomartr](https://github.com/HajkD/biomartr#fast-installation-guide) installed before running any `biomartr` functions.

### Computation Time

__Please note that__ performing __Divergence Stratigraphy__ with two large genomes can take
(even on a multicore machine) some time -> __up to several hours__.
On a 4 core machine with 3.4 GHz i7 processors the computation time of generating a divergence map between _A. thaliana_ and _A. lyrata_ was __2.5-3 hours__.

The `comp_cores` argument implemented in the `divergence_stratigraphy()` function allows users to specify the number of cores they would like to use on their machine. The default value is `comp_cores = 1` which might take __10-12h__ to execute. So users need to make sure that they use all cores available on their machine to speed up the computation time.

## Running `divergence_stratigraphy()`

As mentioned earlier the `divergence_stratigraphy()` function is the main function to perform the __Divergence Stratigraphy__ algorithm.

In `divergence_stratigraphy()` the `query_file` and `subject_file` arguments take an
character string storing the path to the corresponding __fasta__ files containing the CDS
sequences of these organisms. Here the previously downloaded CDS sequence files of _A. thaliana_ (= `query_file`) and _A. lyrata_ (= `subject_file`) need to be specified. The `eval` is set to `1E-5` (default ; see Quint et al., 2012 _Nature_) and `BLAST best reciprocal hit` is used for orthology inference (see [Drost et al., 2015](http://mbe.oxfordjournals.org/content/32/5/1221.full)). In case ``orthologr` is running on a multicore machine, users can set the `comp_cores` argument to any number of cores
supported by their machine. The `clean_folders` argument indicates whether or not the internal folder structure should be deleted (cleaned) after processing is finished. In this case all output files generated by `divergence_stratigraphy` (stored in `tempdir()`) will be removed after the `Divergence Map` was returned. The `quiet` argument indicates whether or not a successful interface call should be printed out to the console (`quiet = FALSE`) or not (`quiet = TRUE`).


```{r,eval=FALSE}
library(orthologr)

# compute the divergence map of A. thaliana
Athaliana_DM <- divergence_stratigraphy(
                         query_file      = "path/to/Arabidopsis_thaliana.TAIR10.23.cds.all.fa",
                         subject_file    = "path/to/Arabidopsis_lyrata.v.1.0.23.cds.all.fa",
                         eval            = "1E-5", 
                         ortho_detection = "RBH",
                         comp_cores      = 1, 
                         quiet           = TRUE, 
                         clean_folders   = TRUE )

```

Before running `divergence_stratigraphy()` with two complete genomes, users can first run a test __Divergence Stratigraphy__ with 20 example genes that are stored in the `orthologr` package:

```{r,eval=FALSE}
library(orthologr)

# performing standard divergence stratigraphy
 divergence_stratigraphy(
      query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
      eval            = "1E-5", 
      ortho_detection = "RBH", 
      dnds.threshold  = 2,
      comp_cores      = 1, 
      quiet           = TRUE, 
      clean_folders   = TRUE)

```


```

   divergence_strata    query_id
1                 10 AT1G01010.1
2                  9 AT1G01020.1
3                  5 AT1G01030.1
4                  4 AT1G01040.1
5                  1 AT1G01050.1
6                  9 AT1G01060.3
7                  6 AT1G01070.1
8                  8 AT1G01080.1
9                  2 AT1G01090.1
10                 7 AT1G01110.2
11                 2 AT1G01120.1
12                 3 AT1G01140.3
13                10 AT1G01150.1
14                 8 AT1G01160.1
15                 1 AT1G01170.2
16                 6 AT1G01180.1
17                 7 AT1G01190.1
18                 4 AT1G01200.1
19                 5 AT1G01210.1
20                 3 AT1G01220.1

```

The resulting output is a `Divergence Map` of the 20 example genes.

To save corresponding `Divergenec Maps` to a hard drive users can
pass the resulting `divergence_stratigraphy()` output to a variable
and then use the `write.table()` function implemented in R to store the
`Divergenec Map` as `*.csv` file.

```{r,eval=FALSE}
Athaliana_DM <- divergence_stratigraphy( ... )

write.table(x         = Athaliana_DM, 
            file      = "Ath_Aly_DivergenceMap.csv", 
            sep       = ";", 
            col.names = TRUE, 
            row.names = FALSE, 
            quote     = FALSE)
```

This way `write.table()` will store the `Divergence Map` to the users current working directory (= `getwd()`).

### Specifying the arguments in `divergence_stratigraphy()`

Several argument combinations can be specified in `divergence_stratigraphy()` (see `Arguments` in `?divergence_stratigraphy`). This section introduces additional output options of `divergence_stratigraphy()`.


#### Example: `blast_path`

Sometimes the machine users are working on does not allow them to install __BLAST__ in the default execution path `usr/local/bin`. For this purpose the `blast_path` argument is implemented in `divergence_stratigraphy()`. This argument takes an character string specifying the `PATH` to the user's `blastp` execution file that is stored in a different place than `usr/local/bin`.

The following example shows a possible specification of `blast_path`.

```{r,eval=FALSE}
library(orthologr)

# performing standard divergence stratigraphy
 divergence_stratigraphy(
      query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
      blast_path      = "here/the/path/to/blastp", 
      eval            = "1E-5", 
      ortho_detection = "RBH", 
      dnds.threshold  = 2, 
      comp_cores      = 1, 
      quiet           = TRUE, 
      clean_folders   = TRUE)

```

```

   divergence_strata    query_id
1                 10 AT1G01010.1
2                  9 AT1G01020.1
3                  5 AT1G01030.1
4                  4 AT1G01040.1
5                  1 AT1G01050.1
6                  9 AT1G01060.3
7                  6 AT1G01070.1
8                  8 AT1G01080.1
9                  2 AT1G01090.1
10                 7 AT1G01110.2
11                 2 AT1G01120.1
12                 3 AT1G01140.3
13                10 AT1G01150.1
14                 8 AT1G01160.1
15                 1 AT1G01170.2
16                 6 AT1G01180.1
17                 7 AT1G01190.1
18                 4 AT1G01200.1
19                 5 AT1G01210.1
20                 3 AT1G01220.1

```

#### Example: `ds.values`

As defined earlier, a `Divergence Map` stores the `divergence strata` for protein coding genes of a `query` organism. However, `divergence strata` are based on `dNdS` values that were categorized into deciles. For this reason it is not possible to map a `divergence stratum` value to the exact initial `dNdS` value. So in case users are interested in the the exact `dNdS` value of protein coding genes, they can specify the `ds.values` argument in `divergence_stratigraphy()` allowing them to retrieve a `dNdS Map` instead of a `Divergence Map`. For this purpose users need to set `ds.values = FALSE`.


```{r,eval=FALSE}
library(orthologr)

# performing standard divergence stratigraphy
# but receive a dNdS Map instead of a Divergence Map
 divergence_stratigraphy(
      query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
      eval            = "1E-5", 
      ortho_detection = "RBH",
      ds.values       = FALSE, 
      dnds.threshold  = 2,
      comp_cores      = 1, 
      quiet           = TRUE, 
      clean_folders   = TRUE)

```


```
      dNdS    query_id
1  0.41950 AT1G01010.1
2  0.38790 AT1G01020.1
3  0.11850 AT1G01030.1
4  0.11560 AT1G01040.1
5  0.00000 AT1G01050.1
6  0.39670 AT1G01060.3
7  0.17280 AT1G01070.1
8  0.32170 AT1G01080.1
9  0.04174 AT1G01090.1
10 0.26620 AT1G01110.2
11 0.02317 AT1G01120.1
12 0.04324 AT1G01140.3
13 0.64120 AT1G01150.1
14 0.37310 AT1G01160.1
15 0.00000 AT1G01170.2
16 0.16830 AT1G01180.1
17 0.17730 AT1G01190.1
18 0.11370 AT1G01200.1
19 0.13420 AT1G01210.1
20 0.10230 AT1G01220.1

```

The corresponding output now stores `dNdS` values instead of `DS` values in the first column.

#### Example: `subject.id`

Although the `Divergence Map` [standard](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd#getting-started) is specified as storing `DS` values in the first column and `GeneIDs` in the second column, in some cases it is important to store the GeneIDs of orthologous genes in the `subject` organism. The `subject.id` argument implemented in `divergence_stratigraphy()` allows users to retrieve the GeneIDs of the orthologous genes of the `subject` organism. For this purpose users need to specify `subject.id = TRUE`.


```{r,eval=FALSE}

# receive a Divergence Map with DS | query GeneID | orthologous subject GeneID 
divergence_stratigraphy(
      query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
      eval            = "1E-5",
      ortho_detection = "RBH",
      comp_cores      = 1,
      quiet           = TRUE,
      clean_folders   = TRUE,
      subject.id      = TRUE)

```

```
   DS    query_id            subject_id
1  10 AT1G01010.1 333554|PACid:16033839
2   9 AT1G01020.1 470181|PACid:16064328
3   5 AT1G01030.1 470180|PACid:16054974
4   4 AT1G01040.1 333551|PACid:16057793
5   1 AT1G01050.1 909874|PACid:16064489
6   9 AT1G01060.3 470177|PACid:16043374
7   6 AT1G01070.1 918864|PACid:16052578
8   8 AT1G01080.1 909871|PACid:16053217
9   2 AT1G01090.1 470171|PACid:16052860
10  7 AT1G01110.2 333544|PACid:16034284
11  2 AT1G01120.1 918858|PACid:16049140
12  3 AT1G01140.3 470161|PACid:16036015
13 10 AT1G01150.1 918855|PACid:16037307
14  8 AT1G01160.1 918854|PACid:16044153
15  1 AT1G01170.2 311317|PACid:16052302
16  6 AT1G01180.1 909860|PACid:16056125
17  7 AT1G01190.1 311315|PACid:16059488
18  4 AT1G01200.1 470156|PACid:16041002
19  5 AT1G01210.1 311313|PACid:16057125
20  3 AT1G01220.1 470155|PACid:16047984
```

The resulting output now shows DS values, query GeneIDs, and orthologous subject GeneIDs.

A similar output can be generated for `dNdS` values instead of `DS` values by specifying `ds.values = FALSE` and `subject.id = TRUE`.

```{r,eval=FALSE}

# receive a dNdS Map with dNdS | query GeneID | orthologous subject GeneID 
divergence_stratigraphy(
      query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
      eval            = "1E-5",
      ortho_detection = "RBH",
      comp_cores      = 1,
      ds.values       = FALSE,
      quiet           = TRUE,
      clean_folders   = TRUE,
      subject.id      = TRUE)

```

```
      dNdS    query_id            subject_id
1  0.41950 AT1G01010.1 333554|PACid:16033839
2  0.38790 AT1G01020.1 470181|PACid:16064328
3  0.11850 AT1G01030.1 470180|PACid:16054974
4  0.11560 AT1G01040.1 333551|PACid:16057793
5  0.00000 AT1G01050.1 909874|PACid:16064489
6  0.39670 AT1G01060.3 470177|PACid:16043374
7  0.17280 AT1G01070.1 918864|PACid:16052578
8  0.32170 AT1G01080.1 909871|PACid:16053217
9  0.04174 AT1G01090.1 470171|PACid:16052860
10 0.26620 AT1G01110.2 333544|PACid:16034284
11 0.02317 AT1G01120.1 918858|PACid:16049140
12 0.04324 AT1G01140.3 470161|PACid:16036015
13 0.64120 AT1G01150.1 918855|PACid:16037307
14 0.37310 AT1G01160.1 918854|PACid:16044153
15 0.00000 AT1G01170.2 311317|PACid:16052302
16 0.16830 AT1G01180.1 909860|PACid:16056125
17 0.17730 AT1G01190.1 311315|PACid:16059488
18 0.11370 AT1G01200.1 470156|PACid:16041002
19 0.13420 AT1G01210.1 311313|PACid:16057125
20 0.10230 AT1G01220.1 470155|PACid:16047984
```

#### Example: `dnds.threshold`

`Divergence Strata` are obtained by categorizing dNdS values into deciles. For decilation the range of dNdS values is important. The `dnds.threshold` defines the upper level cut off of dNdS values. Since dNdS values are in the range [0, +Inf] a upper threshold needs to be specified.
The default value for `dnds.threshold` in `divergence_stratigraphy()` is `dnds.threshold = 2` due to the interpretation of dNdS values for predicting sequence evolution (dNdS < 1 -> negative selection; dNdS = 1 -> neutral selection; dNdS > 1 -> positive selection). Hence, all dNdS values `> 1` predict positive selection. In my experience of computing dNdS values between hundreds of pairwise species comparisons covering all evolutionary distances, dNdS values of orthologous genes rarely take values > 2. Nevertheless, in case you wish to extend or reduce the upper threshold for dNdS values, you can specify the `dnds.threshold` in `divergence_stratigraphy()`.


```{r,eval=FALSE}
# upper threshold for dNdS: dnds.threshold = 5
 divergence_stratigraphy(
      query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
      eval            = "1E-5", 
      ortho_detection = "RBH",
      ds.values       = FALSE, 
      dnds.threshold  = 5,
      comp_cores      = 1, 
      quiet           = TRUE, 
      clean_folders   = TRUE)

```

#### Example: `ortho_detection`

According to [Drost et al., 2015 _Mol. Biol. Evol._](http://mbe.oxfordjournals.org/content/32/5/1221.full) the __Divergence Stratigraphy__ algorithm performs __BLAST__ best reciprocal hit (RBH) as orthology inference method. Despite this convention, the `ortho_detection` argument allows users to perform orthology inference within the __Divergence Stratigraphy__ algorithm that is based on any orthology inference method implemented in `orthologr` (see `?orthologs` or [Orthology Inference Vignette](https://github.com/HajkD/orthologr/blob/master/vignettes/orthology_inference.Rmd) for details). For example in Quint et al., 2012 _Nature_ instead of using __BLAST__ best reciprocal hit, the method __BLAST__ best hit (BH) was used to perform orthology inference within the __Divergence Stratigraphy__ algorithm.


```{r,eval=FALSE}
# orthology inference method: ortho_detection = "BH"
 divergence_stratigraphy(
      query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
      subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
      eval            = "1E-5", 
      ortho_detection = "BH",
      ds.values       = TRUE, 
      dnds.threshold  = 2,
      comp_cores      = 1, 
      quiet           = TRUE, 
      clean_folders   = TRUE)

```

```
   DS    query_id
1  10 AT1G01010.1
2   9 AT1G01020.1
3   5 AT1G01030.1
4   4 AT1G01040.1
5   1 AT1G01050.1
6   9 AT1G01060.3
7   6 AT1G01070.1
8   8 AT1G01080.1
9   2 AT1G01090.1
10  7 AT1G01110.2
11  2 AT1G01120.1
12  3 AT1G01140.3
13 10 AT1G01150.1
14  8 AT1G01160.1
15  1 AT1G01170.2
16  6 AT1G01180.1
17  7 AT1G01190.1
18  4 AT1G01200.1
19  5 AT1G01210.1
20  3 AT1G01220.1
```


## Skip `Divergence Stratigraphy` and Download Already Published `Divergence Maps`

Users can find a detailed list of [published Phylostratigraphic Maps and Divergence Maps](https://github.com/HajkD/published_phylomaps) by following the link. This way
the computation time of 3-4 h on a local machine for 2 genome comparisions can be skipped.


## Combine a Divergence Map with Gene Expression Data

`Divergence Maps` can be used for a wide range of analyses. One example is to combine `Divergence Maps` with gene expression data to capture evolutionary signals in developmental transcriptomes (= `Phylotranscriptomics`; see [Drost et al., 2015 _Mol. Biol. Evol._](http://mbe.oxfordjournals.org/content/32/5/1221.full)). Performing phylotranscriptomic analyses based on an existing `Divergence Map` can easily be done by using the `myTAI` package. You can consult the [Introduction to the myTAI package Vignette](https://github.com/HajkD/myTAI/blob/master/vignettes/Introduction.Rmd)
for more details.




