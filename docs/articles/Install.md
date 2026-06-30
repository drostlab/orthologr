# Installing Prerequisite Tools

## Installing Prerequisite Tools

Most `orthologr` functions are interface functions that pass data to
common bioinformatics tools, internally call the corresponding tool, and
read their output as R object. For this purpose, when using interface
functions in `orthologr` users need to install the underlying
bioinformatics tools to obtain accurate results.

The following sections provide step by step instructions or guidance on
installing all bioinformatics tools for which R interface functions are
implemented in `orthologr`.

**Some tools are not trivial to install, so please read the
corresponding sections carefully and execute test cases that are
presented in each section.**

## Programming Languages

The following bioinformatics tools you are going to install are based on
the these programming languages:

- [**R**](https://cran.r-project.org) \>= 3.1.1

- [**C++11**](https://isocpp.org/about)

- [**Perl**](https://www.perl.org) \>= 5.12

Please make sure these programming languages are installed and
executable on the machines you are going to run `orthologr` on.

## Pairwise Sequence Alignment Tools

The `orthologr` package provides interfaces to the pairwise alignment
tools, `BLAST` and `DIAMOND v2`. We recommend the use of `DIAMOND v2` as
it saves time whilst being as sensitive as `BLAST`.

### Install `BLAST`

[**BLAST**](https://www.ncbi.nlm.nih.gov/guide/howto/run-blast-local/%5D)
(= Basic Local Alignment Search Tool) finds regions of similarity
between biological sequences and is also used as underlying paradigm of
most orthology inference methods.

1.  Go to
    <https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/> and
    download the system specific BLAST program.

2.  Install BLAST :

- On a Windows machine (see [installation manual:
  Windows](https://www.ncbi.nlm.nih.gov/books/NBK52637/)) -\> Please
  carefully read the `Environment Variables` section of the
  `installation manual: Windows` and make sure the execution `PATH`
  variable is set correctly.
- On a Unix machine (see [installation manual:
  Unix](https://www.ncbi.nlm.nih.gov/books/NBK52640/)) -\> Please
  carefully read the `Configuration` section of the
  `installation manual: Unix` and make sure the execution `PATH`
  variable is set correctly to `usr/local/bin`.

Example on Linux systems for BLAST+ version 2.2.31

For example for Linux systems open the `Terminal` application and run
(Thanks to Alexander Gabel):

    # download BLAST+ version 2.2.31
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz

    # extract the compiled version of BLAST
    tar zxvpf ncbi-blast+2.2.31+-x64-linux.tar.gz

    # copy BLAST files to `usr/local/bin`
    cp ncbi-blast-2.2.31+/bin/* usr/local/bin

Alternatively users can set the system call to the BLAST programs by
specifying the `PATH` variable (this is useful, because it allows an
easier update of BLAST versions instead of deleting all BLAST programs
from `usr/local/bin`):

    # open vim text editor
    vi .bash_profile

    # type 'Shift' then I to edit the file .bash_profile
    # and specify the export PATH
    export PATH=${PATH}:/path/to/downloaded/blast/folder/ncbi-blast-2.2.31+/bin

    # type 'ESC' then ':' then 'w' then 'q' to save and quit the .bash_profile file 

    # log out from your server with
    exit

    # log in again and type
    blastp -version

Now users should see the BLAST command line options.

#### Some tips

Based on our personal experience the installation of BLAST works best
when copy/pasting the BLAST executables to the path `usr/local/bin`. In
detail you can run the following steps to copy/paste the BLAST
executables to `usr/local/bin` (on Unix systems). However, updating
BLAST will then need to manually delete all previous BLAST programs from
`usr/local/bin` :

Open the Terminal application on your system and type:

    open /usr/local/bin 

Next, copy/paste the `blastp`, `makeblastdb`, etc files (BLAST
executables) from your BLAST folder to `/usr/local/bin`. To do so you
will need to enter the system password to allow the copy process.

Additional tips for macOS

While one can install BLAST+ via `brewer` and other package managers, we
recommend to install BLAST+ manually. Select `aarch64-macosx` for Apple
Silicon (M1/M2/M3 etc.) and `x64-macosx` for Intel-based Macs.

If you cannot run BLAST+, e.g. `blastp`, due to the message “Apple could
not verify “blastp” is free of malware that may harm your Mac or
compromise your privacy”, consider navigating to the directory where
`blastp` is located using the terminal and running a clear the
quarantine flag `xattr -r -d com.apple.quarantine .`. There may be other
solutions.

In addition, you may want to add BLAST+ to the .Renviron file to make
sure that BLAST+ is available in R sessions. To do so, open the R and
type:

``` r

library(usethis)
# open the .Renviron file for editing
usethis::edit_r_environ()
```

In the .Renviron file add the following line (make sure to change the
path to the BLAST+ bin folder):

    # e.g., for BLAST+ version 2.17.0+ installed in the Downloads folder of the user "user"
    PATH="/Users/User/Downloads/ncbi-blast-2.17.0+/bin:${PATH}"

After installing the BLAST program you can open an R session and type
the following command to check whether or not BLAST can be executed from
R.

``` r

# test whether blastp is correctly installed on your machine
system("blastp -version")
```

    blastp: 2.2.31+
    Package: blast 2.2.31, build Oct 27 2014 17:10:51

**You should see this output if BLAST was installed correctly.**

In case you find the following output:

    sh: blastp: command not found

You should return to `step 2)` and install BLAST so that it can be
executed from the default execution `PATH`.

These interface functions to BLAST+ are implemented in `orthologr`:

- [`blast()`](https://drostlab.github.io/orthologr/reference/blast.md) :
  Interface function to BLAST+
- [`blast_best()`](https://drostlab.github.io/orthologr/reference/blast_best.md)
  : Perform a BLAST+ best hit search
- [`blast_rec()`](https://drostlab.github.io/orthologr/reference/blast_rec.md)
  : Perform a BLAST+ reciprocal best hit (RBH) search
- [`set_blast()`](https://drostlab.github.io/orthologr/reference/set_blast.md)
  : Preparing the parameters and databases for subsequent BLAST+
  searches

### Install `DIAMOND2`

[**DIAMOND2**](https://github.com/bbuchfink/diamond) (= Double Index
alignment of Next-generation sequencing data) finds, like `BLAST`,
regions of similarity between biological sequences. Unlike `BLAST` it is
much much faster (up to 10 000X faster in the default `fast` mode and
over 80X faster in the `ultra-sensitive` mode, which is as sensitive as
`BLAST`). Thus, `DIAMOND2` facilitates *even faster* orthology
inference.

1.  Go to the download site in the [`DIAMOND2`
    wiki](https://github.com/bbuchfink/diamond/wiki/2.-Installation) and
    follow the instructions for installation. `DIAMOND2` is supported on
    Linux, macOS and Windows.

2.  Check the installation of `DIAMOND2` by running the command

&nbsp;

    diamond --version

Check the `Additional tips for macOS` section in the `BLAST`
installation instructions for macOS users, because similar issues may
arise when installing `DIAMOND2` on macOS.

3.  After installing the `DIAMOND2` program you can open an R session
    and type the following command to check whether or not `DIAMOND2`
    can be executed from R.

``` r

# test whether diamond is correctly installed on your machine
system("diamond --version")
```

    diamond version 2.1.8

**You should see this output if `DIAMOND2` was installed correctly.**

In case you find the following output:

    sh: diamond: command not found

You should return to `step 1)` and install `DIAMOND2` so that it can be
executed from the default execution `PATH`.

These interface functions to `DIAMOND2` are implemented in `orthologr`,
akin to the interface functions to `BLAST+`:

- [`diamond()`](https://drostlab.github.io/orthologr/reference/diamond.md)
  : Interface function to DIAMOND2
- [`diamond_best()`](https://drostlab.github.io/orthologr/reference/diamond_best.md)
  : Perform a diamond best hit search
- [`diamond_rec()`](https://drostlab.github.io/orthologr/reference/diamond_rec.md)
  : Perform a diamond reciprocal best hit (RBH) search
- [`set_diamond()`](https://drostlab.github.io/orthologr/reference/set_diamond.md)
  : Preparing the parameters and databases for subsequent diamond
  searches

Furthermore, the following functions use `DIAMOND2` by default, though
the use of BLAST can be specified through the parameter
`aligner = "blast"`:

- [`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md) :
  Compute dNdS values for two organisms
- [`divergence_stratigraphy()`](https://drostlab.github.io/orthologr/reference/divergence_stratigraphy.md)
  : Perform ‘Divergence Stratigraphy’

## Multiple Sequence Alignment Tools

The `orthologr` package also provides interfaces to the following
Multiple Alignment Tools. Nevertheless, non of them have to be installed
if the corresponding interface functions are not used.

### Install `ClustalW2`

To install `ClustalW2` please download via a package manager
(e.g. brewer, conda) since the ClustalW homepage is no longer available.
However, if you are somehow able to obtain the `clustalw2` program that
matches your operating system, you can install it by following the
instructions below.

After downloading and unpacking the `clustalw2` program, please go to
the clustalw-2.1 folder and open a `Terminal` application to type (in
this example for Mac OS X):

    # copy clustalw2 files to `usr/local/bin`
    cp clustalw2 usr/local/bin

### Install `T-Coffee`

To install `T-Coffee` please go to the
[T-Coffee](https://tcoffee.readthedocs.io/en/latest/index.html) homepage
and download the corresponding [T-Coffee
program](https://tcoffee.readthedocs.io/en/latest/index.html) matching
your operating system.

### Install `MUSCLE`

- [**MUSCLE**](https://www.drive5.com/muscle/) : Fast and accurate
  multiple alignment tool of nucleic acid and protein sequences

### Install `ClustalO`

1.  Download the [argtable](https://argtable.sourceforge.io/) program.

2.  Unzip the file.

3.  Run within the argtable folder:

&nbsp;

    ./configure

    make

    make check

    sudo make install

4.  Download ClustalO via package manager (e.g. brewer, conda).

### Install `MAFFT`

- [**MAFFT**](https://mafft.cbrc.jp/alignment/software/) : A tool for
  multiple sequence alignment and phylogeny

In `orthologr` the function
[`multi_aln()`](https://drostlab.github.io/orthologr/reference/multi_aln.md)
provides interfaces to all of these multiple alignment tools as well as
an pairwise alignment interface to the
[Biostrings](https://www.bioconductor.org/packages/release/bioc/html/Biostrings.html)
package performing a [Needleman-Wunsch
algorithm](https://www.sciencedirect.com/science/article/pii/0022283670900574).

## Codon Alignment Tools

The codon alignment tool `Pal2Nal` is already integrated in the
`orthologr` package and doesn’t need to be installed.

- [**Pal2Nal**](https://bio.tools/pal2nal/)

You don’t need to worry about downloading and installing **PAL2NAL**, it
is already included in the `orthologr` package. The corresponding
function
[`codon_aln()`](https://drostlab.github.io/orthologr/reference/codon_aln.md)
takes a protein alignment and the corresponding coding sequences and
returns a codon alignment by calling **Pal2Nal** from inside of the
`orthologr` package.

## dNdS Estimation Methods

dNdS estimation is a method to quantify the selection pressure acting on
a specific protein sequence determined by pairwise comparisons of amino
acid substitutions between two protein sequences and their corresponding
codon alignments. Different models have been proposed to estimate this
ratio quantifying selection pressure on proteins. The `orthologr`
package includes the most common dNdS estimation methods.

Starting with an codon alignment returned by
[`codon_aln()`](https://drostlab.github.io/orthologr/reference/codon_aln.md)
the function
[`dNdS()`](https://drostlab.github.io/orthologr/reference/dNdS.md)
computes the the dN, dS, and dNdS values of pairs of proteins.

Based on implementations provided by `gestimator`, `ape`, and
[KaKs_Calculator](https://code.google.com/archive/p/kaks-calculator),
the following dNdS Estimation Methods are available in `orthologr`:

- [Li](https://link.springer.com/article/10.1007/BF02407308) : Li’s
  method (1993) -\> provided by the `ape package`

- [Comeron](https://link.springer.com/article/10.1007/BF00173196) :
  Comeron’s method (1995)

- [NG](https://academic.oup.com/mbe/content/3/5/418.short) : Nei, M. and
  Gojobori, T. (1986)

- [LWL](https://academic.oup.com/mbe/article/2/2/150/1220392) : Li,
  W.H., et al. (1985)

- [MLWL](https://academic.oup.com/mbe/article/21/12/2290/1071055)
  (Modified LWL), MLPB (Modified LPB): Tzeng, Y.H., et al. (2004)

- [YN](https://academic.oup.com/mbe/article/17/1/32/975527) : Yang, Z.
  and Nielsen, R. (2000)

- [MYN](https://link.springer.com/article/10.1186/1471-2148-6-44)
  (Modified YN): Zhang, Z., et al. (2006)

For this purpose you need to have **KaKs_Calculator** installed on your
system and executable from your default `PATH`, e,g, `/usr/local/bin/`.

### Install KaKs_Calculator (For Linux/Unix/OS)

Please go to the [KaKs_Calculator
homepage](https://code.google.com/archive/p/kaks-calculator/downloads)
and download KaKs_Calculator.

E.g.

``` shell
# download KaKs_Calculator
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/kaks-calculator/KaKs_Calculator1.2.tar.gz
# unzip
gzip -d KaKs_Calculator1.2.tar.gz
tar -xf KaKs_Calculator1.2.tar
# install
cd KaKs_Calculator1.2/src
sudo make
sudo cp KaKs_Calculator /usr/local/bin/
```

Now you should be able to run KaKs_Calculator via `KaKs_Calculator -h`
in your bash or as `system("KaKs_Calculator -h")` in R.

The most recent version `KaKs_Calculator2.0` can be found
[here](https://sourceforge.net/projects/kakscalculator2/).
