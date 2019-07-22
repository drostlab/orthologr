
```sh
conda create --name orthologr R perl gcc clang
conda activate orthologr
conda install -c bioconda orthofinder
conda install -c bioconda blast
conda install -c biobuilds clustalw 
conda install -c bioconda clustalo
conda install -c biobuilds t-coffee
conda install -c bioconda muscle
conda install -c bioconda mafft
conda install -c r r-xml
conda install -c bioconda bioconductor-rtracklayer
```

```r
# open R
install.packages("BiocManager")
BiocManager::install()
BiocManager::install("remotes")
# Install package dependencies
BiocManager::install("Biostrings")
# install orthologr from GitHub
BiocManager::install("HajkD/orthologr")
```

```sh
mkdir tools
cd tools
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/kaks-calculator/KaKs_Calculator1.2.tar.gz
gzip -d KaKs_Calculator1.2.tar.gz
tar -xf KaKs_Calculator1.2.tar
# install
cd KaKs_Calculator1.2/src
sudo make
sudo cp KaKs_Calculator /usr/local/bin/
```

- In the folder `KaKs_Calculator1.2/src` add to file `base.h` the line `#include<cstring>`.
- In the folder `KaKs_Calculator1.2/src` remove from file `makefile` the line `CFLAGS=` and all variables `$(CFLAGS)`.

```sh
conda env export > orthologr.yml
conda deactivate
```

When reading from `*.yml` file:

```sh
conda env create -f orthologr.yml --name orthologr
```



