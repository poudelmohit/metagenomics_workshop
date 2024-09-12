# Downloadig data first:


# Getting softwares ready:

    mkdir software && cd $_ # creating SOFTWARE directory:

    ## libreoffice:
    wget https://appimages.libreitalia.org/LibreOffice-still.basic-x86_64.AppImage
    chmod +x LibreOffice-still.*
    ./LibreOffice-still.*

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    conda activate metagenome
    conda install bioconda::fastqc -y
    conda install bioconda::kraken2 -y
    conda install bioconda::trimmomatic -y
    conda install bioconda::krona -y
    conda install bioconda::maxbin2 -y
    conda install bioconda::spades -y
    conda install bioconda::kraken-biom -y
    conda install bioconda::checkm-genome -y

##
    bash /home/mp067823/anaconda3/envs/metagenome/opt/krona/updateTaxonomy.sh
    cd software                        
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf taxdump.tar.gz 
    mkdir .taxonkit
    cp names.dmp nodes.dmp delnodes.dmp merged.dmp ./.taxonkit
    rm *dmp readme.txt taxdump.tar.gz gc.prt
    cd ..

#### installing R packages:
    R

    install.packages("ggplot2")

    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")

    BiocManager::install("phyloseq")


# downloading data:
    
    mkdir data && cd $_
    wget https://zenodo.org/records/7010950/files/dc_workshop.zip
    wget https://zenodo.org/records/7010950/files/MGRAST_MetaData_JP.xlsx

    unzip dc_workshop.zip
    ls dc_workshop/data






    

