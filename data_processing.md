## Creating a directories to work on the project:
    by moving to the directory for this project (in my case 'metagenomics_workshop'):

    # creating directories:
    mkdir -p data/raw_reads results tools

## downloading raw fastq files:
    
    Though the originial research has 40 total raw fastq reads, I am using just 6 for now
    For that, First I will need SRA ids of the samples I intend to work with,
    Secondly, I need to download the files using those ids and a tool (described below):
    
### 1.Saving the SRA ids for the samples in a text file:
    
    touch data/sra_ids.txt
    echo -e "ERS1949765\nERS1949772\nERS1949791\nERS1949788\nERS1949792\nERS1949793" >> data/sra_ids.txt
    cat data/sra_ids.txt

### 2.Downloading files using fastq-dump:

    # Installing SRAtoolkit:
    cd tools
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
    tar -vxzf sratoolkit.*.tar.gz
    chmod +x sratoolkit.*/bin/*
    cd ../data

    # downloading fastq files:
    cat sra_ids.txt | xargs -n 1 ../tools/sratoolkit.*/bin/fastq-dump --split-files
    mv *.fastq raw_reads

## To Remeber:
    # Controls are the samples ending with 58 and 59
    # Unenriched samples are ending with 67 and 68
    # Fertilized samples are ending with 69 and 70


## Checking the last read of a control sample:
    tail -n 4 raw_reads/ERR*58_?.fastq 


## Installing tools in dedicated conda environment: 

    conda env list
    conda create -n metagenome -y
    conda activate metagenome

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    conda install -y fastqc=0.11.9 trimmomatic=0.39 kraken2=2.1.2 krona=2.8.1 maxbin2=2.2.7 spades=3.15.2 kraken-biom=1.2.0 checkm-genome=1.2.1

    # checking installation:
    which run_MaxBin.pl # for maxbin
    which spades
    which ktImportText # for krona
    which spades
    which checkm
    which kraken-biom

    
## Quality control with fastqc:
    
    fastqc -h
    cd raw_reads
    ls -lh

    fastqc *.fastq

    mkdir -p ../fastqc_results
    mv *.zip ../fastqc_results
    mv *.html ../fastqc_results

    cd ../fastqc_results && ls
    
    # At this point, I am copying the html reports to results folder:
    mkdir ../../results/initial_fastqc_reports && cp *.html ../../results/initial_fastqc_reports
    

### Trimming low quality reads:
    cd ../raw_reads && ls
    trimmomatic # just checking the parameters and options

#### trimming all files using for loop:

    for filename in *_1.fastq;do   
        r1=$filename
        r2=$(echo "$filename" | sed 's/_1/_2/')
        trimmed_r1=$(echo "$r1" | sed 's/^/trimmed_'/)
        trimmed_r2=$(echo "$r2" | sed 's/^/trimmed_'/)
        failed_r1=$(echo "$r1" | sed 's/^/trim_failed_'/)
        failed_r2=$(echo "$r2" | sed 's/^/trim_failed_'/)
    
        trimmomatic PE $r1 $r2 $trimmed_r1 $failed_r1 $trimmed_r2 \
        $failed_r2 SLIDINGWINDOW:4:25 MINLEN:35 ILLUMINACLIP:../TruSeq3-PE.fa:2:40:15
        
    done

    # Adapter sequences are saved in data directory in file TrueSeq3-PE.fa
    # Sliding window of 4 is used that removes bases if phred score is below 25
    # Also, any trimmed reads with less than 35 bases left is discarded

#### moving trimmed files to different directory:
    
    # currently, all trimmed files are inside raw_reads directory, which doesnot make sense.
    mkdir ../trimmed_files && mv ./trim*.fastq ../trimmed_files/
    ls ../trimmed_files

### Rerunning fastqc to check quality of trimmed files:
    cd ../trimmed_files
    for file in trimmed_*.fastq;do
    fastqc $file
    done
    
    # moving fastqc reports to different folder:
    mkdir ../fastqc_after_trimming && mv *_fastqc.* $_
    
    # I am satisifed with the fastqc reports after trimming, just moving on:
    
    # Also, copying the fastqc reports to results folder:
    mkdir -p ../../results/post_trimming_fastqc_reports && cp ../fastqc_after_trimming/*.html $_
    
### Metagenomic assembly:

    mkdir ../assembled_files/
    
    # zipping fasq files, because apparently metaspades work only with compressed files:
    gzip trimmed_*_?.fastq 
    
    for r1 in trimmed_*_1.fastq.gz;do
        r2=$(echo $r1 | sed 's/_1/_2/')
        assembled=$(echo $r1 | sed 's/^[^_]*_//' | sed 's/_1.fastq.gz/_assembled/')
        echo $assembled
        metaspades -1 "$r1" -2 "$r2" -o ../assembled_files/$assembled
    done

metaspades.py -1 trimmed_ERR2143758_1.fastq -2 trimmed_ERR2143758_2.fastq -o ../assembled_files/assembly_temp

    # assembling both samples one by one:
    
    metaspades.py -1 trimmed_JC1A_R1.fastq.gz -2 trimmed_JC1A_R2.fastq.gz -o ../assembled_files/assembly_JC1A
    metaspades.py -1 trimmed_JP4D_R1.fastq.gz -2 trimmed_JP4D_R2.fastq.gz -o ../assembled_files/assembly_JP4D

    cd ../assembled_files && ls

    ls ./assembly_JC1A
    ls ./assembly_JP4D/

    # out of these many files, the contigs.fasta and scaffolds.fasta are the ones with assembly
    # file with extension *.gfa holds information to visualize the assembly using programs like [Bandage](https://github.com/rrwick/Bandage)

### Metagenomic binning:
    
#### for JC1A:
    cd assembly_JC1A
    mkdir MAXBIN/JC1A
    run_MaxBin.pl -thread 8 -contig contigs.fasta -reads ../../trimmed_files/trimmed_JC1A_R1* -reads2 ../../trimmed_files/trimmed_JC1A_R2* -out MAXBIN/JC1A
    # it didn't run because the medium of marker gene number <= 1. So, program stopped.

#### for JP4D: 

    cd assembly_JP4D
    ls ../../trimmed_files
    mkdir MAXBIN
    run_MaxBin.pl -thread 8 -contig contigs.fasta \
    -reads ../../trimmed_files/trimmed_JP4D_R1* -reads2 ../../trimmed_files/trimmed_JP4D_R2* -out MAXBIN/JP4D
    ls MAXBIN

##### checking binned summary and quality:

    cat MAXBIN/JP4D.summary
    
    mkdir CHECM
    checkm taxonomy_wf domain Bacteria -x fasta MAXBIN/ CHECKM/

    # to save in a tsv format:
    checkm qa CHECKM/Bacteria.ms CHECKM/ --file CHECKM/quality_JP4D.txv --tab_table -o 2
    



### taxonomic assignment (trying with kraken2):

    # will be using trimmed but unassembled fastq files for taxonomic assignment:

    cd $metagenomics_tutorial/data/dc_workshop/data/assembled_files && ls
    
    ls assembly_JP4D
    mkdir ../taxonomy_files
    

    #installing kraken db (make sure you have ~100 GB of space for this):
   
    mkdir ../kraken2_db && cd $_
    wget -c https://refdb.s3.climb.ac.uk/kraken2-microbial/hash.k2d 
    # -c means continue to allow resuming without requiring to start over, incase download gets intruppted.
    wget https://refdb.s3.climb.ac.uk/kraken2-microbial/opts.k2d
    wget https://refdb.s3.climb.ac.uk/kraken2-microbial/taxo.k2d
    
    # searching against the created db:

    cd ../trimmed_files/ && ls
    kraken2 --db ../kraken2_db --threads 4 --paired trimmed_JP4D_R1.fastq.gz trimmed_JP4D_R2.fastq.gz --output ../taxonomy_files/JP4D.kraken --report ../taxonomy_files/JP4D.report
    dmesg | grep -i kill
    
    ## It didn't work due to less RAM availablity (memory issue)
    ## Deleting the unused kraken database: 
    rm -r ../kraken2_db ../kraken_db
    

### taxonomic assignment with metaphlan (because Kraken2 db needs >100GB and is constantly failing):
    
    cd .. && ls
    mkdir -p with_metaphlan/trimmed_raw_files
    cp trimmed_files/trimmed_*.gz with_metaphlan/trimmed_raw_files
    cd with_metaphlan/trimmed_raw_files && ls

    gunzip *.gz
    
    cd .. && ls
    
#### working for JP4D sample first:

    metaphlan trimmed_raw_files/trimmed_JP4D_R1.fastq,trimmed_raw_files/trimmed_JP4D_R2.fastq --bowtie2out JP4D_metagenome.bowtie2.b2 --nproc 7 --input_type fastq -o JP4D_profiled_metagenome.txt

    cat JP4D_profiled_metagenome.txt

    # Converting the output into taxonomy-based-profile:
    sgb_to_gtdb_profile.py -i JP4D_profiled_metagenome.txt -o JP4D_gtdb.txt
    cat JP4D_gtdb.txt

#### similarly for JC1A samples:

    metaphlan trimmed_raw_files/trimmed_JC1A_R1.fastq,trimmed_raw_files/trimmed_JC1A_R2.fastq --bowtie2out JC1A_metagenome.bowtie2.b2 --nproc 7 --input_type fastq -o JC1A_profiled_metagenome.txt

    
    cat JC1A_profiled_metagenome.txt
    # No microbial species were detected.    


### Interpreting the results:

    Here, unclassified sequences are discarded as the way code is set up.
    And then, out of the classified ones, all are of Bacterial domain (d stands for domain here)
    Within bacteria, 2 phyla (p) are detected: Proteobacteria and bacteroidota

    Similarly, o stands for order, g stands for genus, and s stands for species.

### Visualizing with krona:

#### for JP4D samples:
    
    metaphlan2krona.py --profile JP4D_profiled_metagenome.txt -k JP4D.krona.txt
    ktImportText JP4D.krona.txt -o JP4D_krona.html

    # saving the plot for use:
    mkdir ../../results
    cp JP4D_krona.html ../../results

#### for JC1A samples (delete later as there was no detection of microbial samples in JC1A):
    
    metaphlan2krona.py --profile JC1A_profiled_metagenome.txt -k JC1A.krona.txt
    ktImportText JC1A.krona.txt -o JC1A_krona.html

    # saving the plot for use:
    mkdir ../../results
    cp JC1A_krona.html ../../results











