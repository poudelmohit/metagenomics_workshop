    cd $metagenomics_tutorial

### unzipping raw reads:

   gunzip data/dc_workshop/data/untrimmed_fastq/*.gz

   ls data/dc_workshop/data/untrimmed_fastq

### checking the last read of JP4D_R1.fastq:

    tail -n 4 data/dc_workshop/data/untrimmed_fastq/JP4D_R1.fastq | head -n


### activating conda environment: 

    conda env list
    conda activate metagenome

### quality control with fastqc:
    
    fastqc -h
    cd data/dc_workshop/data/untrimmed_fastq
    ls -lh

    fastqc *.fastq

    mkdir -p ../fastqc_results
    mv *.zip ../fastqc_results
    mv *.html ../fastqc_results

    cd ../fastqc_results && ls
    
#### unizpping the fastqc results:

    for i in ls *.zip;do
        unzip $i
    done

    ls -F
    less JC1A_R1_fastqc/summary.txt

#### summarizing all summary into a docs:

    mkdir -p ../../docs
    cat */summary.txt > ../../docs/fastqc_summaries.txt
    cat ../../docs/fastqc_summaries.txt
    grep FAIL ../../docs/fastqc_summaries.txt


### Trimming low quality reads:

    trimmomatic # just checking the parameters and options
    cd ../untrimmed_fastq # navigating to untrimmed files directory to trim unwanted reads
    ls
    gzip *.fastq # zipping before running trimmomatic


#### trimming all files using for loop:

    for filename in *_R1.fastq.gz;do   
        r1=$filename
        r2=$(echo "$filename" | sed 's/_R1/_R2/')
        trimmed_r1=$(echo "$r1" | sed 's/^/trimmed_'/)
        trimmed_r2=$(echo "$r2" | sed 's/^/trimmed_'/)
        failed_r1=$(echo "$r1" | sed 's/^/trim_failed_'/)
        failed_r2=$(echo "$r2" | sed 's/^/trim_failed_'/)
    
        trimmomatic PE $r1 $r2 $trimmed_r1 $failed_r1 $trimmed_r2 \
        $failed_r2 SLIDINGWINDOW:4:20 MINLEN:35 ILLUMINACLIP:TruSeq3-PE.fa:2:40:15
    done

    # Adapter sequences are saved in current directory in file TrueSeq3-PE.fa
    # Sliding window of 4 is used that removes bases if phred score is below 20
    # Also, any trimmed reads with less than 35 bases left is discarded

#### moving trimmed files to different directory:

    # currently, all trimmed files are inside untrimmed_fastq directory, which doesnot make sense.
    mkdir -p ../trimmed_files && mv ./trim*.gz ../trimmed_files/
    cd ../trimmed_files && ls

### Metagenomic assembly:

    mkdir ../assembled_files
    
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
    



### taxonomic assignment:

    # will be using trimmed but unassembled fastq files for taxonomic assignment:
    cd ../../trimmed_files
    
    mkdir ../taxonomy_files
    
    free -h
    kraken2 --help

    #installing kraken db (make sure you have ~100 GB of space for this):
    kraken2-build --standard --threads 8 --db kraken_db

    # searching against the created db:
    kraken2 --db kraken_db --threads 8 --paired trimmed_JP4D_R1.fastq.gz trimmed_JP4D_R2.fastq.gz --output ../taxonomy_files/JP4D.kraken --report ../taxonomy_files/JP4D.report
    ls
    echo $KRAKEN_DEFAULT_DB



   
