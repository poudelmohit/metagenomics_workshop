cd $metagenomics_tutorial

### unzipping raw reads:

   gunzip data/dc_workshop/data/untrimmed_fastq/*.gz
   ls data/dc_workshop/data/untrimmed_fastq

### checking the last read of JP4D_R1.fastq:
    tail -n 4 data/dc_workshop/data/untrimmed_fastq/JP4D_R1.fastq | head -n


### activating conda env: 
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
    unzip $i;done
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


   
