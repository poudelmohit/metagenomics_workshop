## Creating a directory to work on the project:
    by moving to the directory for this project (in my case 'metagenomics_workshop'):

    # creating directories:
    mkdir -p data/raw_reads results tools

## downloading raw fastq files:
    
    Though the originial research has 40 total raw fastq reads, I am using just 8 (randomly) for now
    For that, First I will need SRA ids of the samples I intend to work with,
    Secondly, I need to download the files using those ids and a tool (described below):
    
### 1.Saving the SRA ids for the samples in a text file:
    
    touch data/sra_ids.txt
    
    echo -e "ERR2143758\nERR2143759\nERR2143768\nERR2143767\nERR2143784\nERR2143785\nERR2143794\nERR2143795" > data/sra_ids.txt
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
    # Fertilized samples are ending with 84, 85, 94 and 95

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
    which checkm
    which kraken-biom

## Additional Setup:
    
    bash ~/anaconda3/envs/metagenome/opt/krona/updateTaxonomy.sh 
    cd ~/extra/metagenomics_workshop/tools/ && ls                          
    wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar -xzf taxdump.tar.gz
    mkdir ~/.taxonkit && cp names.dmp nodes.dmp delnodes.dmp merged.dmp $_
    rm *dmp readme.txt taxdump.tar.gz gc.prt

## Quality assesment with fastqc:
    
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
    
## Trimming low quality reads:

    cd ../raw_reads && ls
    trimmomatic # just checking the parameters and options

    echo -e ">PrefixPE/1\nTACACTCTTTCCCTACACGACGCTCTTCCGATCT\n>PrefixPE/2\nGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" > ../TruSeq3-PE.fa
    # I obtained this file from the original data repository (zenodo) itself

#### trimming all files using for loop:

    for filename in *_1.fastq;do   
        r1=$filename
        r2=$(echo "$filename" | sed 's/_1/_2/')
        trimmed_r1=$(echo "$r1" | sed 's/^/trimmed_'/)
        trimmed_r2=$(echo "$r2" | sed 's/^/trimmed_'/)
        failed_r1=$(echo "$r1" | sed 's/^/trim_failed_'/)
        failed_r2=$(echo "$r2" | sed 's/^/trim_failed_'/)
    
        trimmomatic PE $r1 $r2 $trimmed_r1 $failed_r1 $trimmed_r2 $failed_r2 SLIDINGWINDOW:4:25 MINLEN:35 ILLUMINACLIP:../TruSeq3-PE.fa:2:40:15
    done
    
    # Adapter sequences (to trim) are saved in data directory in file TrueSeq3-PE.fa
    # Sliding window of 4 is used that removes bases if phred score is below 25
    # Also, any trimmed reads with less than 35 bases left is discarded

#### moving trimmed files to different directory:
    
    # currently, all trimmed files are inside raw_reads directory, which doesnot make sense.
    mkdir ../trimmed_files && mv ./trim*.fastq ../trimmed_files/
    ls ../trimmed_files

## Rerunning fastqc to check quality of trimmed files:
    cd ../trimmed_files
    for file in trimmed_*.fastq;do
        fastqc $file
    done
    
    # moving fastqc reports to different folder:
    mkdir ../fastqc_after_trimming && mv *_fastqc.* $_
    
    # I am satisifed with the fastqc reports after trimming, just moving on:
    
    # Also, copying the fastqc reports to results folder:
    mkdir -p ../../results/post_trimming_fastqc_reports && cp ../fastqc_after_trimming/*.html $_
    
## Metagenomic assembly:

    mkdir ../assembled_files/
    cd ../trimmed_files && ls

### zipping fastq files, because apparently metaspades work only with compressed files:
    gzip trimmed_*_?.fastq
    

### for loop for assembly

    for r1 in trimmed_*_1.fastq.gz;do
        r2=$(echo $r1 | sed 's/_1/_2/')
        assembled=$(echo $r1 | sed 's/^[^_]*_//' | sed 's/_1.fastq.gz/_assembled/')
        echo $assembled
        metaspades -1 "$r1" -2 "$r2" -o ../assembled_files/$assembled
    done

    cd ../assembled_files && ls

    ls ERR*68_assembled
    
    # out of these many files, the contigs.fasta and scaffolds.fasta are the ones with assembly
    # file with extension *.gfa holds information to visualize the assembly using programs like [Bandage](https://github.com/rrwick/Bandage)

## Metagenomic binning:

    # Binning of all samples using a for loop:
    # for that, I need contigs.fasta (from assembly) and trimmed reads:
    
    gunzip ../trimmed_files/*.gz
    
    mkdir ../binned_files
    
    for file in *_assembled;do
       
        sample_id=$(echo "$file" | sed 's/_assembled//')
        contigs="$file/contigs.fasta"
        r1="../trimmed_files/trimmed_${sample_id}_1.fastq"
        r2="../trimmed_files/trimmed_${sample_id}_2.fastq"
        out_file="../binned_files/${sample_id}"
       
        mkdir -p $out_file
       
        run_MaxBin.pl -thread 8 -contig $contigs -reads $r1 -reads2 $r2 -out $out_file >> ../binned_files/binning.log 2>&1
    
        echo "Processed sample: $sample_id, output saved to $out_file" >> ../binned_files/binning.log
    done

    ls -lh ../binned_files/*.summary

    ## At this point, I see summary of only 6 samples out of initial 8. 
    ## It is because, for the remaining 2 samples, the number of marker genes is less than 1, read more at ../binned_files/binning.log

    # lets move the binned fasta files into respective directory by creating subdirectories: 
    cd ../binned_files

    for file in ERR*;do
        dir="${file%%.*}"
        mkdir -p "$dir/"
        mv "$file" "$dir/"
    done

### copying the summary files into results folder:
    mkdir -p ../../results/summary_binning/
    find . -type f -iname "*.summary" -exec cp {} ../../results/summary_binning/ \;

    ls -lh */*.fasta | wc -l

### Checking the quality of the assembled genomes

#### only 6 samples (that passed binning) has .fasta output. So:

    passed_folders="ERR2143759 ERR2143767 ERR2143784 ERR2143785 ERR2143794 ERR2143795"

    for folder in $passed_folders;do
        mkdir -p "$folder"/checkm_outputs
        checkm taxonomy_wf domain Bacteria -x fasta "$folder" "$folder"/checkm_outputs/
    done

    # here, due to insufficinet memory, I have specified the marker at the domain level, specific only for bacteria also mentioned that binned files are in fasta format in folder: binned_files
    # the output will be saved in checkm_outputs


# ####################################################

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
    
    conda deactivate
    conda create -n metaphlan -y
    conda activate metaphlan
    conda install metaphlan -y

    mkdir ../../tools/tx_db
    metaphlan --install --bowtie2db ../../tools/tx_db

    cd .. && ls
    mkdir -p tx_assignment/trimmed_raw_files
    cp trimmed_files/trimmed_*.fastq tx_assignment/trimmed_raw_files/
    cd tx_assignment/trimmed_raw_files && ls
    
    cd .. && ls
    
#### for loop for taxonomic assignment:

    cd trimmed_raw_files

    for file in trimmed_*_1.fastq; do
        sample_id=$(basename "$file" _1.fastq | sed 's/trimmed_//')
        reverse_read=$(echo "$file" | sed 's/_1/_2/')
        
        # echo "$sample_id"
        # echo "$file"
        # echo "$reverse_read"
    
        metaphlan "$file","$reverse_read" --bowtie2db ../../../tools/tx_db/ --bowtie2out ../"$sample_id".bowtie2.bz2 --nproc 7 --input_type fastq -o ../"$sample_id"_profiled_metagenome.txt
    done

    cd ..

    cat ERR*795_*.txt

### Converting the output into taxonomy-based-profile:

    mkdir tx_profile

    for file in *_metagenome.txt;do
        output_file=$(echo "$file" | sed s'/_profiled_metagenome/_gtdb/')
        echo $output_file
        echo "$file"
        sgb_to_gtdb_profile.py -i "$file" -o tx_profile/"$output_file"
    done
 
### Interpreting the results:

    Here, unclassified sequences are discarded as the way code is set up.
    And then, out of the classified ones, all are of Bacterial domain (d stands for domain here)
    Within bacteria, 2 phyla (p) are detected: Proteobacteria and bacteroidota

    Similarly, o stands for order, g stands for genus, and s stands for species.

### Visualizing with krona:

    mkdir krona_files/

    for file in *_profiled_metagenome.txt;do
        output=$(echo "$file" | sed s'/profiled_metagenome/krona/') 
        metaphlan2krona.py --profile "$file" -k ./krona_files/"$output"
    done

    ls krona_files
    cat krona_files/ERR2143759_krona.txt

### Converting to graphical format (html)

    conda activate metagenome
    
    cd krona_files
    
    for file in *.txt;do
        output=$(echo "$file" | sed 's/.txt/.html/')
        ktImportText "$file" -o "$output"
    done

## saving the files in results folder:
    
    mkdir ../../../results/abundance_graphs && cp *.html $_
    
    mkdir ../../results
    cp JP4D_krona.html ../../results

# Next Step ?

### alpha and beta Diversity Calculation !!
