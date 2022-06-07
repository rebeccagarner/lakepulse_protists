# Overview of workflow for 18S rRNA gene ASV inference

## Download and format read files
#### 1.a Download raw fastq.gz files
#### 1.b Check md5 sums 
#### 2. Change raw fastq file names
    Rscript 01_rename_fastq.R

## Trim primers and sequencing adapters
#### 3. Create list of samples
    ls *_R1.fastq.gz | cut -f1 -d "_" > samples
#### 4.a Run Cutadapt
    sh 03_cutadapt.sh
#### 4.b Summarize the fraction of reads and bp retained in each sample
    paste <(grep "passing" lp2017-2019watercolumn18s_cutadapt_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" lp2017-2019watercolumn18s_cutadapt_stats.txt | cut -f3 -d "(" | tr -d ")") > lp2017-2019watercolumn18s_cutadapt_summary.txt

#### 5. Run DADA2
    Rscript 04_dada2.R