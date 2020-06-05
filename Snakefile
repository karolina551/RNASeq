SAMPLES, = glob_wildcards("data/{sample}_R1_001.fastq.gz")
READS=['R1_001', 'R2_001']

rule all:
    input:
         "fastqc_raw/multiqc_report.html",
         expand("clean/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=READS),
         "fastqc_clean/multiqc_report.html",
          directory('STAR'),
         expand("mapped/{sample}.bam", sample=SAMPLES),
         expand("mapped/{sample}.bam.bai", sample=SAMPLES),
         expand("counts/{sample}.txt", sample=SAMPLES),
         "counts/matrix.txt"
         
         
rule fastqc_raw:
    input: expand("data/{sample}_{read}.fastq.gz", sample=SAMPLES, read=READS)
    output: 
        "fastqc_raw/{sample}_{read}_fastqc.html",
        "fastqc_raw/{sample}_{read}_fastqc.zip"
    shell:'''
    fastqc -o fastqc_raw {input}
    '''

rule multiqc_raw:
    input: expand("fastqc_raw/{sample}_{read}_fastqc.html", sample=SAMPLES, read=READS)
    output: "fastqc_raw/multiqc_report.html"
    shell:'''
    multiqc -o fastqc_raw fastqc_raw
    '''

rule bbduk:
    input:"data/{sample}_{read}.fastq.gz"
    output: "clean/{sample}_{read}_trimmed.fastq.gz"
    shell:'''
    bbduk.sh  in={input} out={output} ref=resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=r trimq=10
    '''

rule fastqc_clean:
    input: expand("clean/{sample}_{read}_trimmed.fastq.gz", sample=SAMPLES, read=READS)
    output:
        "fastqc_clean/{sample}_{read}_trimmed_fastqc.html",
        "fastqc_clean/{sample}_{read}_trimmed_fastqc.zip"
    shell:'''
    fastqc -o fastqc_clean {input}
    '''

rule multiqc_clean:
    input: expand("fastqc_clean/{sample}_{read}_trimmed_fastqc.html", sample=SAMPLES, read=READS)
    output: "fastqc_clean/multiqc_report.html"
    shell:'''
    multiqc -o fastqc_clean fastqc_clean
    '''
rule index:
    input:
         fa='chr19_20Mb.fa',
         gtf='chr19_20Mb.gtf' 
    output:
            directory('STAR') # you can rename the index folder
    threads: 20 # set the maximum number of available cores
    shell:
            'mkdir {output} && '
            'STAR --runThreadN {threads} '
            '--runMode genomeGenerate '
            '--genomeDir {output} '
            '--genomeFastaFiles {input.fa} '
            '--sjdbGTFfile {input.gtf} '
            '--sjdbOverhang 100'

rule map_star:
    input:
        R1="clean/{sample}_R1_001_trimmed.fastq.gz",
        R2="clean/{sample}_R2_001_trimmed.fastq.gz",
        index=directory("STAR")
    output:
        bam='mapped/{sample}.bam'
    params:
        prefix = 'mapped/{sample}.',
        starlogs = 'mapped/starlogs'
    threads: 16
    shell:
        r'''
        STAR --runThreadN {threads}\
             --genomeDir {input.index}\
             --outFileNamePrefix {params.prefix} --readFilesIn {input.R1},{input.R2}\
             --outSAMtype BAM SortedByCoordinate\
             --readFilesCommand zcat\
             --limitBAMsortRAM 50000000000 &&\
             mv {params.prefix}Aligned.sortedByCoord.out.bam {output.bam} &&\
             mkdir -p {params.starlogs} &&\
             mv {params.prefix}Log.final.out {params.prefix}Log.out {params.prefix}Log.progress.out {params.starlogs}
        '''


rule sam_index:
    input:
        'mapped/{sample}.bam'
    output:
        'mapped/{sample}.bam.bai'
    shell:
        'samtools index {input}'    
        
        
        
rule feature_counts:
    input:
        bam="mapped/{sample}.bam",
        bai="mapped/{sample}.bam.bai",
        gtf="chr19_20Mb.gtf"
    output:
        "counts/{sample}.txt"
    
    threads: 5
    shell:"""
	featureCounts -T 5 -p -t exon -g gene_id -a {input.gtf} -o {output} {input.bam}
    """

rule matrix_counts:
    input:
        bam=expand("mapped/{sample}.bam", sample=SAMPLES),
        gtf="chr19_20Mb.gtf"
    output:
        "counts/matrix.txt"
    
    threads: 5
    shell:"""
	featureCounts -T 5 -p -t exon -g gene_id -a {input.gtf} -o {output} {input.bam}
    """
     
          
