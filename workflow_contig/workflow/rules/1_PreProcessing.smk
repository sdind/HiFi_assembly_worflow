# The initial task involves examining the quality of the read sequences and generating GenomeScope plots, which will provide a more comprehensive understanding of our target genome



rule long_qc:
    """run longqc to check reads quality"""

    input:
        hifi_fastq = config['hifi']
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/1_PreProcessing/output_LongQC'))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/1_PreProcessing/long_qc/long_qc.log')
    conda:
        '../envs/hifi_QC.yaml'
    threads:4
    resources:
        mem_mb = 90000
    params:
        runtime = '05:00:00'
    shell:
        """
        (python3 /work/FAC/FBM/DEE/rwaterho/evofun/sagane/packages/LongQC/longQC.py sampleqc -x pb-hifi -o {output.outdir} {input.hifi_fastq}) 2> {log}
        """


rule merylCount:
    """run meryl to create k-mer count hist"""
    input:
        hifi_fastq = config['hifi']
    output:
        meryl_db = directory(os.path.join(config['snakemake_dir_path'], 'results/1_PreProcessing/output_meryl/output.meryl')),
        hist = os.path.join(config['snakemake_dir_path'], 'results/1_PreProcessing/output_meryl/meryl.hist')
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/1_PreProcessing/meryl/merylCount.log')
    conda:
        '../envs/meryl_merqury.yaml'
    threads:12
    resources:
        mem_mb = 100000
    params:
        runtime = '03:00:00',
        kmer=config['kmer_len']
    shell:
        """
        (meryl count k={params.kmer} threads={threads} {input.hifi_fastq} output {output.meryl_db}) 2> {log}
        meryl histogram {output.meryl_db} > {output.hist}
        """

rule GenomeScope2:
    """ Run GenomeScope2 to obtain a comprehensive understanding of our genome's profile"""

    input: 
        hist = os.path.join(config['snakemake_dir_path'],'results/1_PreProcessing/output_meryl/meryl.hist')
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/1_PreProcessing/output_GenomeScope2'))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/1_PreProcessing/GenomeScope2/GenomeScope2.log')
    conda:
        '../envs/GenomeScope2.yaml'
    params:
        runtime = '04:00:00',
        kmer=config['kmer_len'],
        ploidy = config['ploidy']
    threads: 1
    resources:
        mem_mb = 10000
    shell:
        """
        (genomescope2 -i {input.hist} -o {output.outdir} -p {params.ploidy} -k {params.kmer} --max_kmercov=1000000) 2> {log}
        """

rule smudgeplot:
    """ Currently, smudgeplot works only with jellyfish or KMC kmer db ... """

    input:
        hifi_fastq = config['hifi']
    output:
        kmer_jf = os.path.join(config['snakemake_dir_path'],'results/1_PreProcessing/smudgeplot/kmer_counts.jf'),
        hist = os.path.join(config['snakemake_dir_path'],'results/1_PreProcessing/smudgeplot/jelly.hist'), 
        kp_cov = os.path.join(config['snakemake_dir_path'],'results/1_PreProcessing/smudgeplot/kmer_pairs_coverages.tsv')
    log:
        os.path.join(config['snakemake_dir_path'],'logs/1_PreProcessing/Smudgeplot/Smudgeplot.log')
    conda:
        '../envs/smudgeplot.yaml'
    params:
        runtime = '09:00:00',
        kmer=config['kmer_len'],
        out_dir = os.path.join(config['snakemake_dir_path'],'results/1_PreProcessing/smudgeplot')
    threads: 8
    resources:
        mem_mb = 110000
    shell:
        """
        (jellyfish count -C -m {params.kmer} -s 1000000000 -t {threads} <(zcat {input.hifi_fastq}) -o {output.kmer_jf}) 2> {log}
        cd {params.out_dir}
        (jellyfish histo {output.kmer_jf} > {output.hist}) 2>> {log}
        L=$(smudgeplot.py cutoff {output.hist} L)
        U=$(smudgeplot.py cutoff {output.hist} U)
        (echo $L $U) 2>> {log} 
        (jellyfish dump -c -L $L -U $U {output.kmer_jf} | smudgeplot.py hetkmers -o kmer_pairs) 2>> {log}
        (smudgeplot.py plot {output.kp_cov} -o my_genome) 2>> {log}
        """


