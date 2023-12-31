# 3nd task is to eliminate contaminants from the assembled genome



rule blast:
    """Perform a BLAST search on the assembled contigs to filter out potential contaminants by identifying taxonomic matches."""

    input:
        asm_p = config['prim_asm']
    output:
        out_blast = os.path.join(config['snakemake_dir_path'],'results/3_decontam/asm.vs.nt.max10.1e25.blastn.out')
    log:
        os.path.join(config['snakemake_dir_path'],'logs/3_decontam/blast/blast.log')
    conda:
        '../envs/blast.yaml'
    threads: 40
    resources: 
        mem_mb = 200000
    params: 
        runtime = '72:00:00',
        nt_dir = os.path.join(config['snakemake_dir_path'], 'results/3_decontam')
    shell:
        """
        cd {params.nt_dir}
        (update_blastdb.pl --decompress nt --num_threads {threads}) 2> {log}
        (blastn -query {input.asm_p} -db nt -outfmt '6 qseqid staxids bitscore evalue std sscinames sskingdoms stitle' -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads {threads} -out {output.out_blast}) 2>> {log}
        """


rule cov_estimate:
    """Align the reads to the assembled genome using Minimap2, and subsequently extract coverage information with the map2cov utility from BlobTools """

    input:
        asm_p = config['prim_asm'],
        hifi_fastq = config['hifi']
    output: 
        align_asm_hifi = os.path.join(config['snakemake_dir_path'],'results/3_decontam/align_asm_hifi.sorted.bam'),
        map2cov = os.path.join(config['snakemake_dir_path'],'results/3_decontam/map2cov_blob.align_asm_hifi.sorted.bam.cov')
    log:
        os.path.join(config['snakemake_dir_path'],'logs/3_decontam/cov_estimate/cov_estimate.log')
    conda:
        '../envs/blobtool.yaml'
    threads: 20
    resources:
        mem_mb = 100000
    params:
        runtime = '10:00:00',
        direc = os.path.join(config['snakemake_dir_path'], 'results/3_decontam')

    shell:
        """
        cd {params.direc}
        (minimap2 -t {threads} -ax map-hifi {input.asm_p} {input.hifi_fastq} | samtools sort -o {output.align_asm_hifi}) 2> {log}
        (samtools index -b {output.align_asm_hifi}) 2>> {log}
        (blobtools map2cov -i {input.asm_p} -b {output.align_asm_hifi} -o map2cov_blob) 2>> {log}
        """


rule blobtools:
    """ creates a BlobTools database using the assembled genome, BLAST output, and coverage information, then generates a blob plot and a summary table """

    input:
        out_blast = os.path.join(config['snakemake_dir_path'],'results/3_decontam/asm.vs.nt.max10.1e25.blastn.out'),
        asm_p = config['prim_asm'],
        map2cov = os.path.join(config['snakemake_dir_path'],'results/3_decontam/map2cov_blob.align_asm_hifi.sorted.bam.cov')
    output:
        blob_json = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/blob_output.blobDB.json'),
        blob_table = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/table.blob_output.blobDB.table.txt'),
        blob_png = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/plot.blob_output.blobDB.json.bestsumorder.phylum.p8.count.100.blobplot.cov0.png') 
    log:
        os.path.join(config['snakemake_dir_path'],'logs/3_decontam/blobtools/blobtools.log')
    conda:
        '../envs/blobtool.yaml'
    threads: 20
    resources:
        mem_mb = 20000
    params:
        runtime = '70:00:00',
        direc= os.path.join(config['snakemake_dir_path'], 'results/3_decontam'),
        db = os.path.join(config['snakemake_dir_path'], 'workflow/files/nodesDB.txt')
    shell:
        """
        cd {params.direc}
        (blobtools create -i {input.asm_p} -t {input.out_blast} --db {params.db} -c {input.map2cov} -x bestsumorder -o blob_output) 2> {log}
        (blobtools plot -i blob_output.blobDB.json --sort count --hist count -x bestsumorder -o plot) 2>> {log}
        (blobtools view -i blob_output.blobDB.json --hits --rank all -x bestsumorder -o table) 2>> {log}
        rm nt.*
        """




rule filter_asm:
    """ Remove contaminant contigs from the assembly """

    input:
        asm_p = config['prim_asm'],
        blob_table = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/table.blob_output.blobDB.table.txt')
    output:
        contam_scaff = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/contam_scaff.txt'),
        decont_asm = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/asm_decontaminated.fasta')
    log:
        os.path.join(config['snakemake_dir_path'],'logs/3_decontam/filter_asm/filter_asm.log')
    conda:
        '../envs/filter_asm.yaml'
    threads: 20
    resources:
        mem_mb = 5000
    params:
        runtime = '70:00:00'
    shell:
        """
        (gawk 'FNR > 10 {{if ($10 != "Arthropoda" && $10 !="no-hit" && $10 != "superkingdom.hits.9") {{print $1}}}}' {input.blob_table} > {output.contam_scaff}) &> {log}
        (filterbyname.sh in={input.asm_p} names={output.contam_scaff} out={output.decont_asm} include=f) &>> {log}
        """



rule busco_1:
    """ Execute BUSCO to assess the completeness of the assembled genome based on conserved genes."""
        
    input: 
        decont_asm = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/asm_decontaminated.fasta')
    output: 
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/3_decontam/output_busco'))
    log:
        os.path.join(config['snakemake_dir_path'],'logs/3_decontam/busco/busco.log')
    conda:
        '../envs/busco5.yaml'
    threads: 10
    resources:
        mem_mb = 30000
    params:
        out_path =  os.path.join(config['snakemake_dir_path'], 'results/3_decontam'),
        out_name = 'output_busco',
        lineage = config['busco_phylum'],
        runtime = '4:00:00',
    shell:
        """
        cd {params.out_path}
        (busco -f -m genome -i {input.decont_asm} -o {params.out_name} --out_path {params.out_path} -l {params.lineage} -c {threads}) 2> {log}
        rm -r busco_downloads/
        """


rule quast_1:
    """ run quast to evaluate the assembly"""

    input:
        decont_asm = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/asm_decontaminated.fasta')
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/3_decontam/output_quast'))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/3_decontam/output_quast/quast.log')
    conda:
        '../envs/quast.yaml'
    threads: 2
    resources: 
        mem_mb = 10000
    params: 
        runtime = '5:00:00'
    shell:
        '''
        (quast -o {output.outdir} -t {threads} {input.decont_asm}) 2> {log}
        '''



