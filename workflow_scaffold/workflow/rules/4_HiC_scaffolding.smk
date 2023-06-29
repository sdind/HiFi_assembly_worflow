
rule filter_5: 
    """Index the reference genome, convert FASTQ files to BAM format, and filter reads based on the 5' end (from Arima pipeline)"""

    input:
        decont_asm = config['prim_asm'],
        hic_1 = config['hic_1'],
        hic_2 = config['hic_2']
    output:
        hic_1_bam = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_1.bam'),
        hic_2_bam = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_2.bam'),
        hic_1_filt = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filtered_1.bam'),
        hic_2_filt = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filtered_2.bam')
    priority: 50
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/4_HiC_scaff/filter_5/filter_5.log')
    conda:
        '../envs/filter_5.yaml'
    threads:20
    params:
        runtime = '30:00:00',
        filt_pl = os.path.join(config['snakemake_dir_path'], 'workflow/scripts/filter_five_end.pl')
    resources:
        mem_mb = 80000
    shell:
        """
        (bwa index -a bwtsw {input.decont_asm}) 2> {log}
        (bwa mem -t {threads} {input.decont_asm} {input.hic_1} | samtools view -@ {threads} -Sb - > {output.hic_1_bam}) 2>> {log}
        (bwa mem -t {threads} {input.decont_asm} {input.hic_2} | samtools view -@ {threads} -Sb - > {output.hic_2_bam}) 2>> {log}
        (samtools view -h {output.hic_1_bam} | perl {params.filt_pl} | samtools view -Sb - > {output.hic_1_filt}) 2>> {log}
        (samtools view -h {output.hic_2_bam} | perl {params.filt_pl} | samtools view -Sb - > {output.hic_2_filt}) 2>> {log}
        """


rule mappingQC_filter:
    """  HiC Pair reads & mapping quality filter and Add read group"""

    input: 
        hic_1_filt = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filtered_1.bam'),
        hic_2_filt = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filtered_2.bam'),
        decont_asm = config['prim_asm']
    output:
        tmp = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/tmp/tmp_combHiC.bam'),
        paired = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/paired_groups.bam')
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/4_HiC_scaff/mappingQC_filter/mappingQC_filter.log')
    conda:
        '../envs/picard_samtools.yaml'
    params:
        runtime = '30:00:00',
        comb = os.path.join(config['snakemake_dir_path'], 'workflow/scripts/two_read_bam_combiner.pl'),
        name = config['name']
    threads:20
    resources:
        mem_mb = 50000
    shell:
        """
        (perl {params.comb} {input.hic_1_filt} {input.hic_2_filt} samtools 10 | samtools view -bS -t {input.decont_asm} - | samtools sort -@ {threads} -o {output.tmp} -) 2> {log}
        (picard AddOrReplaceReadGroups --INPUT {output.tmp} --OUTPUT {output.paired} -ID {params.name} -LB {params.name} -SM rep_1 -PL ILLUMINA -PU none) 2>> {log}
        """


rule mark_dup:
    """Mark duplicates"""

    input:
        paired = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/paired_groups.bam')
    output:
        dedup = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filt_dedup.bam'),
        metrics = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filt_metrics.txt'),
        stats = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filt_dedup.bam.stats')
    params:
        runtime = '30:00:00',
        tmp = directory(os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/tmp')),
        stats_pl = os.path.join(config['snakemake_dir_path'], 'workflow/scripts/get_stats.pl')
    threads:2
    resources:
        mem_mb = 50000
    log:
        'logs/4_HiC_scaff/mark_dup/mark_dup.log'
    conda:
        '../envs/picard_samtools.yaml'
    shell:
        """
        (picard MarkDuplicates --INPUT {input.paired} --OUTPUT {output.dedup} --METRICS_FILE {output.metrics} --TMP_DIR {params.tmp} --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true) 2> {log}
        (samtools index {output.dedup}) 2>> {log}
        (perl {params.stats_pl} {output.dedup} > {output.stats}) 2>> {log}
        """


rule raw_HiC_QC:
    """ From PhaseGenomics, QC method for Hi-C libraries, based on reads in a BAM file aligned to some genome/assembly"""

    input:
        decont_asm = config['prim_asm'],
        hic_1 = config['hic_1'],
        hic_2 = config['hic_2']
    output: 
        aln_bam = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/QC_raw_HiC/QC_aln.bam')
    priority: 1
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/4_HiC_scaff/QC_raw_HiC/QC_raw_HiC.log')
    conda:
        '../envs/QC_raw_HiC.yaml'
    threads:8
    params:
        runtime = '40:00:00',
        QC_script = os.path.join(config['snakemake_dir_path'], 'workflow/scripts/hic_qc.py')
    resources:
        mem_mb = 50000
    shell:
        """
        bwa mem -5SP -t {threads} {input.decont_asm} {input.hic_1} {input.hic_2} | samblaster | samtools view -S -h -b -F 2316 > {output.aln_bam}
        cd results/4_HiC_scaff/QC_raw_HiC
        (python3 {params.QC_script} -b {output.aln_bam} -o QC_raw) 2>> {log}
        """



rule yahs:
    """ run yahs for scaffolding"""

    input:
        decont_asm = config['prim_asm'],
        dedup = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filt_dedup.bam')
    output:
        dedup_sort = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC/HiC_filt_dedup_sorted.bam'),
        yahs = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out_scaffolds_final.fa'),
        yahs_bin = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out.bin'),
        yahs_agp = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out_scaffolds_final.agp')
    params:
        runtime = '08:00:00',
        dir_filt = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/filter_HiC'),
        dir_yahs = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs')
    threads:12
    resources:
        mem_mb = 50000
    log:
        os.path.join(config['snakemake_dir_path'],'logs/4_HiC_scaff/yahs/yahs.log')
    conda:
        '../envs/yahs.yaml'
    shell:
        """
        cd {params.dir_filt}
        (samtools faidx {input.decont_asm}) 2> {log}
        (samtools sort -n -o HiC_filt_dedup_sorted.bam -@ {threads} {input.dedup}) 2>> {log}
        cd {params.dir_yahs}        
        (yahs {input.decont_asm} {output.dedup_sort}) 2>> {log}
        """



rule HiC_map:
   """ Run Juicer to scaffold the assembled genome using Hi-C """

    input: 
        yahs_bin=os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out.bin'),
        yahs = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out_scaffolds_final.fa'),
        yahs_agp = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out_scaffolds_final.agp'),
        decont_asm = config['prim_asm']
    output:
        chrom_size =  os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/output_juicer/yahs_scaffolds_final.chrom.sizes')
    params: 
        runtime = '40:00:00',
        outdir_juicer = directory('results/4_HiC_scaff/output_juicer'),
        #decont_asm_ind = os.path.join(config['snakemake_dir_path'], 'results/3_decontam/asm_decontaminated.fasta.fai'),
        decont_asm_ind=expand("{prim_asm}.fai", prim_asm=config["prim_asm"]),

        yahs_fa_ind = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out_scaffolds_final.fa.fai'),
        juicer = config['path_juicer']
    threads:24
    log:
       os.path.join(config['snakemake_dir_path'],'logs/4_HiC_scaff/hic_map/hic_map.log')
    conda:
        '../envs/hic_map.yaml'
    resources:
        mem_mb = 30000
    shell:
       """
        (samtools faidx {input.decont_asm}) 2> {log}
        (samtools faidx {input.yahs}) 2>> {log}
        (cut -f1-2 {params.yahs_fa_ind} > {output.chrom_size}) 2>> {log}
        cd {params.outdir_juicer}
        ((juicer pre {input.yahs_bin} {input.yahs_agp} {params.decont_asm_ind} | sort -k2,2d -k6,6d -T ./ --parallel={threads} -S50G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)) 2>> {log}
        ((java -jar -Xmx100G {params.juicer} pre --threads {threads} alignments_sorted.txt contact_matrix.hic.part {output.chrom_size}) && (mv contact_matrix.hic.part contact_matrix.hic)) 2>> {log}
        """



rule busco:
    """ Execute BUSCO to assess the completeness of the assembled genome based on conserved genes."""

    input: 
        asm_p = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out_scaffolds_final.fa')
    output: 
        outdir =  directory(os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/busco_after_scaff'))
    log:
        os.path.join(config['snakemake_dir_path'],'logs/4_HiC_scaff/busco_after_scaff/busco_after_scaff.log')
    conda:
        '../envs/busco5.yaml'
    threads: 10
    resources:
        mem_mb = 30000
    params:
        out_path = 'results/4_HiC_scaff',
        out_name = 'busco_after_scaff',
        lineage = config['busco_phylum'],
        runtime = '4:00:00'
    shell:
        """
        (busco -f -m genome -i {input.asm_p} -o {params.out_name} --out_path {params.out_path} -l {params.lineage} -c {threads}) 2> {log}
        rm -r busco_downloads
        """



rule quast:
    """ run quast to evaluate the scaffolded assembly"""

    input:
        asm_p = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out_scaffolds_final.fa')
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/quast_after_scaff'))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/4_HiC_scaff/quast/quast.log')
    conda:
        '../envs/quast.yaml'
    threads: 2
    resources: 
        mem_mb = 10000
    params: 
        runtime = '3:00:00'
    shell:
        '''
        (quast -o {output.outdir} -t {threads} {input.asm_p}) 2> {log}
        '''



rule merqury:
    """ Execute Merqury to evaluate the quality and accuracy of the assembled genome."""

    input:
        asm_p = os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/yahs/yahs.out_scaffolds_final.fa'),
        meryl_db = config["meryl_db"]
    output:
        os.path.join(config['snakemake_dir_path'], 'results/4_HiC_scaff/merqury_after_scaff/after_scaff_merqOutput.yahs.out_scaffolds_final.qv')
    log:
        os.path.join(config['snakemake_dir_path'],'logs/4_HiC_scaff/merqury/merqury.log')
    conda:
        '../envs/meryl_merqury.yaml'
    threads: 2
    resources: 
        mem_mb = 5000
    params:
        out_dir = 'results/4_HiC_scaff/merqury_after_scaff',
        out_prefix = 'after_scaff_merqOutput',
        runtime = '10:00:00'
    shell:
        '''
        cd {params.out_dir}
        (merqury.sh {input.meryl_db} {input.asm_p} {params.out_prefix}) 2> {log}
        '''

