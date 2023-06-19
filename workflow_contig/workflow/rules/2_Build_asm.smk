# 2nd task is to build the primary and alternate assembly and purge haploitigs if polyploid




rule hifiasm:
    """Execute Hifiasm to generate primary and alternative assemblies."""

    input:
        hifi_fastq = config['hifi']
    output:
        os.path.join(config['snakemake_dir_path'],"results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.p_ctg.gfa"),
        config['snakemake_dir_path'] + '/results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.a_ctg.gfa' if config['ploidy'] > 1 else []
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_Build_asm/hifiasmL/hifiasmL{purgel}.log')
    conda:
        '../envs/hifiasm.yaml'
    threads: 20
    resources:
        mem_mb = 70000
    params:
        runtime = '40:00:00',
        ploidy = config['ploidy']
    shell:
        """
        (hifiasm -o results/2_Build_asm/hifiasm/purgeL{wildcards.purgel}/hifiasmL{wildcards.purgel}.asm --primary -t {threads} -f 0 -l {wildcards.purgel} {input.hifi_fastq}) 2> {log}
        gzip results/2_Build_asm/hifiasm/purgeL{wildcards.purgel}/*bin
        """


rule hicanu:
    """Execute HiCanu to generate primary and alternative assemblies."""

    input:
        hifi_fastq = config['hifi']
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hicanu'))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_Build_asm/hicanu/hicanu.log')
    conda:
        '../envs/canu.yaml'
    threads: 20
    resources:
        mem_mb = 70000
    params:
        runtime = '40:00:00',
        genomeSize = config['genome_sie_estimate']
    shell:
        """
        (canu -p asm -d {output.outdir} genomeSize={params.genomeSize} -pacbio-hifi {input.hifi_fastq}) 2> {log}
        """



rule gfa_to_fasta:
    """Convert the GFA output from Hifiasm into a FASTA format."""
    input:
        asm_p = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.p_ctg.gfa'),
        asm_a = config['snakemake_dir_path'] + '/results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.a_ctg.gfa' if config['ploidy'] > 1 else []
    output:
        asm_p_fa = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.p_ctg.fasta'),
        asm_a_fa = config['snakemake_dir_path'] + '/results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.a_ctg.fasta' if config['ploidy'] > 1 else []
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_Build_asm/gfa_to_fa/gfa_to_fa_L{purgel}.log')
    threads: 1
    resources:
        mem_mb = 1000
    params:
        ploidy = config['ploidy'],
        runtime = '03:00:00'
    run:
        if params.ploidy == 1:
            shell("""awk '/^S/{{print ">"$2;print $3}}' {input.asm_p} > {output.asm_p_fa} 2> {log}""")
        else:
            shell("""awk '/^S/{{print ">"$2;print $3}}' {input.asm_p} > {output.asm_p_fa} 2> {log};
                awk '/^S/{{print ">"$2;print $3}}' {input.asm_a} > {output.asm_a_fa} 2>> {log}
                """)


rule merqury_1:
    """ Execute Merqury following Hifiasm to evaluate the quality and accuracy of the assembled genome."""

    input:
        asm_p = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.p_ctg.gfa'), 
        asm_a = config['snakemake_dir_path'] + '/results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.a_ctg.gfa' if config['ploidy'] > 1 else [],
        meryl_db = directory(os.path.join(config['snakemake_dir_path'],'results/1_PreProcessing/output_meryl/output.meryl'))
    output:
        os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/merqury_after_hifiasm/purgeL{purgel}/after_hifiasm_merqOutput.qv'),
        os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/merqury_after_hifiasm/purgeL{purgel}/after_hifiasm_merqOutput.completeness.stats')
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_Build_asm/merqury_after_hifiasm/purgeL{purgel}.log')
    conda:
        '../envs/meryl_merqury.yaml'
    threads: 2
    resources: 
        mem_mb = 5000
    params:
        out_dir = 'results/2_Build_asm/merqury_after_hifiasm/purgeL{purgel}', 
        ploidy = config['ploidy'],
        runtime = '10:00:00',
        out_prefix = 'after_hifiasm_merqOutput'
    shell:
        '''
        if [ {params.ploidy} -eq 1 ]
        then
            cd {params.out_dir}
            (merqury.sh {input.meryl_db} {input.asm_p} {params.out_prefix}) 2> {log}
        else
            cd {params.out_dir}
            (merqury.sh {input.meryl_db} {input.asm_p} {input.asm_a} {params.out_prefix}) 2> {log}
        fi
        '''

rule busco_1:
    """ Execute BUSCO following Hifiasm to assess the completeness of the assembled genome based on conserved genes."""
        
    input: 
        asm_p = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.p_ctg.fasta')
    output: 
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/busco_after_hifiasm/hifiasmL{purgel}'))
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_Build_asm/busco_after_hifiasm/busco_hifiasmL{purgel}.log')
    conda:
        '../envs/busco5.yaml'
    threads: 10
    resources:
        mem_mb = 30000
    params:
        out_path =  os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/busco_after_hifiasm'),
        out_name = 'hifiasmL{purgel}',
        lineage = config['busco_phylum'],
        runtime = '4:00:00',
        dir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/purgeL{purgel}'))
    shell:
        """
        cd {params.dir}
        (busco -f -m genome -i {input.asm_p} -o {params.out_name} --out_path {params.out_path} -l {params.lineage} -c {threads}) 2> {log}
        rm -r busco_downloads/
        """


rule quast_1:
    """ run merqury after hifiasm"""

    input:
        asm_p = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.p_ctg.fasta')
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/quast_after_hifiasm/quastL{purgel}'))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_Build_asm/quast_after_hifiasm/quast_hifiasmL{purgel}.log')
    conda:
        '../envs/quast.yaml'
    threads: 2
    resources: 
        mem_mb = 10000
    params: 
        runtime = '5:00:00'
    shell:
        '''
        (quast -o {output.outdir} -t {threads} {input.asm_p}) 2> {log}
        '''




# rule purge_dups:
#     input:
#         asm_p = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/hifiasm.asm.p_ctg.fasta'),
#         hifi_fastq = config['hifi']
#     output:
#         fa_asm_split = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/purge_dups/o_purge_dups.asm.bp.p_ctg/seqs/o_purge_dups.asm.bp.p_ctg.purged.fa')
#     log:
#         os.path.join(config['snakemake_dir_path'],'logs/2_Build_asm/purge_dups.log')
#     conda:
#         '../envs/runner.yaml'
#     threads: 12
#     resources: 
#         mem_mb = 200000
#     params:
#         pd_config = os.path.join(config['path_purge_dups'], 'scripts/pd_config.py'),
#         purge_dups = os.path.join(config['path_purge_dups'], 'scripts/run_purge_dups.py'),
#         src = os.path.join(config['path_purge_dups'], 'src'),
#         out_dir = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/purge_dups'),
#         runtime = '70:00:00'
#     shell:
#         """
#         cd {params.out_dir}
#         (echo "{input.hifi_fastq}" >> pbfofn) &> {log}
#         (python3 {params.pd_config} -n config.json {input.asm_p} pbfofn) &>> {log}
#         (python3 {params.purge_dups} -p bash config.json {params.src} o_purge_dups) &>> {log}
#         """





rule purge_dups:
    input:
        asm_p = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.p_ctg.fasta'),
        hifi_fastq = config['hifi']
    output:
        fa_asm_split = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/purge_dups/purgeL{purgel}/pri_asm.split'),
        paf_asm_split = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/purge_dups/purgeL{purgel}/pri_asm.split.self.paf.gz'),
        paf_read = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/purge_dups/purgeL{purgel}/output.paf.gz'),
        dups_bed = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/purge_dups/purgeL{purgel}/dups.bed'), 
        output_f = os.path.join(config['snakemake_dir_path'],'results/2_Build_asm/purge_dups/purgeL{purgel}/purged.fa')
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_Build_asm/purge_dups/purge_dupsL{purgel}.log')
    conda:
        '../envs/purge_dups.yaml'
    threads: 32
    resources: 
        mem_mb = 150000
    params: 
        runtime = '70:00:00',
        out_dir = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/purge_dups/purgeL{purgel}'),
        pbstat = os.path.join(config['path_purge_dups'], 'bin/pbcstat'),
        calcuts = os.path.join(config['path_purge_dups'], 'bin/calcuts'),
        hist = os.path.join(config['path_purge_dups'], 'scripts/hist_plot.py'),
        split = os.path.join(config['path_purge_dups'], 'bin/split_fa'),
        purge_dups = os.path.join(config['path_purge_dups'], 'bin/purge_dups'),
        get_seqs = os.path.join(config['path_purge_dups'], 'bin/get_seqs')
    shell:
        """
        cd {params.out_dir}
        (minimap2 -t {threads} -xasm20 {input.asm_p} {input.hifi_fastq} | gzip -c - > {output.paf_read}) 2> {log}
        ({params.pbstat} {output.paf_read}) 2>> {log}
        ({params.calcuts} PB.stat > cutoffs) 2>> {log}
        ({params.hist} -c cutoffs PB.stat purge_dups_cutoff.png) 2>> {log}
        ({params.split} {input.asm_p} > {output.fa_asm_split}) &>> {log}
        (minimap2 -t {threads} -xasm5 -DP {output.fa_asm_split} {output.fa_asm_split} | gzip -c - > {output.paf_asm_split}) 2>> {log}
        ({params.purge_dups} -2 -T cutoffs -c PB.base.cov {output.paf_asm_split} > {output.dups_bed}) 2>> {log}
        ({params.get_seqs} -e {output.dups_bed} {input.asm_p}) 2>> {log}
        """



use rule busco_1 as busco_2 with:
    input:
        asm_p = os.path.join(config['snakemake_dir_path'],'results/2_Build_asm/purge_dups/purgeL{purgel}/purged.fa')
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/busco_after_purge_dups/hifiasmL{purgel}'))
    params:
        out_path = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/busco_after_purge_dups'),
        out_name = 'hifiasmL{purgel}',
        lineage = config['busco_phylum'],
        runtime = '4:00:00',
        dir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/purge_dups/purgeL{purgel}'))
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_Build_asm/busco_after_purge_dups/hifiasmL{purgel}.log')



use rule quast_1 as quast_2 with:
    input:
        asm_p = os.path.join(config['snakemake_dir_path'],'results/2_Build_asm/purge_dups/purgeL{purgel}/purged.fa')
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'],'results/2_Build_asm/quast_after_purge_dups/quastL{purgel}'))
    log:
        os.path.join(config['snakemake_dir_path'], 'logs/2_Build_asm/quast_after_purge_dups/quastL{purgel}.log')


use rule merqury_1 as merqury_2 with:
    input:
        asm_p = os.path.join(config['snakemake_dir_path'],'results/2_Build_asm/purge_dups/purgeL{purgel}/purged.fa'),
        asm_a = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/hifiasm/purgeL{purgel}/hifiasmL{purgel}.asm.a_ctg.fasta' if config['ploidy'] > 1 else ""),
        meryl_db = directory(os.path.join(config['snakemake_dir_path'],'results/1_PreProcessing/output_meryl/output.meryl'))
    output:
        os.path.join(config['snakemake_dir_path'],'results/2_Build_asm/merqury_after_purgedups/purgeL{purgel}/after_purgedups_merqOutput.qv')
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_Build_asm/merqury_after_purgedups/purgeL{purgel}.log')
    params:
        out_dir = 'results/2_Build_asm/merqury_after_purgedups/purgeL{purgel}',
        ploidy = config['ploidy'],
        out_prefix = 'after_purgedups_merqOutput',
        runtime = '10:00:00'



rule prepare_mitohifi:
    """Mitohifi do not process compressed file (and only fasta !) and docker do not allows to manipulate file beside working directory"""
    input:
        hifi_fastq = config['hifi']
    output:
        hifi_fa = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/MitoHifi/hifi.fasta')
    threads: 2
    resources:
        mem_mb = 10000
    params:
        runtime = '4:00:00',
    shell:
        """
        zcat {input.hifi_fastq} | awk '{{if(NR%4==1) {{printf(">%s\\n",substr($0,2));}} else if(NR%4==2) print;}}' > {output.hifi_fa}
        """



rule MitoHifi:
    """ build mito asm from raw reads (better than from asm !), peut etre a deplacer dans decontam comme apres il faut filtrer les reads """

    input:
        hifi_fa = os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/MitoHifi/hifi.fasta')
    output:
        outdir = directory(os.path.join(config['snakemake_dir_path'], 'results/2_Build_asm/MitoHifi/out_ref'))
    log:
        os.path.join(config['snakemake_dir_path'],'logs/2_Build_asm/MitoHifi/MitoHifi.log')
    threads: 10
    resources:
        mem_mb = 20000
    params:
        runtime = '25:00:00',
        email = "sagane.joye@unil.ch",
        species_name = config["species_name"]
    singularity:
        'docker://ghcr.io/marcelauliano/mitohifi:master'
    shell:
        """
        findMitoReference.py --species "{params.species_name}" --email "{params.email}" --outfolder {output.outdir}
        fasta=$(ls {output.outdir}/*.fasta)
        gb=$(ls {output.outdir}/*.gb)
        cd results/2_Build_asm/MitoHifi
        mitohifi.py -r {input.hifi_fa} -f $fasta -g $gb -t {threads} -o 1
#        rm {input.hifi_fa}
        """
