# HiFi Genome Assembly Workflow
This repository contains a Snakemake workflow for genome assembly, consisting of three interconnected workflows: workflow_contig, workflow_decontam, and workflow_scaffold. Each workflow performs specific tasks to facilitate the assembly process and improve the quality of the final genome assembly.

An additional [repository](https://github.com/sdind/genome_annotation_workflow) also offers a snakemake workflow for genome annotation using RNA-seq and protein data. 

<br/>

## Workflow 1: workflow_contig
This workflow focuses on processing the high-fidelity (HiFi) reads and generating multiple assemblies using hifiasm and purge_dups with different purging levels, and evaluate the quality and completeness of the assembled genomes. The workflow also includes the generation of the mitochondrial genome using MitoHifi.
>**Note**: Hifiasm provides multiple levels of purge, ranging from soft to very aggressive purging. This workflow generates assemblies for each purge level and an additional assembly without any purging. The purge_dups tool is then run on all assemblies to identify and remove remaining duplicated sequences. It is worth noting that the evaluation of the assemblies before and after purge_dups has revealed that in some cases, the assembly using only the internal hifiasm purge function outperforms the assemblies with purge_dups. Therefore, it is essential to thoroughly assess the quality of each assembly before making a final decision
>

The workflow includes the following steps:
1. **Preprocessing and Reads Quality Assessment**:
  * Quality Control (QC): Perform [LongQC](https://github.com/yfukasawa/LongQC) analysis on HiFi reads to assess their quality.
  * [GenomeScope2](https://github.com/tbenavi1/genomescope2.0): Provides a comprehensive understanding of the genome’s profile.
  * [Smudgeplot](https://github.com/KamilSJaron/smudgeplot): Assess genome heterozygosity.
2. **Genome Assembly and Evaluation**:
  * Assembly Generation: Utilize [Hifiasm](https://github.com/chhylp123/hifiasm) to build several assemblies using different purge levels.
  * Purge Duplicates: Run [purge_dups](https://github.com/dfguan/purge_dups) to remove remaining duplicated sequences from the assemblies.
  * Assembly QC: Perform quality assessment on all assemblies using [Quast](https://quast.sourceforge.net/), [Merqury](https://github.com/marbl/merqury), and [Busco](https://busco.ezlab.org/) to evaluate the quality, accuracy and completeness of the assembled genomes.
  * Mitochondrial Genome Assembly: Construct the mitochondrial genome using the [MitoHifi](https://github.com/marcelauliano/MitoHiFi) tool.

**Manually choose the best assembly as input for the next workflow based on the QC results.**

<br/>

## Workflow 2: workflow_decontam
The second workflow focuses on decontaminating the selected assembly obtained from the previous workflow. It incorporates the following steps:
* BLAST-Based Decontamination: Blast the raw hifi reads against the blastn database to identify potential contaminants.
* Coverage Estimation: Estimate the coverage of the reads to identify potential contamination sources.
* [Blobtools](https://github.com/DRL/blobtools) Analysis: Utilize Blobtools to generate plots and identify reads associated with contaminants.
* Contig Filtering: Filter out contaminant contigs from the assembly based on the analysis results.

<br/>

## Workflow 3: workflow_scaffold
The third workflow is dedicated to scaffolding, specifically when Hi-C reads are available. It involves the following steps:
* Hi-C Reads QC: Perform quality control on hifi reads using the [hic_qc.py](https://github.com/phasegenomics/hic_qc) script from phasegenomics and FastQC
* Hi-C Reads Mapping and Filtering: Mapp and Filter Hi-C reads using the [ArimaGenomics mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline).
* Scaffolding: Employ [yahs](https://github.com/c-zhou/yahs) for the scaffolding process to connect contigs using the Hi-C information.
* Hi-C Map Generation: Create a Hi-C contact map using [Juicer](https://github.com/aidenlab/juicer) to visualize the scaffolding results.
* Scaffolding QC: Evaluate the scaffolded assembly using [Quast](https://quast.sourceforge.net/), [Merqury](https://github.com/marbl/merqury), and [Busco](https://busco.ezlab.org/) to ensure its quality.


### Directory Structure
All workflows share the same directory structure:
```
.
├── config.yaml    # Configuration file specifying input data and parameters
├── logs        # Log files for each step
├── results      # Directory containing output files for each step
├── README.md     # Readme file 
└── workflow      # workflow directory
  ├── Snakefile   # Main Snakemake file for the workflow
  ├── envs      # Directory for Conda environments required tools
    ├── first_task.yaml  # Conda environment file
    ├── second_task.yaml  
  ├── rules     # Directory for snakefiles implementing rules
    ├── Snakefile_1.smk  # File implementing a subset of rules
    ├── Snakefile_2.smk  
  └── scripts    # Custom scripts (if any)
  └── files    # Files required for the worflow (e.g protein dataset) (if any)
```

### Configuration File
The behavior of the assembly workflow can be customized by modifying the config.yaml file. This file contains various parameters and paths that can be adjusted to suit your specific sequencing data and analysis requirements (it may also be necessary to modify memory usage settings within Snakemake rules depending of your genome size).
Before running the workflow, make sure to update the values in the config.yaml file according to your setup.

<br/>

## Workflow_contig Configuration file
The following parameters can be customized in the configuration file:
* `hifi`: The path to the input HiFi reads file in compressed FASTQ format (`.fastq.gz`).
* `snakemake_dir_path`: The directory path where the Snakemake workflow files are located.
* `path_purge_dups`: The directory path where the purge_dups tool is installed. Please note that purge_dups needs to be manually installed and its path should be provided in this parameter. You can download and install purge_dups from the [official repository](https://github.com/dfguan/purge_dups).
* `species_name`: The complete name of the species being assembled.
* `name`: The short name of your species or assembly run
* `ploidy`: The ploidy level of the organism being assembled.
* `genome_size_estimate`: An estimate of the genome size in megabases (M). For example, ‘400m’ represents 400 megabases.
* `busco_phylum`: The BUSCO database identifier for the phylum of the organism being assembled. For example, ‘hymenoptera_odb10’ represents the hymenoptera phylum.
* `kmer_len`: The length of the k-mer to use.

Please note that the other required tools mentioned in the workflow, will be automatically installed via conda using the provided YAML file during the workflow execution. However, the purge_dups tool requires manual installation. You will need to download and install purge_dups separately from the official repository or any other reliable source. Once purge_dups is installed, provide the directory path where it is installed in the path_purge_dups parameter in the configuration file.

## Workflow_decontam Configuration file
The following parameters should be customized in the configuration file:
* `snakemake_dir_path`: The directory path where the Snakemake workflow files for the decontamination workflow are located.
* `hifi`: The path to the input Hifi reads file in compressed FASTQ format (`.fastq.gz`).
* `prim_asm`: The path to the primary assembly file obtained from the previous workflow (workflow_contig).
* `name`: The short name of your species or assembly run

## Workflow_scaffold Configuration file
The following parameters should be customized in the configuration file:
* `snakemake_dir_path`: The directory path where the Snakemake workflow files for the scaffolding workflow are located.
* `hifi`: The path to the input HiFi reads file in compressed FASTQ format (`.fastq.gz`).
* `prim_asm`: The path to the decontaminated primary assembly file obtained from the previous workflows (workflow_contig + workflow_decontam).
* `hic_1`: The path to the first Hi-C reads file. The file must ends with `_1` suffix to indicate reads 1 in the paired-end data.
* `hic_2`: The path to the second Hi-C reads file.The file must ends with `_2` suffix to indicate reads 2 in the paired-end data.
* `meryl_db`: The path to the Meryl database generated from the previous workflow (workflow_contig).
* `path_juicer`: The path to the Juicer Tools JAR file. Please note that Juicer Tools needs to be manually installed. You can download and install Juicer Tools from the [official repository](https://github.com/aidenlab/JuicerTools).
* `name`: The short name of your species or assembly run
* `genome_size_estimate`: An estimate of the genome size in megabases (MB). For example, 2000 represents 2000 megabases.
* `busco_phylum`: The BUSCO database identifier for the phylum of the organism being assembled. For example, ‘hymenoptera_odb10’ represents the hymenoptera phylum.
* `kmer_len`: The length of the k-mer to use.