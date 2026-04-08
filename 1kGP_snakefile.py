"""
This snakefile aligns high coverage 1KGP samples to the HG002 primary sequence. 
It starts from cram files that are aligned to GRCh38.

This snakefile will not work out of the box. File paths must be changed appropriately. 

This snakefile also requires the HG002 sdust and trf bed file. That is created in the "T2T_snakefile.py"

"""
import pdb
import numpy as np
import pandas as pd

chroms = range(1,23)

def get_theta_MSMC2(wildcards):
    if wildcards.theta == 'default':
        return ''
    else:
        return f'-m {float(wildcards.theta)/2}'
def get_input_cram(wildcards):
    return f'/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/crumble/GRCh38/{samples_accession[wildcards.sample]}.cram.crumble'



metadata = pd.read_csv('/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/metadata/1000G_2504_high_coverage.sequence.index',comment='#',sep='\t',header=None)
samples = ['HG01882',
 'NA19625',
 'HG03006',
 'HG02373',
 'NA12718',
 'NA18530',
 'HG00443',
 'HG01250',
 'HG03515',
 'HG00266',
 'HG00118',
 'NA20845',
 'HG02568',
 'HG01783',
 'HG03977',
 'NA18939',
 'HG02113',
 'NA19017',
 'HG03212',
 'NA19648',
 'HG02285',
 'HG03234',
 'HG01171',
 'HG03753',
 'NA20752',
 'NA18488']
samples_accession = {}
for sample in samples:
    accession = metadata[metadata.iloc[:,14]==sample].iloc[:,2].item()
    samples_accession[sample] = accession

rule all:
    input:    
        # expand(
            # "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.markdup.bam.bai",
            # sample = samples
        # )
        # expand(
        #     '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.260218_4beds.chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.mhs',
        #     sample = samples,
        #     chrom = chroms, 
        #     m = [20], 
        #     sdust_w = [64],
        #     sdust_t = [20],
        # )
        expand(
            '/home/tc557/rds/hpc-work/MSMC2_1KGP/inference260221_alignedtoHG002/MSCM2/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}_m{m}_sdust{sdust_w}.{sdust_t}.final.txt',
            sample = samples,
            chrom = chroms,
            m = [20],
            sdust_w = [64],
            sdust_t = [20],
            D = [64],
            iterations = [40],
            theta = ['default'],
            rhoovermu = [0.25]
        )

ruleorder: collate > align_r1_and_r2 > align_singletons > samtools_merge > sort1_bam > samtools_fixmate > samtools_coord_sort > samtools_markdup > samtools_index > get_mean_coverage_chr20 > convert_bam_to_vcf_and_bed > write_variant_summary > write_mhs > run_MSMC2_1KGP

rule collate:
    input:
        cramfile = get_input_cram,
        reference = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/HGDP/downloaded_220412/GRCh38_full_analysis_set_plus_decoy_hla.fa'
    output:
        r1 = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1.fq.gz",
        r2 = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R2.fq.gz",
        singletons = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/singletons.fq.gz",
    params:
        other = lambda wildcards: f"/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{wildcards.sample}/other.fq.gz"
    resources:
        partition="icelake",
        time="02:00:00",
        account="DURBIN-SL2-CPU",
        mem="200G",
        cores=8
    shell:
        "samtools collate -@ 8 -O --reference {input.reference} {input.cramfile}  -T /home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/TMP/{wildcards.sample} | samtools fastq -@ 8 -1 {output.r1} -2 {output.r2} -s {output.singletons} -0 {params.other} -n -"


rule align_r1_and_r2:
    input:
        r1 = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1.fq.gz",
        r2 = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R2.fq.gz",
        reference = "/home/tc557/rds/hpc-work/T2Tsequences/HG002.pri.fasta"
    output:
        aligned_r1r2bam = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2.bam"
    resources:
        partition="icelake",
        time="06:00:00",
        account="DURBIN-SL2-CPU",
        mem="40G",
        cores=16
    params:
    shell:
        "minimap2 -t {resources.cores} -ax sr --MD {input.reference} {input.r1} {input.r2} | samtools sort -@ {resources.cores} -o {output.aligned_r1r2bam} -"
        
rule align_singletons:
    input:
        aligned_r1r2bam = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2.bam",
        singletons = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/singletons.fq.gz",
        reference = "/home/tc557/rds/hpc-work/T2Tsequences/HG002.pri.fasta"
    output:
        aligned_singletonsbam = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/singletons.bam"
    resources:
        partition="icelake",
        time="06:00:00",
        account="DURBIN-SL2-CPU",
        mem="40G",
        cores=16
    params:
    shell:
        "minimap2 -t {resources.cores} -ax sr --MD {input.reference} {input.singletons} | samtools sort -@ {resources.cores} -o {output.aligned_singletonsbam} -"

rule samtools_merge:
    input:
        aligned_singletonsbam = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/singletons.bam",
        aligned_r1r2bam = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2.bam"
    output:
        aligned_merged = temp("/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.bam")
    resources:
        partition="icelake",
        time="02:00:00",
        account="DURBIN-SL2-CPU",
        mem="20G",
        cores=8
    params:
    shell:
        "samtools merge -@ {resources.cores} -o {output.aligned_merged} {input.aligned_singletonsbam} {input.aligned_r1r2bam}"

rule sort1_bam:
    input:
        aligned_merged = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.bam"
    output:
        aligned_merged_sorted = temp("/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.bam")
    resources:
        partition="icelake",
        time="02:00:00",
        account="DURBIN-SL2-CPU",
        mem="20G",
        cores=8
    params:
    shell:
        "samtools sort -n -@ {resources.cores} -o {output.aligned_merged_sorted} {input.aligned_merged}"

rule samtools_fixmate:
    input:
        aligned_merged_sorted = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.bam"
    output:
        fixmate = temp("/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.bam")
    resources:
        partition="icelake",
        time="02:00:00",
        account="DURBIN-SL2-CPU",
        mem="20G",
        cores=8
    params:
    shell:
        "samtools fixmate -m -@ {resources.cores} {input.aligned_merged_sorted} {output.fixmate}"

rule samtools_coord_sort:    
    input:
        fixmate = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.bam"
    output:
        coord = temp("/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.bam")
    resources:
        partition="icelake",
        time="02:00:00",
        account="DURBIN-SL2-CPU",
        mem="20G",
        cores=8
    params:
    shell:
        "samtools sort -@ {resources.cores} -o {output.coord} {input.fixmate}"

rule samtools_markdup:
    input:
        coord = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.bam"
    output:
        markdup = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.markdup.bam"
    resources:
        partition="icelake",
        time="02:00:00",
        account="DURBIN-SL2-CPU",
        mem="20G",
        cores=8
    params:
    shell:
        "samtools markdup -@ {resources.cores} {input.coord} {output.markdup}"

rule samtools_index:
    input:
        markdup = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.markdup.bam"
    output:
        markdup_idx = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.markdup.bam.bai"
    resources:
        partition="icelake",
        time="02:00:00",
        account="DURBIN-SL2-CPU",
        mem="20G",
        cores=8 # testing
    params:
    shell:
        "samtools index -@ {resources.cores} {input.markdup}"


rule get_mean_coverage_chr20:
    input:
        bam = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.markdup.bam",
        bai = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.markdup.bam.bai"
    output:
        meancovchr20 = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260221meancovchr20.txt"
    resources:
        partition="icelake",
        time="00:30:00",
        account="DURBIN-SL2-CPU",
        mem="4G",
        cores=1
    params:
    shell:
        "samtools depth -r chr20_MATERNAL {input.bam} | awk '{{sum += $3}} END {{print sum / NR}}' > {output.meancovchr20}"

def get_coverage(wildcards):
    filename = f"/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{wildcards.sample}/260221meancovchr20.txt"
    with open(filename,'r') as f:
        lines = f.readlines()
    return lines[0].split('\n')

rule convert_bam_to_vcf_and_bed:
    input:
        bam = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.markdup.bam",
        bai = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/R1R2_singletons_merged.sorted.fixmate.coord.markdup.bam.bai",
        reference = "/home/tc557/rds/hpc-work/T2Tsequences/HG002.pri.fasta",
        meancovchr20 = "/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260221meancovchr20.txt"
    output: 
        vcffile = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.vcf.gz',
        bedfile = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.bed.gz',
    resources:
        partition="icelake",
        time="00:45:00",
        account="DURBIN-SL2-CPU",
        mem="4G",
        cores=1
    params:
        coverage = lambda wildcards: get_coverage(wildcards)
    shell:
        'module load bcftools; echo loaded bcftools!; bcftools mpileup -q 20 -Q 20 -C 50 -r chr{wildcards.chrom}_MATERNAL -f {input.reference} {input.bam} | bcftools call -c | /home/tc557/T2T_PSMC/T2Tsequences/bamCaller.py {params.coverage} {output.bedfile} | bgzip -c > {output.vcffile}'

rule write_variant_summary:
    input:
        vcffile = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.vcf.gz',
    output:
        variant_summary = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.SNP_summary.txt.gz',
        indel_bed = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.indel.bed.gz'
    resources:
        partition="icelake",
        time="00:20:00",
        account="DURBIN-SL2-CPU",
        mem="16G",
        cores=1
    params:
        indel_bed = lambda wildcards: f'/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{wildcards.sample}/260220chr{wildcards.chrom}.indel.bed.gz'
    shell:
        "module load bcftools; echo \"loaded bcftools\"; bcftools index {input.vcffile} ; bcftools view {input.vcffile} -H -r chr{wildcards.chrom}_MATERNAL -i 'GT==\"0/1\"' -v snps | awk '{{print $2,$4$5}}' | gzip -c > {output.variant_summary}; "
        "bcftools view {input.vcffile} -H -r chr{wildcards.chrom}_MATERNAL -v indels | awk '{{print \"chr{wildcards.chrom}\"\"\\t\"$2-5\"\\t\"$2+5}}' | gzip -c > {output.indel_bed}"

rule write_mhs:
    input:
        variant_summary = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.SNP_summary.txt.gz',
        # indel_bed = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.indel.bed.gz',
        bed = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.bed.gz',
        trf_bed = '/home/tc557/rds/hpc-work/T2Tsequences/HG002.pri.chr{chrom}.fasta.trf.bed',
        sdust_bed = '/home/tc557/rds/hpc-work/T2Tsequences/HG002.pri.chr{chrom}.fasta.sdust.w{sdust_w}.t{sdust_t}dat',
    output: 
        mhsfile = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.mhs',
        log = '/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{sample}/260220chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.log.txt'
    params:
        indel_bed = lambda wildcards: f'/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{wildcards.sample}/260220chr{wildcards.chrom}.indel.bed.gz'
        # centromere_start = lambda wildcards: int(centromeres[centromeres[0]==f'chr{wildcards.chrom}'][1]) - centromere_pad,
        # centromere_end = lambda wildcards: int(centromeres[centromeres[0]==f'chr{wildcards.chrom}'][2]) + centromere_pad,
    resources:
        partition="icelake",
        time='00:30:00',
        account = "DURBIN-SL2-CPU",
        mem="20G",
        cores=1
    # shell:
        # "/home/tc557/T2T_PSMC/T2Tsequences/251228_make_mhs.py -iv {input.SNP_pos_summary} -negative_beds {input.indel_bed} {input.trf_bed} {input.sdust_bed} -positive_beds {input.bed} -chrom {wildcards.chrom} -truncate 500000 -centromere_start {params.centromere_start} -centromere_end {params.centromere_end} -m {wildcards.m} -o {output.mhsfile} | tee {output.log}"
    shell:
        """
        python /home/tc557/T2T_PSMC/T2Tsequences/251228_make_mhs.py \
            -iv {input.variant_summary} \
            -negative_beds {params.indel_bed} {input.trf_bed} {input.sdust_bed} \
            -positive_beds {input.bed} \
            -chrom {wildcards.chrom} \
            -truncate 50000 \
            -m {wildcards.m} \
            -o {output.mhsfile} | tee {output.log}
        """

def get_mhs_files(wildcards):
    return [f'/home/tc557/rds/rds-aegis-2-XMbJZDSpoJk/projects/human/1000Genomes_30X/sequencing_reads/{wildcards.sample}/260220chr{chrom}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.mhs' for chrom in range(1,23)]

rule run_MSMC2_1KGP:
    input:
        mhsfiles = get_mhs_files,
    output:
        outfile = '/home/tc557/rds/hpc-work/MSMC2_1KGP/inference260221_alignedtoHG002/MSCM2/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}_m{m}_sdust{sdust_w}.{sdust_t}.final.txt',
        outlog = '/home/tc557/rds/hpc-work/MSMC2_1KGP/inference260221_alignedtoHG002/MSCM2/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}_m{m}_sdust{sdust_w}.{sdust_t}.log',
    resources:
        partition="icelake-himem",
        time='01:30:00',
        account = "DURBIN-SL2-CPU",
        mem="180G",
        cores=20
    params:
        outfile = lambda wildcards: f'/home/tc557/rds/hpc-work/MSMC2_1KGP/inference260221_alignedtoHG002/MSCM2/D_{wildcards.D}/iterations_{wildcards.iterations}/theta_{wildcards.theta}/rhoovermu_{wildcards.rhoovermu}/{wildcards.sample}_m{wildcards.m}_sdust{wildcards.sdust_w}.{wildcards.sdust_t}',
        ztheta = lambda wildcards: get_theta_MSMC2(wildcards)
    shell:
        "/home/tc557/msmc2_installation_220407/msmc2-2.1.3/build/release/msmc2 -i {wildcards.iterations} -t {resources.cores} -p {wildcards.D}*1 {params.ztheta} -r {wildcards.rhoovermu} -o {params.outfile} {input.mhsfiles}"

