configfile: "/home/tc557/T2T_PSMC/T2Tsequences/260127snakefile/config.yaml"
"""
This snakefile starts with a haplotype resolved T2T assembly, aligns the alternate to the primary, calls variants, filters them, then runs PSMC (as in MSMC2).
Download the human T2T sequence (HG002) from: https://www.biorxiv.org/content/10.1101/2025.09.21.677443v1
Download the chimpanzee, bonobo, gorilla, Bornean orangutan and Sumatran orangutan from: https://www.nature.com/articles/s41586-025-08816-3

This snakefile will not work out of the box. File paths must be changed appropriately. 

This processes 4 beds into mhs (the multihetsep file, the summary of genetic variants and masking), these are :
    -indel bed (+ or - 5bps around indel)
    -bed from dipcall
    -bed from tandem repeat masker
    -bed from sdust (low complexity)

"""
import numpy as np
import pandas as pd
import pdb

def get_chrom_suffix(wildcards,hsa_mPanPan1,hsa_mPanTro3,hsa_mGorGor1,hsa_mPonAbe1,hsa_mPonPyg2):
    if wildcards.sample=='HG002':
        return '_MATERNAL'
    elif wildcards.sample=='mPanPan1':
        return f'_{hsa_mPanPan1[wildcards.chrom]}'
    elif wildcards.sample=='mPanTro3':
        return f'_hap1_hsa{hsa_mPanTro3[int(wildcards.chrom)]}'
    elif wildcards.sample=='mGorGor1':
        return '_' + hsa_mGorGor1[int(wildcards.chrom)]
    elif wildcards.sample=='mPonAbe1':
        return f'_hap1_hsa{hsa_mPonAbe1[int(wildcards.chrom)]}'
    elif wildcards.sample=='mPonPyg2':
        return f'_hap1_hsa{hsa_mPonPyg2[int(wildcards.chrom)]}'
    else:
        print(f'ERROR IN get_chrom_suffix')
        return None 

def get_theta(wildcards):
    if wildcards.theta == 'default':
        return ''
    else:
        return f'-theta {wildcards.theta}'

def get_theta_MSMC2(wildcards):
    if wildcards.theta == 'default':
        return ''
    else:
        return f'-m {float(wildcards.theta)/2}'

def chroms_for_sample(wildcards):
    return config["chromosomes"][wildcards.sample]

hsa_mPanPan1 = {
    '1':'pat_hsa1',
    '2':'mat_hsa3',
    '3':'pat_hsa4',
    '4':'mat_hsa5',
    '5':'pat_hsa6',
    '6':'mat_hsa7',
    '7':'pat_hsa8',
    '8':'mat_hsa10',
    '9':'pat_hsa11',
    '10':'mat_hsa12',
    '11':'mat_hsa9',
    '12':'mat_hsa2a',
    '13':'mat_hsa2b',
    '14':'pat_hsa13',
    '15':'mat_hsa14',
    '16':'pat_hsa15',
    '17':'mat_hsa18',
    '18':'pat_hsa16',
    '19':'pat_hsa17',
    '20':'pat_hsa19',
    '21':'mat_hsa20',
    '22':'pat_hsa21',
    '23':'pat_hsa22'
}

hsa_mPanTro3 = {1:1,
    2:3,
    3:4,
    4:5,
    5:6,
    6:7,
    7:8,
    8:10,
    9:11,
    10:12,
    11:9,
    12:'2a',
    13:'2b',
    14:13,
    15:14,
    16:15,
    17:18,
    18:16,
    19:17,
    20:19,
    21:20,
    22:21,
    23:22,
}

hsa_mGorGor1 = {
    1:'pat_hsa1',
    2:'pat_hsa3',
    3:'pat_hsa4',
    4:'pat_hsa17x5',
    5:'mat_hsa6',
    6:'mat_hsa7',
    7:'pat_hsa8',
    8:'pat_hsa10',
    9:'pat_hsa11',
    10:'mat_hsa12',
    11:'mat_hsa2b',
    12:'pat_hsa2a',
    13:'pat_hsa9',
    14:'pat_hsa13',
    15:'pat_hsa14',
    16:'pat_hsa15',
    17:'mat_hsa18',
    18:'pat_hsa16',
    19:'pat_hsa5x17',
    20:'mat_hsa19',
    21:'pat_hsa20',
    22:'mat_hsa21',
    23:'mat_hsa22'
}

hsa_mPonAbe1 = {1:1,
    2:3,
    3:4,
    4:5,
    5:6,
    6:7,
    7:8,
    8:10,
    9:11,
    10:12,
    11:'2b',
    12:'2a',
    13:9,
    14:13,
    15:14,
    16:15,
    17:18,
    18:16,
    19:17,
    20:19,
    21:20,
    22:21,
    23:22}

hsa_mPonPyg2 = {1:1,
    2:3,
    3:4,
    4:5,
    5:6,
    6:7,
    7:8,
    8:10,
    9:11,
    10:12,
    11:'2b',
    12:'2a',
    13:9,
    14:13,
    15:14,
    16:15,
    17:18,
    18:16,
    19:17,
    20:19,
    21:20,
    22:21,
    23:22}

Ds = [64]
iterationss = [40]
rhoovermus = [0.25]
thetas = ['default']
samples = ['HG002' ,'mPanTro3','mGorGor1','mPanPan1','mPonPyg2','mPonAbe1']
ms = [20]

sdust_wts = [(64,20)]
bootstraps = range(0,20)


rule all:
    input:
        # [f'/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.260127_4beds.chr{chrom}.m{m}.w{sdust_wt[0]}.t{sdust_wt[1]}.mhs' for sample in samples for chrom in config["chromosomes"][sample] for m in ms for sdust_wt in sdust_wts]
        # [f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_wt[0]}.t{sdust_wt[1]}.final.txt' for D in Ds for iterations in iterationss for theta in thetas for rhoovermu in rhoovermus for sample in samples for chrom in config["chromosomes"][sample] for m in ms for sdust_wt in sdust_wts] # run MSMC2
        [f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/bootstrap{bootstrap}.{sample}.m{m}.w{sdust_wt[0]}.t{sdust_wt[1]}.final.txt' for bootstrap in bootstraps for D in Ds for iterations in iterationss for theta in thetas for rhoovermu in rhoovermus for sample in samples for chrom in config["chromosomes"][sample] for m in ms for sdust_wt in sdust_wts] # run MSMC2 with bootstraps
        # [f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_wt[0]}.t{sdust_wt[1]}.decoding.chr{chrom}.txt.gz' for D in Ds for iterations in iterationss for theta in thetas for rhoovermu in rhoovermus for sample in samples for m in ms for sdust_wt in sdust_wts for chrom in config["chromosomes"][sample]] # decode
        # [f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_wt[0]}.t{sdust_wt[1]}.decoding.chr{chrom}.bootstrap{bootstrap}.txt.gz' for bootstrap in bootstraps for D in Ds for iterations in iterationss for theta in thetas for rhoovermu in rhoovermus for sample in samples for m in ms for sdust_wt in sdust_wts for chrom in config["chromosomes"][sample]] # decode bootstraps

rule run_dipcall_align:
    input:
        input_seq1 = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.fasta',
        input_seq2 = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.alt.fasta'
    output:
        outvcf = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.vcf.gz',
        outbed = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.bed'
    params:
        prefix = lambda wildcards: f'/home/tc557/rds/hpc-work/T2Tsequences/{wildcards.sample}_all_pap'
    resources:
        partition="icelake",
        time='03:00:00',
        account = "DURBIN-SL2-CPU",
        mem="200G",
        cores=2
    shell:
        "/home/tc557/dipcall.kit/run-dipcall {params.prefix} {input.input_seq1} {input.input_seq1} {input.input_seq2} > {params.prefix}.mak; make -j2 -f {params.prefix}.mak"

rule index_vcf:
    input:
        invcf = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.vcf.gz',
    output:
        outvcf = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.vcf.gz.csi'
    resources:
        partition="icelake",
        time='00:20:00',
        account = "DURBIN-SL2-CPU",
        mem="10G",
        cores=1    
    shell:
        "bcftools index {input.invcf}"

rule index_reference:
    input:
        flag1 = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.vcf.gz.csi',
        fasta = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.fasta'
    output:
        reference = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.mmi'
    resources:
        partition="icelake",
        time='00:30:00',
        account = "DURBIN-SL2-CPU",
        mem="20G",
        cores=1
    params:
    shell:
        "/home/tc557/minimap2-2.30_x64-linux/minimap2 -d {output.reference} {input.fasta}"

rule make_variant_bed_summaries:
    input:
        flag = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.mmi',
        vcf = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.vcf.gz',
        # flag = '/home/tc557/rds/hpc-work/T2Tsequences/trd_sdust_test_mappability_on_{sample}pri.chr{chrom}.sorted.bed'
    output:
        SNP_pos_summary = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.variantsummary.chr{chrom}.0based.txt.gz',
        indel_bed = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.chr{chrom}.indel.bed.gz'
    params:
        chrom_name_suffix = lambda wildcards: get_chrom_suffix(wildcards,hsa_mPanPan1,hsa_mPanTro3,hsa_mGorGor1,hsa_mPonAbe1,hsa_mPonPyg2)
    resources:
        partition="icelake",
        time='00:20:00',
        account = "DURBIN-SL2-CPU",
        mem="20G",
        cores=1
    shell:
        "bcftools view {input.vcf} -H -r chr{wildcards.chrom}{params.chrom_name_suffix} -v snps | awk '{{print $2,$4$5}}' | gzip -c > {output.SNP_pos_summary};"
        "bcftools view {input.vcf} -H -r chr{wildcards.chrom}{params.chrom_name_suffix} -v indels | awk '{{print \"chr{wildcards.chrom}\"\"\\t\"$2-5\"\\t\"$2+5}}' | gzip -c > {output.indel_bed}"

rule run_trf:
    input:
        fasta = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.fasta',
        flag = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.mmi' 
        # flag = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.chr{chrom}.fasta.sdust.dat',
    output:
        dat =  '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.chr{chrom}.fasta.trf.dat',
    resources:
        partition="icelake",
        time='09:00:00',
        account = "DURBIN-SL2-CPU",
        mem="64G",
        cores=3
    params:
        chrom_name_suffix = lambda wildcards: get_chrom_suffix(wildcards,hsa_mPanPan1,hsa_mPanTro3,hsa_mGorGor1,hsa_mPonAbe1,hsa_mPonPyg2)
    shell:
        # Extract single fasta sequence and run TRF
        # awk taken from https://onestopdataanalysis.com/get-sequence-fasta/
        # Note '-l 10' - supports <10Mbp TR arrays, likely to be enough for T2T (centromeres)
        "awk -v seq=\"chr{wildcards.chrom}{params.chrom_name_suffix}\" -v RS='>' '$1 == seq {{print RS $0}}' {input.fasta} | "
        "/home/tc557/TRF/build/trf - 2 6 6 80 10 50 500 -ngs -h -l 10 "
        "> {output.dat}"

rule write_trf_bed:
    input:
        dat =  '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.chr{chrom}.fasta.trf.dat',
    output:
        dat =  '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.chr{chrom}.fasta.trf.bed',
    resources:
        partition="icelake",
        time='00:01:00',
        account = "DURBIN-SL2-CPU",
        mem="4G",
        cores=1
    params:
    shell:
        "cat {input.dat} | awk -v CHROM={wildcards.chrom} '{{print \"chr\"CHROM\"\\t\"$1-1\"\\t\"$2-1 }}' | tail -n +2 > {output.dat}"

rule run_sdust_T2T:
    input:
        flag =  '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.chr{chrom}.fasta.trf.bed',
        fasta = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.fasta'
    output:
        dat = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.chr{chrom}.fasta.sdust.w{sdust_w}.t{sdust_t}dat',
    params:
        chrom_name_suffix = lambda wildcards: get_chrom_suffix(wildcards,hsa_mPanPan1,hsa_mPanTro3,hsa_mGorGor1,hsa_mPonAbe1,hsa_mPonPyg2)
    resources:
        partition="icelake",
        time='08:00:00',
        account = "DURBIN-SL2-CPU",
        mem="12G",
        cores=1
    shell:
        "awk -v seq=\"chr{wildcards.chrom}{params.chrom_name_suffix}\" -v RS='>' '$1 == seq {{print RS $0}}' {input.fasta} | "
        "/home/tc557/sdust/sdust -w {wildcards.sdust_w} -t {wildcards.sdust_t} - "
        "> {output.dat}"

rule write_mhs:
    input:
        SNP_pos_summary = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.variantsummary.chr{chrom}.0based.txt.gz',
        indel_bed = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.chr{chrom}.indel.bed.gz',
        bed = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.bed',
        trf_bed = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.chr{chrom}.fasta.trf.bed',
        sdust_bed = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}.pri.chr{chrom}.fasta.sdust.w{sdust_w}.t{sdust_t}dat',
    output: 
        mhsfile = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.260127_4beds.chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.mhs',
        log = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.260127_4beds.chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.log.txt'
    resources:
        partition="icelake",
        time='00:05:00',
        account = "DURBIN-SL2-CPU",
        mem="20G",
        cores=1
    shell:
        r"""
        if [[ "{wildcards.sample}" == "HG002" && "{wildcards.chrom}" == "23" ]]; then
            echo empty > {output.mhsfile}
            echo empty > {output.log}
        else
            /home/tc557/T2T_PSMC/T2Tsequences/251228_make_mhs.py \
                -iv {input.SNP_pos_summary} \
                -negative_beds {input.indel_bed} {input.trf_bed} {input.sdust_bed} \
                -positive_beds {input.bed} \
                -chrom {wildcards.chrom} \
                -truncate 50000 \
                -m {wildcards.m} \
                -o {output.mhsfile} | tee {output.log}
        fi
        """

def get_mhs_files_260127(wildcards):
    if 'HG002' in wildcards.sample:
        zchroms = range(1,23)
    else:
        zchroms = range(1,24)
    return [f'/home/tc557/rds/hpc-work/T2Tsequences/{wildcards.sample}_all_pap.dip.260127_4beds.chr{chrom}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.mhs' for chrom in zchroms]


rule run_MSMC2:
    input:
        mhsfiles = get_mhs_files_260127
    output:
        outfile = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_w}.t{sdust_t}.final.txt',
        outlog = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_w}.t{sdust_t}.log',
    resources:
        partition="icelake-himem",
        time='01:30:00',
        account = "DURBIN-SL2-CPU",
        mem="180G",
        cores=20
    params:
        outfile = lambda wildcards: f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{wildcards.D}/iterations_{wildcards.iterations}/theta_{wildcards.theta}/rhoovermu_{wildcards.rhoovermu}/{wildcards.sample}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}',
        ztheta = lambda wildcards: get_theta_MSMC2(wildcards)
    shell:
        "/home/tc557/msmc2_installation_220407/msmc2-2.1.3/build/release/msmc2 -i {wildcards.iterations} -t {resources.cores} -p {wildcards.D}*1 {params.ztheta} -r {wildcards.rhoovermu} -o {params.outfile} {input.mhsfiles} | tee {output.outlog}"

rule blockbootstrap_mhs:
    input:
        mhsfile = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.260127_4beds.chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.mhs'
    output:
        bbs_mhsfile = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.260127_4beds.bbs5Mb_{bootstrap}.chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.mhs'
    resources:
        partition="icelake",
        time='00:02:00',
        account = "DURBIN-SL2-CPU",
        mem="5G",
        cores=1
    shell:
        'python /home/tc557/cobraa/block_bootstrap_mhs.py -inmhs {input.mhsfile} -windowsize 5000000 -outmhs {output.bbs_mhsfile}'


def get_mhs_files_260127_bootstrap(wildcards):
    if 'HG002' in wildcards.sample:
        zchroms = range(1,23)
    else:
        zchroms = range(1,24)
    return [f'/home/tc557/rds/hpc-work/T2Tsequences/{wildcards.sample}_all_pap.dip.260127_4beds.bbs5Mb_{wildcards.bootstrap}.chr{chrom}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.mhs' for chrom in zchroms]

rule run_MSMC2_bootstrap:
    input:
        mhsfiles = get_mhs_files_260127_bootstrap
    output:
        outfile = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/bootstrap{bootstrap}.{sample}.m{m}.w{sdust_w}.t{sdust_t}.final.txt',
        outlog = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/bootstrap{bootstrap}.{sample}.m{m}.w{sdust_w}.t{sdust_t}.log',
    resources:
        partition="icelake-himem",
        time='01:30:00',
        account = "DURBIN-SL2-CPU",
        mem="180G",
        cores=20
    params:
        outfile = lambda wildcards: f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{wildcards.D}/iterations_{wildcards.iterations}/theta_{wildcards.theta}/rhoovermu_{wildcards.rhoovermu}/bootstrap{wildcards.bootstrap}.{wildcards.sample}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}',
        ztheta = lambda wildcards: get_theta_MSMC2(wildcards)
    shell:
        "/home/tc557/msmc2_installation_220407/msmc2-2.1.3/build/release/msmc2 -i {wildcards.iterations} -t {resources.cores} -p {wildcards.D}*1 {params.ztheta} -r {wildcards.rhoovermu} -o {params.outfile} {input.mhsfiles} | tee {output.outlog}"

def get_params_from_final_file(wildcards):
    filename_log = f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{wildcards.D}/iterations_{wildcards.iterations}/theta_{wildcards.theta}/rhoovermu_{wildcards.rhoovermu}/{wildcards.sample}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.log'
    with open(filename_log,'r') as f: lines = f.readlines()
    theta = float([i for i in lines if 'mutationRate' in i][0].split(' ')[-1])
    rho = float([i for i in lines if 'recombinationRate' in i][0].split(' ')[-1])
    
    filename_finalparams = f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{wildcards.D}/iterations_{wildcards.iterations}/theta_{wildcards.theta}/rhoovermu_{wildcards.rhoovermu}/{wildcards.sample}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.final.txt'
    data = pd.read_csv(filename_finalparams,sep="\t")
    lambda_ = np.array(data['lambda']*theta)
    lambda_str = "".join([f'{str(i)},' for i in lambda_])[0:-1]
    if wildcards.sample == 'mPonAbe1': # error in MSMC2 if I give the natural file - I suspect a numerical issue because some the numbers are very large and some are very small. See https://github.com/stschiff/msmc2/issues/69
        lambda_str = '18.388186752,11.732607744000001,1.9501360319999999,1.4972413247999998,1.1211273216,0.87359196288,1.06168767744,0.95656556736,0.83748919488,0.79523810496,0.7476692620799998,0.65183591616,0.55690469952,0.53003555904,0.49715417856,0.4683153696,0.48458067839999996,0.54814853184,0.64868491584,0.75359234496,0.8151072384,0.81064117056,0.76490591424,0.71553515904,0.69324118656,0.70562346816,0.7457171904,0.80733179904,0.8867707929599999,0.97038257472,1.03194321984,1.05693067584,1.05319663488,1.03961425152,1.037536656,1.06237043328,1.11802089984,1.1949283008,1.2724363392,1.3316906304,1.3702393536,1.4090227007999998,1.4837973695999997,1.6240555968000001,1.8259495488,2.010375744,2.0275854144,1.7775348864,1.4278278143999998,1.3317727488,1.5853074431999998,1.882153728,2.5943783423999998,4.1076796799999995,5,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,10'
    return theta, rho, lambda_str

rule decode_MSMC2:
    input:
        # inferencefile = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_w}.t{sdust_t}.final.txt',
        # logfile = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_w}.t{sdust_t}.final.txt',
        # mhsfile = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.260127_4beds.chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.mhs'
    output:
        inferencefile = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_w}.t{sdust_t}.decoding.chr{chrom}.txt.gz'
    params:
        # mutrate = get_mutation_rate_from_final_file(input.inferencefile),
        # recombrate = get_recombination_rate_from_final_file(input.inferencefile),
        # lambdavec = get_lambda_from_final_file(input.inferencefile),
        mutrate_recombrate_lambdavec = lambda wildcards: get_params_from_final_file(wildcards),
        mhsfile = lambda wildcards: f'/home/tc557/rds/hpc-work/T2Tsequences/{wildcards.sample}_all_pap.dip.260127_4beds.chr{wildcards.chrom}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.mhs'
    resources:
        partition="icelake-himem",
        time='00:10:00',
        account = "DURBIN-SL2-CPU",
        mem="30G",
        cores=1
    shell:
        "/home/tc557/msmc2_installation_220407/msmc2-2.1.3/build/decode -m {params.mutrate_recombrate_lambdavec[0]} -r {params.mutrate_recombrate_lambdavec[1]} -l {params.mutrate_recombrate_lambdavec[2]} {params.mhsfile} | gzip -c > {output.inferencefile}"

def get_params_from_final_bootstrap_file(wildcards):
    filename_log = f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{wildcards.D}/iterations_{wildcards.iterations}/theta_{wildcards.theta}/rhoovermu_{wildcards.rhoovermu}/bootstrap{wildcards.bootstrap}.{wildcards.sample}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.log'
    
    with open(filename_log,'r') as f: lines = f.readlines()
    theta = float([i for i in lines if 'mutationRate' in i][0].split(' ')[-1])
    rho = float([i for i in lines if 'recombinationRate' in i][0].split(' ')[-1])
    
    filename_finalparams = f'/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{wildcards.D}/iterations_{wildcards.iterations}/theta_{wildcards.theta}/rhoovermu_{wildcards.rhoovermu}/bootstrap{wildcards.bootstrap}.{wildcards.sample}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.final.txt'
    data = pd.read_csv(filename_finalparams,sep="\t")
    lambda_ = np.array(data['lambda']*theta)
    lambda_str = "".join([f'{str(i)},' for i in lambda_])[0:-1]
    # if wildcards.sample == 'mPonAbe1': # error in MSMC2 if I give the natural file - I suspect a numerical issue because some the numbers are very large and some are very small. See https://github.com/stschiff/msmc2/issues/69
    #     lambda_str = '18.388186752,11.732607744000001,1.9501360319999999,1.4972413247999998,1.1211273216,0.87359196288,1.06168767744,0.95656556736,0.83748919488,0.79523810496,0.7476692620799998,0.65183591616,0.55690469952,0.53003555904,0.49715417856,0.4683153696,0.48458067839999996,0.54814853184,0.64868491584,0.75359234496,0.8151072384,0.81064117056,0.76490591424,0.71553515904,0.69324118656,0.70562346816,0.7457171904,0.80733179904,0.8867707929599999,0.97038257472,1.03194321984,1.05693067584,1.05319663488,1.03961425152,1.037536656,1.06237043328,1.11802089984,1.1949283008,1.2724363392,1.3316906304,1.3702393536,1.4090227007999998,1.4837973695999997,1.6240555968000001,1.8259495488,2.010375744,2.0275854144,1.7775348864,1.4278278143999998,1.3317727488,1.5853074431999998,1.882153728,2.5943783423999998,4.1076796799999995,5,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,10'
    return theta, rho, lambda_str

rule decode_MSMC2_bootstrap:
    input:
        inferencefile = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/bootstrap{bootstrap}.{sample}.m{m}.w{sdust_w}.t{sdust_t}.final.txt',
        # logfile = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/bootstrap{bootstrap}.{sample}.m{m}.w{sdust_w}.t{sdust_t}.final.txt',
        mhsfile = '/home/tc557/rds/hpc-work/T2Tsequences/{sample}_all_pap.dip.260127_4beds.bbs5Mb_{bootstrap}.chr{chrom}.m{m}.w{sdust_w}.t{sdust_t}.mhs'
    output:
        inferencefile = '/home/tc557/rds/hpc-work/T2Tsequences/MSCM2inference260127_4beds/D_{D}/iterations_{iterations}/theta_{theta}/rhoovermu_{rhoovermu}/{sample}.m{m}.w{sdust_w}.t{sdust_t}.decoding.chr{chrom}.bootstrap{bootstrap}.txt.gz'
    params:
        # mutrate = get_mutation_rate_from_final_file(input.inferencefile),
        # recombrate = get_recombination_rate_from_final_file(input.inferencefile),
        # lambdavec = get_lambda_from_final_file(input.inferencefile),
        mutrate_recombrate_lambdavec = lambda wildcards: get_params_from_final_bootstrap_file(wildcards),
        mhsfile = lambda wildcards: f'/home/tc557/rds/hpc-work/T2Tsequences/bootstrap{wildcards.bootstrap}.{wildcards.sample}_all_pap.dip.260127_4beds.chr{wildcards.chrom}.m{wildcards.m}.w{wildcards.sdust_w}.t{wildcards.sdust_t}.mhs'
    resources:
        partition="icelake-himem",
        time='00:10:00',
        account = "DURBIN-SL2-CPU",
        mem="30G",
        cores=1
    shell:
        "/home/tc557/msmc2_installation_220407/msmc2-2.1.3/build/decode -m {params.mutrate_recombrate_lambdavec[0]} -r {params.mutrate_recombrate_lambdavec[1]} -l {params.mutrate_recombrate_lambdavec[2]} -s 10000 {input.mhsfile} | gzip -c > {output.inferencefile}"
