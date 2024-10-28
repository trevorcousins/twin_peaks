import pdb
import numpy as np 
import pandas as pd

def split_popsam(wildcards):
    zpop = wildcards.popsam.split('_')[0]
    zsam = wildcards.popsam.split('_')[1]
    return zpop, zsam

def get_mhs_files_allchroms(wildcards):
    zpop, zsam = split_popsam(wildcards)
    return [f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/230213/{zpop}/{zsam}/mhs/chr{zchromz}.mhs' for zchromz in range(1,23)]


def get_mhs_files_allchroms_filter(wildcards):
    zpop, zsam = split_popsam(wildcards)
    return [f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/230213/{zpop}/{zsam}/mhs_{wildcards.filter}/chr{zchromz}.mhs' for zchromz in range(1,23)]

def get_mhs_files_allchroms_bbs(wildcards):
    zpop, zsam = split_popsam(wildcards)
    return [f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/230213/{zpop}/{zsam}/mhs_bootstrap240814/chr{zchromz}_bootstrap_ws5Mb_{wildcards.bbs}.mhs' for zchromz in range(1,23)]

def get_params_from_ouput_file_240813(wildcards):
    # panmictic 
    final_params_file = get_ouput_file_240106(wildcards)
    with open(final_params_file) as f:
        finallines = f.readlines()
    ztheta = float([i for i in finallines if 'theta' in i ][0].split(' ')[-1])
    zrho = float([i for i in finallines if 'rho' in i ][0].split(' ')[-1])
    
    final_params = np.loadtxt(final_params_file)
    lambdaA_parameters = final_params[:,2]*ztheta/4

    logfile = final_params_file.replace('final_parameters.txt','log.txt')
    with open(logfile) as f:
        loglines = f.readlines()
    zlambdaA = [i for i in loglines if 'lambda_A updated' in i][-1].split(' ')[2].split('\n')[0].split('[')[1].split(']')[0]
    zlambdaA_array = np.array([float(i) for i in zlambdaA.split(',')])
    if np.max(np.abs(zlambdaA_array - lambdaA_parameters))>1e-15:
        print(f'PROBLEM; log lambdaA does not seem to be the same as final lambdaA')
    return ztheta, zrho, zlambdaA

pops = ['ACB','ASW','BEB','CDX','CEU','CHB','CHS','CLM','ESN','FIN','GBR','GIH','GWD','IBS','ITU','JPT','KHV','LWK','MSL','MXL','PEL','PJL','PUR','STU','TSI','YRI']
pops_africans = ['ESN','YRI','MSL','GWD','ACB','ASW','LWK']
pops=pops_africans

numsamplesperpop = 1
pops_samples = {}
names_IDs = {}
for pop in pops:
    file_to_samples = f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/230213/231103_samples_name_ID_{pop}.txt'

    with open(file_to_samples,'r') as f:
        lines = f.readlines()
    pops_samples[pop] = [(line.split(' ')[0],line.split(' ')[1]) for line in lines]
    for i in pops_samples[pop]:
        names_IDs[i[0]] = i[1]
pop_and_sample = []
zpop_and_sample = []
for zpop in pops:
    for i in range(0,numsamplesperpop):
        pop_and_sample.append(f'{zpop}_{pops_samples[zpop][i][0]}')
        zpop_and_sample.append([zpop,pops_samples[zpop][i][0]])

Ds = [64]
bs = [100]
spread1s = [0.005,0.01,0.075]
spread2s = [50,100,150] # ,100,150]
# spread1s = [0.075]
# spread1s = [0.01]
spread2s = [50,100,150] 
# spread2=200 breaks the code
iterationss = [30]
threshs = [1]
muoverrs = [1.5] # starting mu over r ratio
fixedthetas = [0.0008]
chroms = range(1,23)
bootstraps = range(1,31)
filters = ['SRM','SGD','CPG']

rule all:
    input:
        # [f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/final_parameters.txt' for theta in fixedthetas for D in Ds for b in bs for spread1 in spread1s for spread2 in spread2s for muoverr in muoverrs for iterations in iterationss for thresh in threshs for popsam in pop_and_sample ] # standard PSMC
        # [f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/chr{chrom}_decode.txt.gz' for chrom in chroms for theta in fixedthetas for D in Ds for b in bs for spread1 in spread1s for spread2 in spread2s for muoverr in muoverrs for iterations in iterationss for thresh in threshs for popsam in pop_and_sample ] # PSMC decoding
        # [f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/230213/{popsam[0]}/{popsam[1]}/mhs_bootstrap240814/chr{chrom}_bootstrap_ws5Mb_{bootstrap}.mhs' for chrom in chroms for bootstrap in bootstraps for popsam in zpop_and_sample]
        # [f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/bbs{bbs}_final_parameters.txt' for bbs in bootstraps for theta in fixedthetas for D in Ds for b in bs for spread1 in spread1s for spread2 in spread2s for muoverr in muoverrs for iterations in iterationss for thresh in threshs for popsam in pop_and_sample ] # standard PSMC with block bootstrap
        [f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/filter{filter}_final_parameters.txt' for filter in filters for theta in fixedthetas for D in Ds for b in bs for spread1 in spread1s for spread2 in spread2s for muoverr in muoverrs for iterations in iterationss for thresh in threshs for popsam in pop_and_sample ] # standard PSMC with filters (SRM, SGD, CpGs)

rule run_cobraa:
    input:
        mhsfiles = get_mhs_files_allchroms,
    output:
        outfile = '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/final_parameters.txt'
    log:
        '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/log.txt'
    resources:
        partition="icelake-himem",
        time='02:00:00',
        account = "DURBIN-SL2-CPU",
        mem="180G",
        cores=20
    params:
        outfile = lambda wildcards: f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{wildcards.D}/b_{wildcards.b}/spread1_{wildcards.spread1}/spread2_{wildcards.spread2}/muoverr_{wildcards.muoverr}/iterations{wildcards.iterations}/thresh_{wildcards.thresh}/thetafixed_{wildcards.theta}/popsample_{wildcards.popsam}/',
        ztheta = '' # default theta (inferred from data)
    shell:
        # 'python /home/tc557/PSMCplus/PSMCplus.py -in {input.mhsfiles} -o {params.outfile} -D {wildcards.D} -b {wildcards.b} -mu_over_rho_ratio {wildcards.muoverr} {params.ztheta} -its {wildcards.iterations} -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2} | tee {log}'
        'python /home/tc557/cobraa/cobraa.py -in {input.mhsfiles} -o {params.outfile} -D {wildcards.D} -b {wildcards.b} -mu_over_rho_ratio {wildcards.muoverr} {params.ztheta} -its {wildcards.iterations} -theta {wildcards.theta} -thresh {wildcards.thresh} -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2}  | tee {log}'


# run decode, with the inferred parameters 
rule decode_cobraa:
    input:
        mhsfile = get_mhs_file_single,
        inferred_params = '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/final_parameters.txt'
    output:
        decode_file = '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/chr{chrom}_decode.txt.gz'
    log:
        '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/chr{chrom}_decode_log.txt'
    resources:
        partition="icelake-himem",
        time='00:05:00',
        account = "DURBIN-SL2-CPU",
        mem="10G",
        cores=1
    params:
        zts = 13, # composite ML pair
        zte = 21, # composite ML pair
        theta_rho_lambdaA = lambda wildcards: get_params_from_ouput_file_240813(wildcards),
        zD = 32,
        zspread_1 = lambda wildcards: f'{wildcards.spread1}',
        zspread_2 = lambda wildcards: f'{wildcards.spread2}'
    shell:
        'python /home/tc557/cobraa/cobraa.py -in {input.mhsfile} -o {output.decode_file} -D {wildcards.D} -b {wildcards.b} -theta {params.theta_rho_lambdaA[0]} -rho {params.theta_rho_lambdaA[1]} -its 1 -thresh 1 -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2} -lambda_A_fg {params.theta_rho_lambdaA[2]} -decode -decode_downsample 10 | tee {log}'

rule blockbootstrap_mhs:
    input:
        mhsfile = '/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/230213/{population}/{sample}/mhs/chr{chrom}.mhs'
    output:
        bbs_mhsfile = '/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/human/1000Genomes_30X/230213/{population}/{sample}/mhs_bootstrap240814/chr{chrom}_bootstrap_ws5Mb_{bootstrap}.mhs'
    resources:
        partition="icelake",
        time='00:02:00',
        account = "DURBIN-SL2-CPU",
        mem="5G",
        cores=1
    shell:
        'python /home/tc557/cobraa/block_bootstrap_mhs.py -inmhs {input.mhsfile} -windowsize 5000000 -outmhs {output.bbs_mhsfile}'

rule run_cobraa_bbs:
    input:
        mhsfiles = get_mhs_files_allchroms_bbs,
    output:
        outfile = '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/bbs{bbs}_final_parameters.txt'
    log:
        '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/bbs{bbs}_log.txt'
    resources:
        partition="icelake-himem",
        time='02:00:00',
        account = "DURBIN-SL2-CPU",
        mem="180G",
        cores=20
    params:
        outfile = lambda wildcards: f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{wildcards.D}/b_{wildcards.b}/spread1_{wildcards.spread1}/spread2_{wildcards.spread2}/muoverr_{wildcards.muoverr}/iterations{wildcards.iterations}/thresh_{wildcards.thresh}/thetafixed_{wildcards.theta}/popsample_{wildcards.popsam}/bbs{wildcards.bbs}_',
        ztheta = '' # default theta (inferred from data)
    shell:
        # 'python /home/tc557/PSMCplus/PSMCplus.py -in {input.mhsfiles} -o {params.outfile} -D {wildcards.D} -b {wildcards.b} -mu_over_rho_ratio {wildcards.muoverr} {params.ztheta} -its {wildcards.iterations} -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2} | tee {log}'
        'python /home/tc557/cobraa/cobraa.py -in {input.mhsfiles} -o {params.outfile} -D {wildcards.D} -b {wildcards.b} -mu_over_rho_ratio {wildcards.muoverr} {params.ztheta} -its {wildcards.iterations} -theta {wildcards.theta} -thresh {wildcards.thresh} -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2}  | tee {log}'

rule run_cobraa_filter:
    input:
        mhsfiles = get_mhs_files_allchroms_filter,
    output:
        outfile = '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/filter{filter}_final_parameters.txt'
    log:
        '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/thetafixed_{theta}/popsample_{popsam}/filter{filter}_log.txt'
    resources:
        partition="icelake-himem",
        time='02:00:00',
        account = "DURBIN-SL2-CPU",
        mem="180G",
        cores=20
    params:
        outfile = lambda wildcards: f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813/D_{wildcards.D}/b_{wildcards.b}/spread1_{wildcards.spread1}/spread2_{wildcards.spread2}/muoverr_{wildcards.muoverr}/iterations{wildcards.iterations}/thresh_{wildcards.thresh}/thetafixed_{wildcards.theta}/popsample_{wildcards.popsam}/filter{wildcards.filter}_',
        ztheta = '' # default theta (inferred from data)
    shell:
        # 'python /home/tc557/PSMCplus/PSMCplus.py -in {input.mhsfiles} -o {params.outfile} -D {wildcards.D} -b {wildcards.b} -mu_over_rho_ratio {wildcards.muoverr} {params.ztheta} -its {wildcards.iterations} -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2} | tee {log}'
        'python /home/tc557/cobraa/cobraa.py -in {input.mhsfiles} -o {params.outfile} -D {wildcards.D} -b {wildcards.b} -mu_over_rho_ratio {wildcards.muoverr} {params.ztheta} -its {wildcards.iterations} -theta {wildcards.theta} -thresh {wildcards.thresh} -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2}  | tee {log}'

