import pdb
import os
import pandas as pd
import numpy as np


base_path = '/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/greatape'
files_list = os.listdir(f'{base_path}/vcfs/')
species = [file.split('.vcf')[0] for file in files_list if file.endswith('.vcf.gz')]

chroms = [i for i in range(1,23)]

def get_mhs_files(wildcards):
    mhsfiles = [f'/home/tc557/rds/rds-durbin-group-8b3VcZwY7rY/projects/greatape/230131/{wildcards.species}/{wildcards.individual}/mhs/chr{zchromz}.mhs'for zchromz in chroms]
    return mhsfiles

def debug_me(wildcards):
    pdb.set_trace()
    return None

theta_pergenus_230426 = {}
theta_pergenus_230426['Gorilla'] = 0.0014
theta_pergenus_230426['Pan_paniscus'] = 0.001
theta_pergenus_230426['Pan_troglodytes'] = 0.001
theta_pergenus_230426['Pongo_pygmaeus'] = 0.0014
theta_pergenus_230426['Pongo_abelii']= 0.0014


species_individuals = [('Pongo_pygmaeus', 'Pongo_pygmaeus-A939_Nonja'), ('Pongo_pygmaeus', 'Pongo_pygmaeus-A940_Temmy'), ('Pongo_pygmaeus', 'Pongo_pygmaeus-A941_Sari'), ('Pongo_pygmaeus', 'Pongo_pygmaeus-A943_Tilda'), ('Pongo_pygmaeus', 'Pongo_pygmaeus-A944_Napoleon'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Akwaya_Jean'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Banyo'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Basho'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Damian'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Julie'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Kopongo'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Koto'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Paquita'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Taweh'), ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Tobi'), ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-100037_Vincent'), ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-100040_Andromeda'), ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-9729_Harriet'), ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-A910_Bwambale'), ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-A911_Kidongo'), ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-A912_Nakuu'), ('Pan_troglodytes', 'Pan_troglodytes_troglodytes-A957_Vaillant'), ('Pan_troglodytes', 'Pan_troglodytes_troglodytes-A958_Doris'), ('Pan_troglodytes', 'Pan_troglodytes_troglodytes-A959_Julie'), ('Pan_troglodytes', 'Pan_troglodytes_troglodytes-A960_Clara'), ('Pan_troglodytes', 'Pan_troglodytes_verus-9668_Bosco'), ('Pan_troglodytes', 'Pan_troglodytes_verus-9730_Donald'), ('Pan_troglodytes', 'Pan_troglodytes_verus-A956_Jimmie'), ('Pan_troglodytes', 'Pan_troglodytes_verus-Clint'), ('Pan_troglodytes', 'Pan_troglodytes_verus-X00100_Koby'), ('Gorilla', 'Gorilla_beringei_graueri-9732_Mkubwa'), ('Gorilla', 'Gorilla_beringei_graueri-A929_Kaisi'), ('Gorilla', 'Gorilla_beringei_graueri-Victoria'), ('Gorilla', 'Gorilla_gorilla_dielhi-B646_Nyango'), ('Gorilla', 'Gorilla_gorilla_gorilla-9749_Kowali'), ('Gorilla', 'Gorilla_gorilla_gorilla-9750_Azizi'), ('Gorilla', 'Gorilla_gorilla_gorilla-9751_Bulera'), ('Gorilla', 'Gorilla_gorilla_gorilla-9752_Suzie'), ('Gorilla', 'Gorilla_gorilla_gorilla-9753_Kokomo'), ('Gorilla', 'Gorilla_gorilla_gorilla-A930_Sandra'), ('Gorilla', 'Gorilla_gorilla_gorilla-A931_Banjo'), ('Gorilla', 'Gorilla_gorilla_gorilla-A932_Mimi'), ('Gorilla', 'Gorilla_gorilla_gorilla-A933_Dian'), ('Gorilla', 'Gorilla_gorilla_gorilla-A934_Delphi'), ('Gorilla', 'Gorilla_gorilla_gorilla-A936_Coco'), ('Gorilla', 'Gorilla_gorilla_gorilla-A937_Kolo'), ('Gorilla', 'Gorilla_gorilla_gorilla-A962_Amani'), ('Gorilla', 'Gorilla_gorilla_gorilla-B642_Akiba_Beri'), ('Gorilla', 'Gorilla_gorilla_gorilla-B643_Choomba'), ('Gorilla', 'Gorilla_gorilla_gorilla-B644_Paki'), ('Gorilla', 'Gorilla_gorilla_gorilla-B647_Anthal'), ('Gorilla', 'Gorilla_gorilla_gorilla-B650_Katie'), ('Gorilla', 'Gorilla_gorilla_gorilla-KB3782_Vila'), ('Gorilla', 'Gorilla_gorilla_gorilla-KB3784_Dolly'), ('Gorilla', 'Gorilla_gorilla_gorilla-KB4986_Katie'), ('Gorilla', 'Gorilla_gorilla_gorilla-KB5792_Carolyn'), ('Gorilla', 'Gorilla_gorilla_gorilla-KB5852_Helen'), ('Gorilla', 'Gorilla_gorilla_gorilla-KB6039_Oko'), ('Gorilla', 'Gorilla_gorilla_gorilla-KB7973_Porta'), ('Gorilla', 'Gorilla_gorilla_gorilla-X00108_Abe'), ('Gorilla', 'Gorilla_gorilla_gorilla-X00109_Tzambo'), ('Pongo_abelii', 'Pongo_abelii-A947_Elsi'), ('Pongo_abelii', 'Pongo_abelii-A948_Kiki'), ('Pongo_abelii', 'Pongo_abelii-A949_Dunja'), ('Pongo_abelii', 'Pongo_abelii-A950_Babu'), ('Pongo_abelii', 'Pongo_abelii-A952_Buschi'), ('Pan_paniscus', 'Pan_paniscus-9731_LB502'), ('Pan_paniscus', 'Pan_paniscus-A914_Hortense'), ('Pan_paniscus', 'Pan_paniscus-A915_Kosana'), ('Pan_paniscus', 'Pan_paniscus-A917_Dzeeta'), ('Pan_paniscus', 'Pan_paniscus-A918_Hermien'), ('Pan_paniscus', 'Pan_paniscus-A919_Desmond'), ('Pan_paniscus', 'Pan_paniscus-A922_Catherine'), ('Pan_paniscus', 'Pan_paniscus-A923_Kombote'), ('Pan_paniscus', 'Pan_paniscus-A924_Chipita'), ('Pan_paniscus', 'Pan_paniscus-A925_Bono'), ('Pan_paniscus', 'Pan_paniscus-A926_Natalie'), ('Pan_paniscus', 'Pan_paniscus-A927_Salonga'), ('Pan_paniscus', 'Pan_paniscus-A928_Kumbuka')]
# some of these have obvious data quality issues, so I removed them

species_individuals_nonbad = [
    ('Pongo_pygmaeus', 'Pongo_pygmaeus-A939_Nonja'), 
    ('Pongo_pygmaeus', 'Pongo_pygmaeus-A940_Temmy'), 
    ('Pongo_pygmaeus', 'Pongo_pygmaeus-A941_Sari'), 
    ('Pongo_pygmaeus', 'Pongo_pygmaeus-A943_Tilda'), 
    ('Pongo_pygmaeus', 'Pongo_pygmaeus-A944_Napoleon'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Akwaya_Jean'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Banyo'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Basho'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Damian'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Julie'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Kopongo'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Koto'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Paquita'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Taweh'), 
    ('Pan_troglodytes', 'Pan_troglodytes_ellioti-Tobi'), 
    ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-100037_Vincent'), 
    ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-100040_Andromeda'), 
    ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-9729_Harriet'), 
    ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-A910_Bwambale'), 
    ('Pan_troglodytes', 'Pan_troglodytes_schweinfurthii-A912_Nakuu'), 
    ('Pan_troglodytes', 'Pan_troglodytes_troglodytes-A957_Vaillant'), 
    ('Pan_troglodytes', 'Pan_troglodytes_troglodytes-A958_Doris'), 
    ('Pan_troglodytes', 'Pan_troglodytes_verus-9730_Donald'), 
    ('Pan_troglodytes', 'Pan_troglodytes_verus-A956_Jimmie'), 
    ('Pan_troglodytes', 'Pan_troglodytes_verus-Clint'), 
    ('Pan_troglodytes', 'Pan_troglodytes_verus-X00100_Koby'), 
    ('Gorilla', 'Gorilla_beringei_graueri-9732_Mkubwa'), 
    ('Gorilla', 'Gorilla_beringei_graueri-A929_Kaisi'), 
    ('Gorilla', 'Gorilla_beringei_graueri-Victoria'), 
    ('Gorilla', 'Gorilla_gorilla_dielhi-B646_Nyango'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-9750_Azizi'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-9751_Bulera'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-9752_Suzie'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-A930_Sandra'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-A931_Banjo'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-A933_Dian'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-A934_Delphi'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-A936_Coco'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-A937_Kolo'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-A962_Amani'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-B642_Akiba_Beri'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-B643_Choomba'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-B644_Paki'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-B647_Anthal'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-B650_Katie'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-KB3782_Vila'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-KB3784_Dolly'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-KB4986_Katie'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-KB5792_Carolyn'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-KB5852_Helen'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-KB7973_Porta'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-X00108_Abe'), 
    ('Gorilla', 'Gorilla_gorilla_gorilla-X00109_Tzambo'), 
    ('Pongo_abelii', 'Pongo_abelii-A947_Elsi'), 
    ('Pongo_abelii', 'Pongo_abelii-A948_Kiki'), 
    ('Pongo_abelii', 'Pongo_abelii-A949_Dunja'), 
    ('Pongo_abelii', 'Pongo_abelii-A950_Babu'), 
    ('Pongo_abelii', 'Pongo_abelii-A952_Buschi'), 
    ('Pan_paniscus', 'Pan_paniscus-A914_Hortense'), 
    ('Pan_paniscus', 'Pan_paniscus-A915_Kosana'), 
    ('Pan_paniscus', 'Pan_paniscus-A917_Dzeeta'), 
    ('Pan_paniscus', 'Pan_paniscus-A918_Hermien'), 
    ('Pan_paniscus', 'Pan_paniscus-A919_Desmond'), 
    ('Pan_paniscus', 'Pan_paniscus-A922_Catherine'), 
    ('Pan_paniscus', 'Pan_paniscus-A923_Kombote'), 
    ('Pan_paniscus', 'Pan_paniscus-A924_Chipita'), 
    ('Pan_paniscus', 'Pan_paniscus-A925_Bono'), 
    ('Pan_paniscus', 'Pan_paniscus-A926_Natalie'), 
    ('Pan_paniscus', 'Pan_paniscus-A928_Kumbuka')]

spec_ind_nonbad = {}
for spe in species:
    spec_ind_nonbad[spe] = [i[1] for i in species_individuals_nonbad if i[0]==spe]

Ds = [64]
bs = [100]
spread1s = [0.005,0.01,0.05,0.075,0.1]
spread2s = [30,50,70,100,150] 
iterationss = [30]
threshs = [1]
muoverrs = [1.5] # starting mu over r ratio
chroms = range(1,23)
bootstraps = range(1,31)

rule all:
    input:
        [f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813_ape/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/{species}/{individual}/final_parameters.txt' for D in Ds for b in bs for spread1 in spread1s for spread2 in spread2s for muoverr in muoverrs for iterations in iterationss for thresh in threshs for (species,individual) in species_individuals_nonbad]

rule run_cobraa:
    input:
        mhsfiles = get_mhs_files,
    output:
        outfile = '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813_ape/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/{species}/{individual}/final_parameters.txt'
    log:
        '/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813_ape/D_{D}/b_{b}/spread1_{spread1}/spread2_{spread2}/muoverr_{muoverr}/iterations{iterations}/thresh_{thresh}/{species}/{individual}/log.txt'
    resources:
        partition="icelake-himem",
        time='02:00:00',
        account = "DURBIN-SL2-CPU",
        mem="180G",
        cores=20
    params:
        outfile = lambda wildcards: f'/home/tc557/rds/hpc-work/PSMCplus_analysis_240813/240813_ape/D_{wildcards.D}/b_{wildcards.b}/spread1_{wildcards.spread1}/spread2_{wildcards.spread2}/muoverr_{wildcards.muoverr}/iterations{wildcards.iterations}/thresh_{wildcards.thresh}/{wildcards.species}/{wildcards.individual}/',
        ztheta = lambda wildcards: theta_pergenus_230426[wildcards.species]
    shell:
        # 'python /home/tc557/PSMCplus/PSMCplus.py -in {input.mhsfiles} -o {params.outfile} -D {wildcards.D} -b {wildcards.b} -mu_over_rho_ratio {wildcards.muoverr} {params.ztheta} -its {wildcards.iterations} -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2} | tee {log}'
        'python /home/tc557/cobraa/cobraa.py -in {input.mhsfiles} -o {params.outfile} -D {wildcards.D} -b {wildcards.b} -mu_over_rho_ratio {wildcards.muoverr} -theta {params.ztheta} -its {wildcards.iterations} -thresh {wildcards.thresh} -spread_1 {wildcards.spread1} -spread_2 {wildcards.spread2}  | tee {log}'
