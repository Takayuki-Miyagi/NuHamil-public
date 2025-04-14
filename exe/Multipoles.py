#!/usr/bin/env python3
#-*- coding:utf-8 -*-
# script to submit the job
# execute ./NuHamil.py
import os
import sys
import subprocess
from collections import OrderedDict
import time

exe_file = 'NuHamil.exe'
n_nodes=1               # Basically, you can use mpi, but if there is no NNrel file, it will stop. I recommend you to try without mpi for a first run.
"""
 This is for cedar
"""
account="rrg-holt"      # Do not change
n_omp_threads=1            # threads per node
mem = "125G"            # memory per node
walltime = "0-01:00"    # request time
HOME=os.path.expanduser("~")
mpi = False
if(n_nodes > 1): mpi=True

scheduler = 'terminal'
#scheduler = 'local'
#scheduler = 'slurm'
#scheduler = 'pbs'

if(scheduler == 'terminal' or scheduler == 'local'):
    header = ""
elif(scheduler == 'pbs'):
    header  = "#!/bin/bash \n"
    header += "#PBS -q oak \n"
    header += "#PBS -j oe \n"
    header += "#PBS -l mem=250gb \n"
    header += "#PBS -l walltime=72:00:00 \n"
    if(not mpi): header += f"#PBS -l nodes=1:ppn={n_omp_threads}\n"
    if(mpi): header += f"#PBS -l nodes={n_nodes}:ppn={n_omp_threads}\n"
    header += "cd $PBS_O_WORKDIR\n"
elif(scheduler == 'slurm'):
    header =  "#!/bin/bash\n"
    header += f"#SBATCH --account={account}\n"
    header += f"#SBATCH --nodes={n_nodes}\n"
    header += f"#SBATCH --ntasks-per-node=1\n"
    header += f"#SBATCH --cpus-per-task={n_omp_threads}\n"
    header += f"#SBATCH --mem={mem}\n"
    header += f"#SBATCH --time={walltime}\n\n"

operators = "hamil,AxialV_Tz1-N2LO"
operators = "GamowTeller"
operators = "AxialV_Tz1-N2LO"
ext_2bme = ".me2j.gz" # For two-body interaction (gzip)
ext_2bme = ".snt" # For two-body interaction (gzip)
#ext_2bme = ".op.me2j.gz" # For two-body operator (gzip)

# List of available two-body interactions
NNs = [ "N3LO_EM500",\
        "LO_EMN500", "NLO_EMN500", "N2LO_EMN500", "N3LO_EMN500", "N4LO_EMN500", \
        "LQCD1", "LQCD2", "LQCD3", "LQCD4", "LQCD5",\
        "Idaho", "CDBonn", "AV18", "N2LO_opt", "N2LO_sat",\
        "N3LO_EGM_450_500", "N3LO_EGM_600_600", "N3LO_EGM_550_600", "N3LO_EGM_450_700", "N3LO_EGM_600_700",\
        "DN2LOGO394","DN2LOGO450","DNLOGO450","N3LO_EM500_PWA","DN2LO394",\
        "DNLO394"]
NNF = 'N3LO_EM500'
operators = ""
operators = f"Tmag_2B_J1_Tz0_Q100-NLO"

def set_input(params, hw=20, emax=14, e2max=28):
    params['rank'] = 2
    params['hw'] = hw
    params['renorm'] = 'bare'
    params['Operators'] = operators
    params['emax'] = emax
    params['e2max'] = e2max
    params['NNInt'] = NNF

    params['pmax2'] = 8
    params['NMeshMultiPole'] = 20
    params['NMesh2'] = 100
    params['jmax2'] = 8
    params['J2maxLab'] = 13
    params['Lcm2Max'] = 13
    set_file_name_nn(params)

def get_script_name(params):
    fbase = f"vHamil_MECalc_A{params['rank']}_hw{params['hw']}"
    ms = ""
    fbase += ms
    fbase = f"vHamil_MECalc_A2_{params['file_name_nn'].split(ext_2bme)[0]}"
    return fbase

def gen_script(params, batch):
    exe = exe_file
    if(scheduler == 'slurm'): exe = "srun " + exe_file
    fbase = get_script_name(params)
    file_input = "Input_" + fbase + ".dat"
    file_log   = "log_" + fbase + ".dat"
    fsh = "run_" + fbase + ".sh"
    prt = header
    if(scheduler == 'local'): prt += f'export OMP_NUM_THREADS={n_omp_threads}\n'
    prt += 'echo "start ' +fsh + '..."\n'
    prt += "cat > "+file_input + " <<EOF\n"
    prt += "&input\n"
    for key, value in params.items():
        if(isinstance(value, str)):
            prt += f'  {key} = "{value}"\n'
            continue
        if(isinstance(value, list)):
            prt += f"  {key} = "
            for x in value[:-1]:
                prt += f"{x},"
            prt += f"{value[-1]}\n"
            continue
        prt += f"  {key} = {value}\n"
    prt += "&end\n"
    prt += "EOF\n"
    if(scheduler == 'terminal'):
        prt += exe + " " + file_input + "\n"
        prt += "rm " + file_input + "\n"
    else:
        prt += exe + " " + file_input + " > " + file_log + " 2>&1\n"
        prt += "rm " + file_input + "\n"
    f = open(fsh, "w")
    f.write(prt)
    f.close()
    os.chmod(fsh, 0o755)
    return fsh

def set_file_name_nn(params):
    try:
        f = f'TwBME-HO_NN-only_{params["NNInt"]}_{params["renorm"]}'
        f = f'{params["renorm"]}'
        if(params["renorm"] != 'Bare' and params["renorm"] != 'bare'): f += str(params["lambda"])
        f += f'_hw{params["hw"]}_emax{params["emax"]}_e2max{params["e2max"]}'
        if("coul" in params and not params["coul"]): f += '_wocoul'
        if("spin_tensor_decomposition" in params):
            if(params['spin_tensor_decomposition']==0): f += '_C'
            if(params['spin_tensor_decomposition']==1): f += '_LS'
            if(params['spin_tensor_decomposition']==2): f += '_T'
            if(params['spin_tensor_decomposition']==3): f += '_C+LS'
            if(params['spin_tensor_decomposition']==4): f += '_C+T'
        f += ext_2bme
    except:
        f = "default"
    params["file_name_nn"] = f

def main():
    params = OrderedDict()
    for hw in [20]:
        set_input(params,hw,emax=2, e2max=4)
        #set_input(params,hw,emax=0,e2max=8)
        fsh = gen_script(params, scheduler)
        if(scheduler == 'terminal'): cmd = "./" + fsh
        elif(scheduler == 'local'): cmd = "nohup ./" + fsh + " &"
        elif(scheduler == "pbs"): cmd = "qsub " + fsh
        elif(scheduler == 'slurm'): cmd = "sbatch " + fsh

        subprocess.call(cmd,shell=True)
        time.sleep(1)


if(__name__ == "__main__"):
    main()
