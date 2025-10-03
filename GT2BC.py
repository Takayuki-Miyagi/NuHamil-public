#!/usr/bin/env python
#-*- coding:utf-8 -*-
# script to submit the job
# execute ./NuHamil.py
# Normalization is different from usual GamowTeller operator, sigma tau
# To combine with sigma tau operator, you need to multiply -gA / sqrt(2) to the result matrix elements
import os
import sys
import subprocess
from collections import OrderedDict
import time

HOME=os.path.expanduser("~")
path_to_nninput = f"{HOME}/NuHamil-public/input_nn_files" # this is default; please modify here if you did not download the code in your home directory.
exe_file = 'NuHamil.exe'
n_nodes=1               # Basically, you can use mpi, but if there is no NNrel file, it will stop. I recommend you to try without mpi for a first run.
n_omp_threads=48        # threads per node
mem = "125G"            # memory per node
walltime = "0-01:00"    # request time
mpi = False             # Do NOT change
batch = 'terminal'
#batch = 'local'
#batch = 'slurm'
#batch = 'pbs'

if(batch == 'terminal' or batch == 'local'):
    header = ""
elif(batch == 'pbs'):
    header  = "#!/bin/bash \n"
    header += "#PBS -q oak \n"
    header += "#PBS -j oe \n"
    header += "#PBS -l mem=250gb \n"
    header += "#PBS -l walltime=72:00:00 \n"
    if(not mpi): header += "#PBS -l nodes=1:ppn="+str(n_omp_threads)+" \n"
    if(mpi): header += "#PBS -l nodes="+str(n_nodes)+":ppn="+str(n_omp_threads)+" \n"
    header += "cd $PBS_O_WORKDIR\n"
elif(batch == 'slurm'):
    header =  "#!/bin/bash\n"
    header += "#SBATCH --account=XXX\n"
    header += "#SBATCH --nodes="+str(n_nodes)+"\n"
    header += "#SBATCH --ntasks-per-node=1\n"
    header += "#SBATCH --cpus-per-task="+str(n_omp_threads)+"\n"
    header += "#SBATCH --mem="+mem+"\n"
    header += "#SBATCH --time="+walltime+"\n\n"
ext_2bme = ".me2j.gz" # For two-body interaction (gzip)

NNF = "N3LO_EM500"
def set_input(params, hw=20, emax=14, e2max=28):
    params['rank'] = 2
    params['hw'] = hw
    params['renorm'] = 'bare'
    params['emax'] = emax
    params['e2max'] = e2max
    params['NNInt'] = NNF
    params['pmax2'] = 8
    params['NMesh2'] = 100
    params['Operators'] = "AxialV_Tz1-N2LO-NonLocal4-394"
    params['input_nn_file'] = f"{path_to_nninput}/{params['NNInt']}_kmax8_N100_Jmax8.bin"
    params['LECs'] = [-8.1e-1, -3.2, 5.4, 1.264, -0.120] # c1, c3, c4, cD, cE
    set_file_name_nn(params)

def get_script_name(params):
    fbase = "vHamil_MECalc_A" + str(params["rank"]) + \
            "_hw" + str(params["hw"])
    ms = ""
    fbase += ms
    fbase = "vHamil_MECalc_A2_" + params["file_name_nn"].split(ext_2bme)[0]
    return fbase

def gen_script(params, batch):
    exe = exe_file
    if(batch == 'slurm'): exe = "srun " + exe_file
    fbase = get_script_name(params)
    file_input = "Input_" + fbase + ".dat"
    file_log   = "log_" + fbase + ".dat"
    fsh = "run_" + fbase + ".sh"
    prt = header

    prt += 'echo "start ' +fsh + '..."\n'
    prt += "cat > "+file_input + " <<EOF\n"
    prt += "&input\n"
    for key, value in params.items():
        if(isinstance(value, str)):
            prt += "  " + str(key) + '= "' + str(value) + '" \n'
            continue
        if(isinstance(value, list)):
            prt += "  " + str(key) + "= "
            for x in value[:-1]:
                prt += str(x) + ", "
            prt += str(value[-1]) + "\n"
            continue
        prt += "  " + str(key) + '=' + str(value) + '\n'
    prt += "&end\n"
    prt += "EOF\n"
    if(batch == 'terminal'):
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
        f = f"c3_{params['LECs'][1]:.2f}_c4_{params['LECs'][2]:.2f}_cD_{params['LECs'][3]:.2f}_{params['renorm']}"
        if(params["renorm"] != 'Bare' and params["renorm"] != 'bare'): f += str(params["lambda"])
        f += '_hw' + str(params["hw"]) + '_emax' + str(params["emax"])
        f += '_e2max' + str(params["e2max"])
        f += ext_2bme
    except:
        f = "default"
    params["file_name_nn"] = f

def main():
    params = OrderedDict()
    for hw in [16]:
        set_input(params,hw,emax=12,e2max=24)
        fsh = gen_script(params, batch)
        if(batch == 'terminal' or batch == 'local'):
            cmd = "./" + fsh
        elif(batch == "pbs"):
            cmd = "qsub " + fsh
        elif(batch == 'slurm'):
            cmd = "sbatch " + fsh
        subprocess.call(cmd,shell=True)
        time.sleep(1)

if(__name__ == "__main__"):
    main()
