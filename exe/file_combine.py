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

n_omp_threads=48        # threads per node
mem = "125G"            # memory per node
walltime = "0-4:00"    # request time
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
    header += "#PBS -l nodes=1:ppn="+str(n_omp_threads)+" \n"
    header += "cd $PBS_O_WORKDIR\n"
elif(batch == 'slurm'):
    header =  "#!/bin/bash\n"
    header += "#SBATCH --account=XXX\n"
    header += "#SBATCH --nodes=1\n"
    header += "#SBATCH --ntasks-per-node=1\n"
    header += "#SBATCH --cpus-per-task="+str(n_omp_threads)+"\n"
    header += "#SBATCH --mem="+mem+"\n"
    header += "#SBATCH --time="+walltime+"\n\n"

def set_input(params):
    params['rank']=2
    params['emax']=4
    params['e2max']=8
    params['e3max']=4
    params['files_combined'] = "file1.me2j.gz,file2.me2j.gz,file3.me2j.gz"
    params['weights_combined'] = "1.0,1.0,1.0"
    params['file_name_nn'] = "out.me2j.gz"

    # NNN
    #params['rank']=3
    #params['emax']=4
    #params['e2max']=8
    #params['e3max']=4
    #params['files_combined'] = ""
    #params['weights_combined'] = ""
    #params['file_name_3n'] = ""

def get_script_name(params):
    if(params['rank']==2): fbase = "vHamil_MECalc_A2_" + params["file_name_nn"].split(ext_2bme)[0]
    if(params['rank']==3): fbase = "vHamil_MECalc_A3_" + params["file_name_3n"].split(ext_3bme)[0]
    return fbase

def gen_script(params):
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

def main()
    params = OrderedDict()
    set_input(params)
    fsh = gen_script(params)
    cmd = "./" + fsh
    subprocess.call(cmd,shell=True)
    time.sleep(1)

if(__name__ == "__main__"):
    main()
