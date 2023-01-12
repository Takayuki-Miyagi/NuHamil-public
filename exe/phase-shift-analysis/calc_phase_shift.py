#!/usr/bin/env python
#-*- coding:utf-8 -*-
# script to submit the job
# execute ./NuHamil.py
import os
import sys
import subprocess
from collections import OrderedDict
import time

HOME = os.path.expanduser("~")
path_to_nninput = f"{HOME}/NuHamil/input_nn_files" # this is default; please modify here if you did not download the code in your home directory.
exe_file = 'NuHamil.exe'

def set_input(params, NNint):
    # basic parameters
    params["is_PhaseShift"] = True
    params['rank'] = 2
    params['renorm'] = 'bare'
    #params['renorm'] = 'srg'
    params['lambda'] = 2.0
    params["renorm_space2"] = "mom"

    params['NNInt'] = NNint
    params['input_nn_file'] = f"{path_to_nninput}/{params['NNInt']}_kmax8_N100_Jmax8.bin"
    params['J2max_NNint'] = 4
    params['Tlab_phase_shift'] = '1,5,10,25,50,75,100,125,150,175,200,225,250,275,300,325,350'
    params['file_name_phase_shift'] = 'phase_shift_' + params['NNInt'] + '_' + params['renorm']
    if(params["renorm"] != 'Bare' and params["renorm"] != 'bare'): params['file_name_phase_shift'] += str(params["lambda"])
def get_script_name(params): return "vHamil_" + str(params["file_name_phase_shift"])

def gen_script(params):
    exe = exe_file
    fbase = get_script_name(params)
    file_input = "Input_" + fbase + ".dat"
    file_log   = "log_" + fbase + ".dat"
    fsh = "run_" + fbase + ".sh"
    prt = ""
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
    prt += exe + " " + file_input + "\n"
    prt += "rm " + file_input + "\n"
    f = open(fsh, "w")
    f.write(prt)
    f.close()
    os.chmod(fsh, 0o755)
    return fsh

def main(machinename=None):
    params = OrderedDict()
    NNlist = ["N3LO_EM500"]
    for NN in NNlist:
        set_input(params,NN)
        fsh = gen_script(params)
        cmd = "./" + fsh
        subprocess.call(cmd,shell=True)
        time.sleep(1)

if(__name__ == "__main__"):
    main()
