#!/usr/bin/env python
#-*- coding:utf-8 -*-
# script to submit the job
# execute ./NuHamil.py 
import os
import sys
import subprocess
from collections import OrderedDict
import time

HOME=os.path.expanduser("~")
path_to_nninput = f"{HOME}/NuHamil-public/input_nn_files" # this is default; please modify here if you did not download the code in your home directory.
exe_file = 'NuHamil.exe'
operators = "hamil,Mmoment_IA,Qmoment_IA,Radius"
NNF = 'N3LO_EM500'

def set_input(params, hw=20):
    # basic parameters
    params['rank'] = 2
    params['hw'] = hw
    params['renorm'] = 'bare'
    params['lambda'] = 1.8
    params['Operators'] = operators

    params["bra"] = [1,1,0,0]
    params["ket"] = [1,1,0,0]
    params['trans2lab'] = False
    params['NNInt'] = "N3LO_EM500"
    params['input_nn_file'] = f"{path_to_nninput}/{params['NNInt']}_kmax8_N100_Jmax8.bin"

def get_script_name(params):
    fbase = "vHamil_H2" + "_hw" + str(params["hw"])
    ms = ""
    pb = "+"
    pk = "+"
    if(params["bra"][1] == -1): pb = "-"
    if(params["ket"][1] == -1): pk = "-"
    ms += "_j"+ str(params["bra"][0]) + "p" + pb + "t" + str(params["bra"][2]) + "tz" + str(params["bra"][3])
    ms += "_j"+ str(params["ket"][0]) + "p" + pk + "t" + str(params["ket"][2]) + "tz" + str(params["ket"][3])
    fbase += ms
    return fbase

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

def main():
    params = OrderedDict()
    for hw in [30]:
        set_input(params,hw)
        fsh = gen_script(params)
        cmd = "./" + fsh
        subprocess.call(cmd,shell=True)
        time.sleep(1)

if(__name__ == "__main__"):
    main()
