#!/usr/bin/env python
#-*- coding:utf-8 -*-
# script to submit the job
# execute ./NuHamil.py
import os
import sys
import subprocess
from collections import OrderedDict
exe = 'NuHamil.exe'
def set_input(params):
    # basic parameters
    params["file_convert"] = True
    # NN file
    params['rank'] = 2
    params['emax'] = 0
    params['e2max'] = 0
    params['emax_convert'] = 0
    params['e2max_convert'] = 0
    params["file_name_nn_original"]=""
    params["file_name_nn_converted"]=""

    # NNN file
    #params['rank'] = 3
    #params['emax'] = 0
    #params['e2max'] = 0
    #params['e3max'] = 0
    #params["file_name_3n_original"]=""
    #params["file_name_3n_converted"]=""

def gen_script(params):
    try:
        if( "file_name_nn_original" in params ): fbase = "Converter_A" + str(params["rank"]) + "_" + os.path.splitext(os.path.basename(params["file_name_nn_original"]))[0]
        if( "file_name_3n_original" in params ): fbase = "Converter_A" + str(params["rank"]) + "_" + os.path.splitext(os.path.basename(params["file_name_3n_original"]))[0]
    except:
        fbase = "Converter"
    file_input = "Input_" + fbase + ".dat"
    file_log   = "log_" + fbase + ".dat"
    fsh = "run_" + fbase + ".sh"
    prt = ""
    prt += 'echo "run ' +fsh + '..."\n'
    prt += "cat > "+file_input + " <<EOF\n"
    prt += "&input\n"
    for key, value in params.items():
        if(isinstance(value, str)):
            prt += str(key) + '= "' + str(value) + '" \n'
            continue
        if(isinstance(value, list)):
            prt += str(key) + "= "
            for x in value[:-1]:
                prt += str(x) + ", "
            prt += str(value[-1]) + "\n"
            continue
        prt += str(key) + '=' + str(value) + '\n'
    prt += "&end\n"
    prt += "EOF\n"
    prt += exe + " " + file_input + "\n"
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

if(__name__ == "__main__"):
    main()

