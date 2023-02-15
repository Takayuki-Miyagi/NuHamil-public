#!/usr/bin/env python
#-*- coding:utf-8 -*-
# script to submit the job
# execute ./NuHamil.py
import os
import sys
import subprocess
from collections import OrderedDict
import time, itertools

HOME=os.path.expanduser("~")
path_to_nninput = f"{HOME}/NuHamil-public/input_nn_files" # this is default; please modify here if you did not download the code in your home directory.
exe_file = 'NuHamil.exe'
operators = "Rp2,Rn2,Rm2,Mmoment_IA"

use_genuine_3NF = True
path_to_Hebeler_files=None
triton = [1,1,1, 1]
helium3= [1,1,1,-1]
NNF = 'N3LO_EM500'

def set_input(params, hw=20, N3max=40):
    # basic parameters
    params['rank'] = 3
    params['hw'] = hw
    params['renorm'] = 'bare'
    params['NNInt'] = NNF
    params['input_nn_file'] = f"{path_to_nninput}/{params['NNInt']}_kmax8_N100_Jmax8.bin"
    params['Operators'] = operators

    # two- or three-body matrix element #
    params["bra"] = triton
    params["ket"] = triton
    params['trans2lab'] = False
    #params['NN_only'] = True # (turnning off 3NF)
    params['N3max'] = N3max
    params['genuine_3bf'] = True
    params["3nf"] = 'N3LOEMlnl'
    if(params['genuine_3bf']): set_three_body_force(params)

def get_script_name(params):
    fbase = "vHamil_A3_hw" + str(params["hw"])
    ms = ""
    ms += "_Nmax" + str(params["N3max"])
    pb = "+"
    pk = "+"
    if(params["bra"][1] == -1): pb = "-"
    if(params["ket"][1] == -1): pk = "-"
    ms += "_j"+ str(params["bra"][0]) + "p" + pb + "t" + str(params["bra"][2]) + "tz" + str(params["bra"][3])
    ms += "_j"+ str(params["ket"][0]) + "p" + pk + "t" + str(params["ket"][2]) + "tz" + str(params["ket"][3])
    fbase += ms
    if(params['genuine_3bf']):
        try:
            fbase += "_" + params["3nf"]
        except:
            pass
    if('3nf' in params): del params["3nf"]
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

def set_three_body_force(params):
    if(params["rank"] < 3): return
    if(not params['genuine_3bf']): return
    if(params["3nf"]=="EM1.8-2.0"): # use with SRG evolved two-body interaction (EM500)
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 4
        params['LECs'] = [-8.1e-1, -3.2, 5.4, 1.264, -0.120]
        params['lambda_3nf_nonlocal'] = 2 * 197.3269788
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="EM2.0-2.0"): # use with SRG evolved two-body interaction (EM500)
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 4
        params['LECs'] = [-8.1e-1, -3.2, 5.4, 1.271, -0.131]
        params['lambda_3nf_nonlocal'] = 2 * 197.3269788
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="EM2.0-2.5"): # use with SRG evolved two-body interaction (EM500)
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 4
        params['LECs'] = [-8.1e-1, -3.2, 5.4, -0.292, -0.592]
        params['lambda_3nf_nonlocal'] = 2.5 * 197.3269788
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="EM2.2-2.0"): # use with SRG evolved two-body interaction (EM500)
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 4
        params['LECs'] = [-8.1e-1, -3.2, 5.4, 1.214, -0.137]
        params['lambda_3nf_nonlocal'] = 2 * 197.3269788
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="EM2.8-2.0"): # use with SRG evolved two-body interaction (EM500)
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 4
        params['LECs'] = [-8.1e-1, -3.2, 5.4, 1.278, -0.078]
        params['lambda_3nf_nonlocal'] = 2 * 197.3269788
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="PWA2.0-2.0"): # use with SRG evolved two-body interaction (EM500)
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 4
        params['LECs'] = [-7.6e-1, -4.78, 3.96, -3.007, -0.686]
        params['lambda_3nf_nonlocal'] = 2 * 197.3269788
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="N2LOsat"): # use with bare two-body interaction (N2LO_sat)
        params['NNInt'] = 'N2LO_sat'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 3
        params['LECs'] = [-1.12152120, -3.92500586, 3.76568716, 0.81680589, -0.03957471]
        params['lambda_3nf_nonlocal'] = 450
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="DN2LOGO394"):
        params['NNInt'] = 'DN2LOGO394'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 4
        params['LECs'] = [-0.74, -3.622246, 2.446123, 0.081, -0.002]
        params['lambda_3nf_nonlocal'] = 394
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="DN2LOGO450"):
        params['NNInt'] = 'DN2LOGO450'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 3
        params['LECs'] = [-0.74, -3.622246, 2.446123, -0.454, -0.186]
        params['lambda_3nf_nonlocal'] = 450
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="DNLOGO450"):
        params['NNInt'] = 'DNLOGO450'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 3
        params['LECs'] = [0.0, -2.972246, 1.486123, 0.0, 0.0]
        params['lambda_3nf_nonlocal'] = 450
        if(params['rank'] > 2): params['renorm'] = 'bare'
        return

    if(params["3nf"]=="N3LO-NoggaA"):
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 2
        params['LECs'] = [-8.1e-1, -3.2, 5.4, -1.11, -0.66]
        params['lambda_3nf_nonlocal'] = 500
        return

    if(params["3nf"]=="N3LO-NoggaB"):
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'NonLocal'
        params['RegulatorPower'] = 2
        params['LECs'] = [-8.1e-1, -3.2, 5.4, 8.14, -2.02]
        params['lambda_3nf_nonlocal'] = 500
        return

    if(params["3nf"]=="3NF400"): # use with EM500 two-body interaction
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'Local'
        params['RegulatorPower'] = 2
        params['LECs'] = [-8.1e-1, -3.2, 5.4, -2.e-1, 9.8e-2]
        params['lambda_3nf_local'] = 400
        return

    if(params["3nf"]=="3NF500-BetaDecay-Local"): # use with EM500 two-body interaction
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'Local'
        params['RegulatorPower'] = 2
        params['LECs'] = [-8.1e-1, -3.2, 5.4, 0.83, -0.052]
        params['lambda_3nf_local'] = 500
        return

    if(params["3nf"]=="N3LOEMlnl"): # use with EM500 two-body interaction
        params['NNInt'] = 'N3LO_EM500'
        params['Regulator'] = 'LNL'
        params['RegulatorPower'] = 2
        params['LECs'] = [-8.1e-1, -3.2, 5.4, 7.e-1, -6.e-2]
        params['lambda_3nf_local'] = 650
        params['lambda_3nf_nonlocal'] = 500
        return

    if(params["3nf"]=="N4LOEMNlnl"): # use with EMN500 two-body interaction
        params['NNInt'] = 'N4LO_EMN500'
        params['Regulator'] = 'LNL'
        params['RegulatorPower'] = 2
        params['LECs'] = [-7.3e-1, -3.38, 1.69, -1.8, -3.1e-1]
        params['lambda_3nf_local'] = 650
        params['lambda_3nf_nonlocal'] = 500
        return

    print("Unknown three-body force keyword" + params["3nf"])
    sys.exit()

def main():
    params = OrderedDict()
    for hw, N3max in itertools.product([30],range(40,42,2)):
        set_input(params,hw,N3max=N3max)
        fsh = gen_script(params)
        cmd = "./" + fsh
        subprocess.call(cmd,shell=True)
        time.sleep(1)

if(__name__ == "__main__"):
    main()
