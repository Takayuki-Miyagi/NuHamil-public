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
n_nodes=1               # Basically, you can use mpi, but if there is no NNrel file, it will stop. I recommend you to try without mpi for a first run.
n_omp_threads=48        # threads per node
mem = "125G"            # memory per node
walltime = "0-01:00"    # request time
mpi = False
if(n_nodes > 1): mpi=True
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

# 3-body SRG space
ramp_space = "ramp40-5-36-7-32-9-28-11-24" # This should be enough up to 40Ca or so. From Phys. Rev. C90, 024325 (2014).
#ramp_space = "ramp40-7-38-9-36" # This should be enough up to A~100. From Phys. Lett. B736, 119 (2014).

# file format
ext_3bme = ".me3j.gz"


"""
NN
"N3LO_EM500", "LO_EMN500", "NLO_EMN500", "N2LO_EMN500", "N3LO_EMN500", "N4LO_EMN500"
"""
NNF = "N3LO_EM500"

"""
3N
"EM1.8-2.0",
"EM2.0-2.0","EM2.0-2.5","EM2.2-2.0","EM2.8-2.0","PWA2.0-2.0",
"N2LOsat","N3LO-NoggaA","N3LO-NoggaB","3NF400",
"3NF500-BetaDecay-Local","N3LOEMlnl",
"N4LOEMNlnl","DN2LOGO450","DNLOGO450","DN2LOGO394"
"""
TNF = "N3LOEMlnl"

def set_input(params, hw=20, hw_target=None, emax=14, e2max=28, e3max=16):
    # basic parameters
    params['rank'] = 3
    params['hw'] = hw
    if(hw_target != None): params['hw_target']=hw_target # hw conversion
    params['renorm'] = 'srg'
    params['lambda'] = 2.0
    params['emax'] = emax
    params['e2max'] = e2max
    params['e3max'] = e3max
    params['NNInt'] = NNF
    params['input_nn_file'] = f"{path_to_nninput}/{params['NNInt']}_kmax8_N100_Jmax8.bin"
    params['jmax3'] = params['e3max'] * 2 + 3
    params['ramp'] = ramp_space
    params['genuine_3bf'] = True
    params["3nf"] = TNF
    #params["path_to_Hebeler_files"]= "" # You need set this properly if you run with the non-locally regulated 3N interaction, such as "EM1.8-2.0" and "DN2LOGO394"
    if(params['genuine_3bf']): set_three_body_force(params)
    if(params['rank']==3): set_file_name_3n(params)

def get_script_name(params):
    fbase = "vHamil_MECalc_A" + str(params["rank"]) + \
            "_hw" + str(params["hw"])
    ms = ""
    fbase += ms
    if(params['rank']>2):
        if(params['genuine_3bf']):
            try:
                fbase += "_" + params["3nf"]
            except:
                pass
    if("3nf" in params): del params["3nf"]
    fbase = "vHamil_MECalc_A3_" + params["file_name_3n"].split(ext_3bme)[0]
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

def set_three_body_force(params):
    warning = False
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

def set_file_name_3n(params):
    try:
        f = ""
        if("only_hf_monopole" in params and params["only_hf_monopole"]): f = "Monopole_"
        if("only_no2b_elements" in params and params["only_no2b_elements"]): f = "NO2B_"
        if("lab_3bme_precision" in params and params["lab_3bme_precision"]!="single"): f += "{:s}_".format(params["lab_3bme_precision"])
        f += "ThBME"
        if(params["3nf"]=="EM1.8-2.0"):
            f += "_EM1.8_2.0"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="EM2.0-2.0"):
            f += "_EM2.0_2.0"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="EM2.0-2.5"):
            f += "_EM2.0_2.5"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="EM2.2-2.0"):
            f += "_EM2.2_2.0"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="EM2.8-2.0"):
            f += "_EM2.8_2.0"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="PWA2.0-2.0"):
            f += "_PWA2.0_2.0"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="N2LOsat"):
            f += "_N2LOsat"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="DN2LOGO394"):
            f += "_DN2LOGO394"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="DN2LOGO450"):
            f += "_DN2LOGO450"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        elif(params["3nf"]=="DNLOGO450"):
            f += "_DNLOGO450"
            if(params["Regulator"] != "NonLocal"):
                print("Error: Regulator="+params["Regulator"])
            if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
            if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
        else:
            if(params["renorm"] == "Bare" or params["renorm"] == "bare"):
                if(params["genuine_3bf"]):
                    if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
                    if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
                    f += '_c1_' + str(params["LECs"][0]) + '_c3_' + str(params["LECs"][1]) + '_c4_' + str(params["LECs"][2]) +\
                            '_cD_' + str(params["LECs"][3]) + '_cE_' + str(params["LECs"][4])
                    f += "_" + params["Regulator"] + str(params["RegulatorPower"])
                    if(params["Regulator"] == "LNL"):
                        f += "_" + str(params["lambda_3nf_local"]) + "_" + str(params["lambda_3nf_nonlocal"])
                    if(params["Regulator"] == "Local"):
                        f += "_" + str(params["lambda_3nf_local"])
                    if(params["Regulator"] == "NonLocal"):
                        f += "_" + str(params["lambda_3nf_nonlocal"])
            else:
                f += '_' + str(params["renorm"]) + str(params["lambda"]) + "_" + params["ramp"]
                f += '_' + str(params["NNInt"])
                if(params["genuine_3bf"]):
                    if("j3max_initial_3nf" in params): f += "_3NFJmax" + str(params["j3max_initial_3nf"])
                    if("jmax3" in params and params['jmax3']<2*params['e3max']+3): f += "_JJmax" + str(params["jmax3"])
                    f += '_c1_' + str(params["LECs"][0]) + '_c3_' + str(params["LECs"][1]) + '_c4_' + str(params["LECs"][2]) +\
                            '_cD_' + str(params["LECs"][3]) + '_cE_' + str(params["LECs"][4])
                    f += "_" + params["Regulator"] + str(params["RegulatorPower"])
                    if(params["Regulator"] == "LNL"):
                        f += "_" + str(params["lambda_3nf_local"]) + "_" + str(params["lambda_3nf_nonlocal"])
                    if(params["Regulator"] == "Local"):
                        f += "_" + str(params["lambda_3nf_local"])
                    if(params["Regulator"] == "NonLocal"):
                        f += "_" + str(params["lambda_3nf_nonlocal"])
        if("hw_target" in params): f += '_IS_hw' + str(params["hw_target"]) + "from"+ str(params["hw"])
        if(not "hw_target" in params): f += '_IS_hw' + str(params["hw"])
        f += '_ms'+ str(params["emax"]) + '_' + str(params["e2max"]) +"_"+ str(params["e3max"])
        f += ext_3bme
    except:
        f = "default"
    params["file_name_3n"] = f


def main():
    params = OrderedDict()
    hw = 30
    for hw_target in [16]:
        set_input(params,hw,hw_target=hw_target,emax=14,e2max=28,e3max=16)
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
