#!/usr/bin/env python3
import sys, os, itertools
HOME = os.path.expanduser("~")
sys.path.append(HOME)
import readline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mylib_python import myplt

def main():

    Lstring = ["S","P","D","F","G"]
    plt.rcParams = myplt.set_style(fontsize=0.15)
    fig, axs = myplt.set_canvas(r=2,c=5,width=14,height=6,shy=False)
    df = pd.read_csv("Nijmegen.csv", sep="\s*,\s*", engine='python')
    NNlist = ["N3LO_EM500"]
    for J in range(5):
        for S, P in itertools.product(range(2), [1,-1]):
            if(S==1): axs[S,J].set_xlabel(r'$T_{\rm lab}$ (MeV)')
            if(J==0): axs[S,J].set_ylabel(r'$\delta$')
            axs[S,J].set_title(r'$^{'+f'{2*S+1}'+"}$"+Lstring[J]+r"$_{"+str(J)+"}$")
            partial_wave = f'{2*S+1}{Lstring[J]}{J}'
            pn = 'pn'
            if(J==0 and P==-1 and S==1):
                axs[S,J].set_title(r'$^{3}$P$_{0}$')
                partial_wave = '3P0'
            Z = 0
            for idx, NN in enumerate(NNlist):
                filename = f"phase_shift_{NN}_bare_uncoupled.txt"
                data = np.loadtxt(filename, comments="#", skiprows=1)
                t, delta = [], []
                for i in range(len(data[:,0])):
                    if(data[i,0] != J): continue
                    if(data[i,1] != S): continue
                    if(data[i,2] != P): continue
                    if(data[i,3] != Z): continue
                    t.append(data[i,5])
                    delta.append(data[i,6]*180/np.pi)
                if(len(t)==0): continue
                axs[S,J].plot(t,delta,lw=0.8,label=NN)
            _ = df[(df['partial wave']==partial_wave) & (df['pn']==pn)]
            if(not _.empty): axs[S,J].errorbar(_['Tlab'], _['delta'], yerr=_['error'], ls='', marker='o', ms=2, mfc='none', c='k', capsize=3)

    axs[0,0].legend()
    fig.tight_layout()
    plt.savefig("uncoupled.pdf")

if(__name__ == "__main__"):
    main()
