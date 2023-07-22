#!/usr/bin/env python3
import numpy as np
from scipy.special import spherical_jn
from scipy.special import gamma
from matplotlib import pyplot as plt

"""
Local projection based on Eq. (9) and (10) in K. A. Wendt, R. J. Furnstahl, and S. Ramanan, Phys. Rev. C 86, 014003 (2012).
<Lp,S,J | V | L,S,J >
"""

filename = 'N3LO_EM500_kmax8_N100_Jmax8.bin'
Lp, L, S, J, Tz = 0, 0, 1, 1, 0 # np 3S1
#Lp, L, S, J, Tz = 0, 2, 1, 1, 0 # np 3S1-3D1
#Lp, L, S, J, Tz = 2, 2, 1, 1, 0 # np 3D1
#Lp, L, S, J, Tz = 0, 0, 0, 0, 0 # np 1S0
r = np.linspace(0,6,100)

fig, ax = plt.subplots(ncols=1, nrows=1)
vmin, vmax = -60, 120 # min and max

def read_nn_file(filename, J=0, Prty=1, S=0, Tz=0, byte_order='little', channel='--'):
    """
    filename: input momentum-space file name
    J: total angular momentum
    S: total spin
    Prty: Parity, 1 or -1
    Tz: -1 (pp), 0 (pn), or 1 (nn)
    channel: only valid for spin triplet coupled channel. '--' (<L=J-1| V | L=J-1>) '-+' (<L=J-1|V|L=J+1>), or '++' (<L=J+1|V|L=J+1>)
    """
    with open(filename, "rb") as fp:
        NMesh = int.from_bytes(fp.read(4), byteorder=byte_order, signed=True)
        Jmax = int.from_bytes(fp.read(4), byteorder=byte_order, signed=True)
        NChan = int.from_bytes(fp.read(4), byteorder=byte_order, signed=True)
        Mesh = np.zeros(NMesh)
        Mesh = np.frombuffer(fp.read(8*NMesh),count=NMesh,dtype=np.float64)
        Weights = np.frombuffer(fp.read(8*NMesh),count=NMesh,dtype=np.float64)
        for ichan in range(NChan):
            J_read = int.from_bytes(fp.read(4), byteorder=byte_order, signed=True)
            Prty_read = int.from_bytes(fp.read(4), byteorder=byte_order, signed=True)
            S_read = int.from_bytes(fp.read(4), byteorder=byte_order, signed=True)
            Tz_read = int.from_bytes(fp.read(4), byteorder=byte_order, signed=True)
            Ndim = int.from_bytes(fp.read(4), byteorder=byte_order, signed=True)
            V = np.zeros(Ndim*Ndim)
            V = np.frombuffer(fp.read(8*Ndim*Ndim), count=Ndim*Ndim, dtype=np.float64)

            if(J_read==J and Prty_read==Prty and S_read==S and Tz_read==Tz): break
        fp.close()
    if(NMesh==Ndim):
        z = V
    else:
        if(channel=='--'):
            z = np.reshape(np.reshape(V, (Ndim,Ndim))[:NMesh,:NMesh], NMesh**2)
        elif(channel=='-+'):
            z = np.reshape(np.reshape(V, (Ndim,Ndim))[:NMesh,NMesh:], NMesh**2)
        elif(channel=='+-'):
            z = np.reshape(np.reshape(V, (Ndim,Ndim))[NMesh:,:NMesh], NMesh**2)
        elif(channel=='++'):
            z = np.reshape(np.reshape(V, (Ndim,Ndim))[NMesh:,NMesh:], NMesh**2)
    return Mesh, Weights, np.reshape(z, (NMesh,NMesh))

def local_projection_0l(r, l, kmesh, weights, mat):
    v = np.zeros(r.size)
    for idx, k in enumerate(kmesh):
        v += k**2 * weights[idx] * spherical_jn(l, k*r) * mat[idx, 0] # note: index 0 is approximation.
    return v

def local_projection_ll(r, lp, l, kmesh, weights, mat):
    norm = 4 / np.sqrt(np.pi) * gamma((l+3)/2) / gamma(l/2)
    if(lp!=l): norm *= -1
    v = np.zeros(r.size)
    for i, kp in enumerate(kmesh):
        for j, k in enumerate(kmesh):
            v += weights[i]*weights[j]*(kp**2/k) * spherical_jn(lp, kp*r) * mat[i,j]
    return v * norm

channel = '--'
if(S==1 and Lp==J-1 and L==J-1): channel='--'
if(S==1 and Lp==J-1 and L==J+1): channel='-+'
if(S==1 and Lp==J+1 and L==J-1): channel='+-'
if(S==1 and Lp==J+1 and L==J+1): channel='++'

k, w, v = read_nn_file(filename, J=J, Prty=(-1)**L, S=S, Tz=Tz, channel=channel)
if(Lp == 0):
    vlocal = local_projection_0l(r, L, k, w, np.transpose(v))
elif(L  == 0):
    vlocal = local_projection_0l(r,Lp, k, w, v)
else:
    vlocal = local_projection_ll(r,Lp, L, k, w, v)
ax.plot(r, vlocal, ls='-')
ax.set_xlim(0,4)
ax.set_ylim(vmin,vmax)
ax.set_xlabel(r'$r$ (fm)')
ax.set_ylabel(r'$V(r)$ (MeV)')
fig.tight_layout()
plt.show()
