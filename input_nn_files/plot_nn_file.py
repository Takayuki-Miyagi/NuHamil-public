#!/usr/bin/env python3
import os, sys
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(ncols=1, nrows=1)
vmin, vmax = -30, 30 # min and max value for colormap

def main(filename, J=0, Prty=1, S=0, Tz=0, byte_order='little', channel='--'):
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
        fp.seek(NMesh*8, 1) # skip weights
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
        elif(channel=='++'):
            z = np.reshape(np.reshape(V, (Ndim,Ndim))[NMesh:,NMesh:], NMesh**2)

    Meshes = np.reshape(np.tile(Mesh,NMesh), (NMesh,NMesh))
    x = Meshes.reshape(NMesh*NMesh)
    y = np.transpose(Meshes).reshape(NMesh*NMesh)
    X=np.arange(np.min(x), np.max(x), 0.05)
    Y=np.arange(np.min(y), np.max(y), 0.05)
    Z=griddata((x,y), z, (X[None,:], Y[:,None]), method='cubic')
    im = ax.imshow(Z, norm=colors.Normalize(vmin=vmin, vmax=vmax), \
            extent=(np.min(x), np.max(x), np.max(y), np.min(y)), cmap="jet")
    divider=make_axes_locatable(ax)
    cax = divider.append_axes('right','5%',pad='3%')
    cbar = fig.colorbar(im,ax=ax, cax=cax)
    cbar.set_label(r"$V(k',k)$ (MeV fm$^{-3}$)")
    ax.set_xlabel("$k$ (fm$^{-1}$)")
    ax.set_ylabel("$k'$ (fm$^{-1}$)")
    plt.show()

if(__name__=='__main__'):
    main(sys.argv[1])

