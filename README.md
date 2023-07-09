# NuHamil

[![DOI](https://zenodo.org/badge/534840092.svg)](https://zenodo.org/badge/latestdoi/534840092)

<img align="center" width="120" height="100" src="LogoNuHamil.png">

This code is to generate the nucleon-nucleon (NN) and three-nucleon (3N) matrix elements.
The generated NN and 3N files can be used as inputs of the many-body solvers such as imsrg++ code (https://github.com/ragnarstroberg/imsrg) and HartreeFock code (https://github.com/Takayuki-Miyagi/HartreeFock).
It is also possible to compute the few body system with the Jacobi coordinate no-core shell model.

# Requirements
* Fortran compiler
* BLAS, LAPACK libraries
* GNU Scientific library
* HDF5 library
* zlib

# Installation
Download the code:
```
cd ~
git clone https://github.com/Takayuki-Miyagi/NuHamil-public.git
```
(Going to your home directory is not mandatory, but recommended.)

Change the directory and download the submodules:
```
cd NuHamil
git submodule init
git submodule update
```
The executable file can be built via make command.
By default, it will be compiled with gfortran without MPI.
If you want to use the intel compiler and/or turn on the MPI parallelization, please edit the Makefile by yourself.
```
make -j
make install
echo 'export PATH=$PATH:$HOME/bin' >> $HOME/.bashrc
source $HOME/.bashrc
```
By default, the symbolic link is created in the HOME/bin directory.
If you want to change it, please change the "INSTLDIR" in the Makefile.

# Some examples
Running a job is managed by a simple python script.
Some example scripts are prepared.
Note that you need to change "path_to_nninput" in the python scripts if you did not download the code in your home directory.
* deuteron properties

This should be finished with in a few seconds.
```
cd $HOME/exe/few-body
python3 deuteron.py
```

* triton properties

As the $N_{\rm max}$ truncation is sufficiently large $N_{\rm max}=40$, this will take a while. (~ half hour with 16 OpenMP threads on MacBook Pro, 2.3 GHz 8-Core Intel Core i9)
```
cd $HOME/exe/few-body
python3 triton_helium3.py
```

* lab-frame HO NN matrix elements

As the $e_{\rm max}$ truncation is reasonably large $e_{\rm max}=14$, this will take a while. (~ 15 min. with 16 OpenMP threads on MacBook Pro, 2.3 GHz 8-Core Intel Core i9)
```
cd $HOME/exe
python3 NuHamil_2BME.py
```


* lab-frame HO 3N matrix elements

This will not work on the local machine because of the computational power.
I recommend submitting the job to a workstation or supercomputer.
You need to edit NuHamil_3BME.py script depending on which job submission system is installed on the machine.
This will take ~ 2 days, using 32 OpenMP threads (without MPI). The job is also memory expensive, and you need ~ 200GB RAM.
```
cd $HOME/exe
python3 NuHamil_3BME.py
```

# Quick benchmarks
Using the NN+3N interaction generated in the above example and the imsrg++ code (https://github.com/ragnarstroberg/imsrg), the ground-state energies of doubly-magic nuclei are the following.
| Nucleus | Interaction| $\lambda_{\rm SRG} \ ({\rm fm}^{-1})$|$e_{\rm max}$ | $E_{\rm 3max}$| $\hbar\omega \ ({\rm MeV})$|$E_{\rm g.s.} \ ({\rm MeV})$ |
|---------------:|---------------------------:|----:|---:|----:|----:|----------:|
|$^{4}{\rm He}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 4  | 16  | 16  | -27.221   |
|$^{4}{\rm He}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 6  | 16  | 16  | -28.221   |
|$^{4}{\rm He}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 8  | 16  | 16  | -28.569   |
|$^{4}{\rm He}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 10 | 16  | 16  | -28.669   |
|$^{4}{\rm He}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 12 | 16  | 16  | -28.690   |
|$^{16}{\rm O}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 4  | 16  | 16  | -115.540  |
|$^{16}{\rm O}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 6  | 16  | 16  | -123.445  |
|$^{16}{\rm O}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 8  | 16  | 16  | -126.270  |
|$^{16}{\rm O}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 10 | 16  | 16  | -127.016  |
|$^{16}{\rm O}$  | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 12 | 16  | 16  | -127.161  |
|$^{40}{\rm Ca}$ | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 4  | 16  | 16  | -284.589  |
|$^{40}{\rm Ca}$ | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 6  | 16  | 16  | -322.582  |
|$^{40}{\rm Ca}$ | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 8  | 16  | 16  | -336.124  |
|$^{40}{\rm Ca}$ | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 10 | 16  | 16  | -340.546  |
|$^{40}{\rm Ca}$ | ${\rm NN-N}^{3}{\rm LO+3N}_{\rm lnl}$| 2.0 | 12 | 16  | 16  | -341.678  |


# Remarks
Unfortunately, non-locally regulated 3N interaction cannot be fully suppported.
Feel free to contact me (miyagi@theorie.ikp.physik.tu-darmstadt.de) if you are interested in such interactions.
I can provide the HO 3N matrix element files in the Jacobi coordinate or in the laboratory frame.

# How to cite
If you use this code in your research, please cite [T. Miyagi, Eur. Phys. J. A 59, 150 (2023)](https://link.springer.com/article/10.1140/epja/s10050-023-01039-y).


# Acknowledgement
**NuHamil** code is greatly inspired by the **manyeff** code by P.Navratil and **VRenormalize** in Computational Environment for Nuclear Structure (CENS) project (https://github.com/ManyBodyPhysics/CENS).
The code uses VODE library by G.D.Byrne and S.Thompson.
I thank P.Navratil, N.Shimizu, and N.Tsunoda for the discussions, optimizations, and parallelizations.
I also thank P.Arthuis, A.Belly, M.Heinz, B.S.Hu, S.R.Stroberg, and A.Tichai for testing the code and useful feedback.

