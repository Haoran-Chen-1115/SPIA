# Stochastic path-integral approach (SPIA)
Stochastic path-integral approach for determing superconducting transition temperature Tc of general systems, including liquids, anharmonic solids and superionic systems.
The main theory can be found in the following references:

&emsp;**[1] Theory in liquids:** Huiying Liu, Ying Yuan, Donghao Liu, Xin-Zheng Li, and Junren Shi. Phys. Rev. Research 2, 013340 (2020). 

&emsp;**[2] Implementations in liquids:** Haoran Chen, Xiao-Wei Zhang, Xin-Zheng Li, and Junren Shi. Phys. Rev. B 104, 184516 (2021).

&emsp;**[3] Theory and implementations in anharmonic solids:** Haoran Chen, Junren Shi. Phys. Rev. B 106, 184501 (2022).

&emsp;**[4] Theory and implementations in superionic systems:** Haoran Chen, Junren Shi. arXiv:2305.04875.

## Interface to VASP
The interface to VASP 5.4.4 and VASP 6.4.1 is uploaded.

The additional subroutines output necessary intermediate variables for re-constructing Hamiltoninians and overlap matrices. Green's functions are calculated using the program in /main_new. The inerface to i-PI is also encoded. To apply the patch, run the following in the top level directory of VASP:
```
patch -p0 < vasp.6.4.1-SPIA.patch
```
Then you need to re-compile VASP.

There are several new flags to control the outputs:
* **LFHAM:** *Logical (default=.FALSE.)* Set to .TRUE. to enable output of intermediate variables for re-constructing Hamiltoninians.
* **NB:** *Integer (default=1)* The bead index in PIMD simulations, will be automatically set when combining with i-pi (see below).
* **LINDEX:** *Logical (default=.TRUE.)* Whether output some public informations in the first step (NB=1,NSTEP=1), including files: INDEX_\*, VKPT_\*, LM_INDEX_\*, SIMPI_\*, PWAV_\*, WAE_\*, WPS_\*, PSPNL_\*, PROJ_\*, ELECT_\*.
* **LKPOINTS:** *Logical (default=.FALSE.)* When set to .FALSE., the program only outputs information at the Gamma point.
* **NCONT:** *Integer (default=0)* The first step index will be NCONT+1, useful when performing a continuation job.
* **NJUMP:** *Integer (default=1)* Files are written every NJUMP steps.

If you are using machine-learning force field (MLFF), and the following flags can be used to post-processing the trajectory in **XDATCAR_old**:
* **LNODYN:** *Logical (default=.FALSE.)* Postions will be directly read from the file XDATCAR_old, and no molecular dynamics will be performed.
* **NCONT:** *Integer (default=0)* Start from NCONT-th configuration stored in XDATCAR_old

In order to perform PIMD simulations using i-pi (https://github.com/i-pi/i-pi), you need to prepare a makefile in the working directory containing:
```
# Makefile for the VASP example
#
# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.

.PHONY: all clean job 
all: job

# Command for execute VASP
VASP:=vasp_gpu
IPI:=i-pi

# Example of an job on a 8-gpu machine
define run_vasp
  for i in `seq 1 $1`; do \
    mkdir -p run_$$i; cp INCAR KPOINTS POTCAR POSCAR run_$$i; cd run_$$i; \
    echo "NB=$$i" >> INCAR; \
    j=$$(( ($$i-1)%8 )); \
    CUDA_VISIBLE_DEVICES=$$j $(VASP) & cd ..; \
  done;
endef


# 16 for the number of beads in PIMD simulations
job:
        $(IPI) input.xml & sleep 5; \
        $(call run_vasp,16) \
        wait

clean:
        rm -rf ./*simulation* ./*out* ./*.log* ./*run_* ./*RESTART* ./EXIT
```
and run the following command:
```
make job
```
Note that the line **echo** **"NB=$$i"** **>>** **INCAR;** \ is important when defining run_vasp!

## Preparations before calculating superconducting properties.
During the PIMD simulations, all outputs are placed in directories ROOT_DIR/run_$i/BEAD_$n, where $i is the bead index, and $n is the step index. You need to move the files to ROOT_DIR/BEAD_$n.
In addition, you need to prepare three more sets of files in the following directories:
* **ROOT_DIR/BEAD_1** One-step calculation using the initial ion configuration (a regular self-consistent calculation while setting LFHAM=.TRUE., LINDEX=.TRUE., NB=1). The equilibrium positions are used to anaylze the symmetry of the system. It provides informations like lattice parameters, Fourier transformation grids, and so on.
* **ROOT_DIR/BEAD_1_sym** One-step calculation using the equilibrium ion configuration (a regular self-consistent calculation while setting LFHAM=.TRUE., LINDEX=.TRUE., NB=1). The equilibrium positions are used to anaylze the symmetry of the system. In the case of a superionic system, you can use a POSCAR containing positions which possess the symmetry you assumed. For example, in Li2MgH16, you may use the solid-state position file. The file can contain different numbers of ions from your simulation.
* **ROOT_DIR/BEAD_1_primitive** One-step calculation using the equilibrium ion configuration in the primitive cell. It is used to calculate the equilibrium Bloch waves, and used to find the relation of k-points between the primitive Brillouin zone and supercell Brillouin zone.


## Main program
The main program is in the directory SPIA/main_new. A matlab script main.m is provided, describing the running procedure.
The program basically runs from 1_GEN_CORE to 6_Collect in steps:

* **0_Public:** Contains functions that are shared by multiple steps. 
Parameters of the calculations are the descriptions can be found in 0_Public/parameter.m.

* **1_GEN_CORE:** Calculate the interpolation table for Eqs. (15) and (A1) in Ref. [3]. 
```
>> GEN_CORE_s
```
* **2_Dens:** The file find_position_new.m is used to find the average positions in the base cell. The base cell is the original `primitive cell' to construct (non-)diagonal supercells. For the non-diagonal supercell technique, refer to  Jonathan H. Lloyd-Williams and Bartomeu Monserrat, Phys. Rev. B 92, 184301 (2015).
```
>> find_position_new
```

* **3_EFERMI:** Calculates EFERMI_AV=\<Ef\> as an initial guess of the Fermi energy.
```
>> CAC_EFERMI
```
      
* **4_Spec:** Calculate the eigen-energies distribution (spectral function) on the non-self-consistent k mesh. On a uniform k-mesh, the program is used to calculate the density of states and Fermi ernergy.
```
>> NCL=1; % Select the cluster number
>> Spec_main;
>> CAC_EFERMI_r;
```
  If L_path is set to true in 0_Public/parameter.m, the program can also be used to calculate the spectral function or band structure on the given k-point path.
        
* **5_SPIA:** The main program of SPIA. It calculates electron Green's functions of different ion configurations, and take the average. Parallelization is excuted by Gbar_ND.m. With the average Green's function, New_basis.m is used to calculates the EPC-renormalized Bloch bases using the average Green's function according to Eq. (2) of Ref. [3]. Finally, Tbar_ND.m calculates the electron-electron pair scateering amplitude according to Eqs. (1) and (5) of Ref. [3].
```
>> NCL=1; % Select the cluster number
>> SPIA;
```      

* **6_Collect:** Solves the Bethe-Salpeter equation to calculate the effective electron-electron interactions, and solve linearized Eliashberg equations to solve Tc. The Lorentz smearing for calculating EPC parameters range from 0.02 *eV* to 1.6 *eV* with a step 0.02 *eV* by default.
```
>> Collectq_NDsc_lor_sym;
```   
