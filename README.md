# Stochastic path-integral approach (SPIA)
Stochastic path-integral approach for determing superconducting transition temperature Tc of anharmonic systems.
The main theory can be found in the following references:

&emsp;**[1] Theory in liquids:** Huiying Liu, Ying Yuan, Donghao Liu, Xin-Zheng Li, and Junren Shi. Phys. Rev. Research 2, 013340 (2020). 

&emsp;**[2] Implementations in liquids:** Haoran Chen, Xiao-Wei Zhang, Xin-Zheng Li, and Junren Shi. Phys. Rev. B 104, 184516 (2021).

&emsp;**[3] Theory and implementations in anharmonic solids:** Haoran Chen, Junren Shi. arXiv:2205.07247.

Instructions are not complete yey.

## Main program
The main program is in the directory SPIA/main_nscf. 
The program basically runs from 1_GEN_CORE to 7_Collect in steps:

* **0_Public:** Contains functions that are shared by multiple steps. 
Parameters of the calculations are the descriptions can be found in 0_Public/parameter.m.

* **1_GEN_CORE:** Calculate the interpolation table for Eqs. (15) and (A1) in Ref. [3]. 
```
>> GEN_CORE_s
```

* **2_GEN_BASIS:** Calculates the Bloch waves of a provided ion configuration. If L_Bloch is set to true, the basis set will be used to expand electron Green's function of different ion configurations (See Eq. (12) of Ref. [3]).
```
>> Gen_BASIS_NDsym_nscf
```

* **3_EFERMI:** Calculates EFERMI_AV=\<Ef\> as an initial guess of the Fermi energy.
```
>> CAC_EFERMI
```       
        
* **4_Gbar:** Calculates electron Green's functions of different ion configurations, and take the average. Parallelization is excuted by Gbar_ND.m.
```
>> NCL=1; % Select the cluster number
>> Gbar_ND
```      

* **5_New_base:** Calculates the EPC-renormalized Bloch bases using the average Green's function according to Eq. (2) of Ref. [3].
```
>> New_basis
```        

* **6_Tbar:** Calculates the electron-electron pair scateering amplitude according to Eqs. (1) and (5) of Ref. [3].
```
>> NCL=1; % Select the cluster number
>> Tbar_ND
```   

* **7_Collect:** Solves the Bethe-Salpeter equation to calculate the effective electron-electron interactions, and solve linearized Eliashberg equations to solve Tc. The gaussian smearing for calculating EPC parameters are set to 0.02 *Ry* by default.
```
>> Collectq_NDsc_nscf
```   

* **8_Tools:** Contains tools to calculate mean squared displacements (MSD) and radial distribuction functions (RDF).

## Interface to VASP
The interface to VASP 5.4.4 is uploaded. The additional subroutines output necessary intermediate variables for re-constructing Hamiltoninians and overlap matrices. Green's functions are calculated using the program in /main_nscf. The inerface to i-pi is also encoded. To apply the patch, run the following in the top level directory of VASP:
```
patch -p0 < vasp.5.4.4-SPIA.patch
```
There are several new flags to control the outputs:
