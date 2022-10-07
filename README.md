# Stochastic path-integral approach (SPIA)
Stochastic path-integral approach for determing superconducting transition temperature Tc of anharmonic systems.
The main theory can be found in the following references:

[1] Theory in liquids: Huiying Liu, Ying Yuan, Donghao Liu, Xin-Zheng Li, and Junren Shi. Phys. Rev. Research 2, 013340 (2020). 

[2] Implementations in liquids: Haoran Chen, Xiao-Wei Zhang, Xin-Zheng Li, and Junren Shi. Phys. Rev. B 104, 184516 (2021).

[3] Theory and implementations in liquids: Haoran Chen, Junren Shi. arXiv:2205.07247.

Instructions are not complete yey.

# Main program
The main program is in the directory SPIA/main_nscf. 
The program basically runs from 1_GEN_CORE to 7_Collect in steps:

0_Public: contains functions that are shared by multiple steps. 
Parameters of the calculations are the descriptions can be found in 0_Public/parameter.m.

1_GEN_CORE: calculate the interpolation table for Eqs. (15) and (A1) in Ref. [3]. 

2_GEN_BASIS: calculate the Bloch waves of a provided ion configuration. If L_Bloch is set to true, the basis set will be used to expand electron Green's function of different ion configurations (See Eq. (12) of Ref. [3]).

3_EFERMI: calculate EFERMI_av=<Ef> as an initial guess of the Fermi energy.

4_Gbar: calculate electron Green's functions of different ion configurations, and take the average. Parallelization is excuted by Gbar_ND.m.

5_New_base: calculate the EPC-renormalized Bloch bases using the average Green's function according to Eq. (2) of Ref. [3].

6_Tbar: calculate the electron-electron pair scateering amplitude according to Eqs. (1) and (5) of Ref. [3].

7_Collect: solve the Bethe-Salpeter equation to calculate the effective electron-electron interactions, and solve linearized Eliashberg equations to solve Tc.

8_Tools: tools to calculate mean squared displacements (MSD) and radial distribuction functions (RDF).

# Interface to VASP
