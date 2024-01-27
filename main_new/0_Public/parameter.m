%% Parameter
%% Constants
kB=8.61733034E-5; % Boltzmann constant
RYTOEV=13.605826;
AUTOA=0.529177249;
HSQDTM=RYTOEV*AUTOA*AUTOA; % hbar^2/2m_e
CITPI=2*pi*1i;

%% Machine info
NGPU=4;
NCLUSTER = 10; % How many  

%% Simulation parameters
ROOT_DIR='../../'; % Directory where PIMD is performed
TEMP=240; % Temperature
NCONF=6040; % Configuration number
NSTART=40; % Ignore pre-equilibrated steps
NBEAD=16; % Bead number (Trotter number)
NSKIP=40; % Number of steps between calculations

ENCUT=350; % Plane-wave energy cutoff

%% What is the system?
L_Liquid = false;
% if true, simply calculate all things under plane waves

%% Methods for calculating overlap matrix
L_OVERL_OLD=false; % should not be changed

L_OVERL_CORE=false; 
% When set to true, the core function <\phi_i^a-tiled{\phi}_i^a|\phi_j^b-\tilde{\phi}_j^b> will be calculated.
% The result will be stored in CQIJ3.mat.
% Must be set true in the very first calculation. 
% In later calculations in systems with the same elements, 
% one can copy CQIJ3.mat to ../ , and set L_OVERL_CORE to false.

%% How to determine the Fermi energy
I_Ef=3;
% I_Ef = 1: Approximate Ef as average over all configurations, Ef_Av=<Ef>.
% I_Ef = 2: Determine Ef from energy distribution of all configurations.
% I_Ef = 3: Same as 2, but use the nscf k-mesh (see N_k_nscf).
%           This require recalculation of all eigenenergies.
% Typically, I_Ef = 1 and 2 should yield similar results

%% Whether use provide positions to expand a common PS space
L_POS_new=true;
% If true, will read POSION_new.mat
% load([ROOT_DIR,'POSION_new.mat']

L_CAC_POS_new=true;
% When set to true, calculate the ion density distribution, 
% and extract the density peaks as the new positions of ions.
% The result will be stored in POSION_new.mat.
% Only set to false if one is performing calculations in a continuation job
% i.e., a system that has already generated POSION_new.mat previously
if ~L_POS_new && L_CAC_POS_new
    disp('No need to calculate new positions when L_POS_new is false!')
    disp('Resetted L_CAC_POS_new to false!')
    L_CAC_POS_new=false;
end

% When calculating POSION_new, all configurations are read from XDATCAR_PI. A larger set of configurations can be used.
% NCONF_POS_new=NCONF_XDATCAR_PI; 
% Program directly reads NCONF from XDATCAR_PI
NSTART_POS_new=0;
NSKIP_POS_new=1;

%% Symmetry
L_IBZ=true; % Whether use only the irreducibe part of Brillouin Zone
L_SYM=true; % Whether we force the system to obey the symmetry of the BEAD_1_sym 
% L_SYM=true will also symmetrize POSION_new
if L_IBZ && L_POS_new
    L_SYM=true; % Force POSION_new to be of the same symmetry as BEAD_1_sym
end
% When L_IBZ=false and L_SYM=true, both G and Gamma will be symmetrized
% This will slightly increase the computational cost, but increase the statistical stability.
% But make sure that the system does have the symmetry you supposed.
% This may be checked using a large set of PIMD configurations.

%% Some information for calculating the Green's function
% NBANDS=NBANDS_mul(_Ef/g)*NELE/2
NBANDS_mul = 3^3; % Bands number to be included when calculating Gamma 
NBANDS_mul_Ef = 3; % Bands number to be calculated when calculating Ef
NBANDS_mul_g = 3; % Bands number to expand the instaneous Green's function

LBloch = false; % Whether Hamiltionians are expanded in terms of initial Bloch waves
if L_POS_new && LBloch
    disp('L_Bloch is not available when L_POS_new is set to true!')
    disp('Resetted L_Bloch to false!')
    LBloch = false;
end

if L_Liquid & (LBloch | L_POS_new | L_SYM | L_IBZ)
    disp('For a liquid, all-electron plane waves are used.')
    disp('And no average positions exist.')
    disp('Presently, we also do not apply symmetrization to force isotropy.')
    disp('Resetted L_Bloch, L_POS_new, L_SYM and L_IBZ to false!')
    LBloch = false
    L_POS_new = false
    L_SYM = false
    L_IBZ = false
end

KOMEGA=16; % omega_Ns in quasi-static approximation

%% Supercell shape and k-mesh
N_sc_base = [2 2 2]; % Base 'primitive' cell for constructing non-diagonal supercell
N_sc_ex = [1 0 0;0 1 0;0 0 1]; % Non-diagonal supercell shape
N_k = [2,2,2]; % of the base cell, self-consistent k-mesh
N_k_nscf = [4,4,4]; % of the base cell, non-self-consistent k-mesh

%% Perform calculations along high-symmetry paths to visualize electronic structures
L_path = false; % Whether along high symmetry path
if L_path
    N_sc_ex = N_sc_ex*diag(N_sc_base);
    N_sc_base = [1 1 1];
end

NK_line = [5 11 9 5 3]; % Number of k points along each path 
VK_hsym = ... % High-symmetry points
    [0.5  , 0.25 , 0.75 ;... % W
     0.5  , 0.5  , 0.5  ;... % L
     0.0  , 0.0  , 0.0  ;... % Gamma
     0.5  , 0.0  , 0.5  ;... % X
     0.5  , 0.25 , 0.75 ;... % W
     0.375, 0.375, 0.75];    % K
 
%% Morel-Anderson Pseudopotential
mustar = 0.10;
