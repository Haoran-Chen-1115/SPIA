%% Parameter
%% Constants
kB=8.61733034E-5; % Boltzmann constant
RYTOEV=13.605826;AUTOA=0.529177249;HSQDTM=RYTOEV*AUTOA*AUTOA; % hbar^2/2m_e
CITPI=2*pi*1i;

%% Machine info
NGPU=4;
NCLUSTER = 1; 

%% Methods for calculating overlap matrix
L_OVERL_OLD=false; % should not be changed

%% Simulation parameters
ROOT_DIR='../../'; % Directory where PIMD is performed
TEMP=450; % Temperature
NCONF=10000; % Configuration number
NSTART=2000; % Pre-equilibrated steps
NBEAD=1; % Bead number (Trotter number)
NSKIP=5; % Number of steps between calculations

ENCUT=300; % Plane-wave energy cutoff

%% Some information for calculating the Green's function
NELE_TYP=[6,11]; % valance electrons for each type
NBANDS_mul = 3^3; % Bands number to be included: NBANDS=NBANDS_mul*NELE/2

LBloch = false; % Whether Hamiltionians are expanded in terms of initial Bloch waves

KOMEGA=16; % omega_Ns in quasi-static approximation

%% Supercell shape and k-mesh
N_sc_base = [2 2 2]; % Base 'primitive' cell for constructing non-diagonal supercell
N_sc_ex = [1 0 0; 0 1 0; 0 0 1]; % Non-diagonal supercell shape
N_k = [2,2,2]; % of the 2*2*2 supercell, self-consistent k-mesh
N_k_nscf = [4,4,4]; % of the 2*2*2 supercell, non-self-consistent k-mesh

%% Morel-Anderson Pseudopotential
mustar = 0.12;
