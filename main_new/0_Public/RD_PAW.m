function [NMAX,RG,SIMPI,LLMAX,WAE_PS,WAE,WPS,...
    LMAX,LMMAXC,CH0,CH1,CH2,CH3,LPS]=RD_PAW(ROOT_DIR,NTYP)
NMAX=cell(1,NTYP);RG=cell(1,NTYP);SIMPI=cell(1,NTYP);LLMAX=cell(1,NTYP);
WAE_PS=cell(1,NTYP);WAE=cell(1,NTYP);WPS=cell(1,NTYP);
for NT=1:NTYP
    FID=fopen([ROOT_DIR,'/BEAD_1/PROJ_1_1_1_',int2str(NT)]);
    LLMAX{NT}=fread(FID,1,'int');
    fclose(FID);
    
    FID=fopen([ROOT_DIR,'/BEAD_1/SIMPI_1_',int2str(NT)]);
    % Read in information needed for transformation matrix T
    % Radial grids and coefficients of Simpson's integration
    NMAX{NT}=fread(FID,1,'int');fread(FID,1,'int');
    RG{NT}=fread(FID,NMAX{NT},'double');
    SIMPI{NT}=fread(FID,NMAX{NT},'double');
    fclose(FID);
    
    FID=fopen([ROOT_DIR,'/BEAD_1/PWAV_1_',int2str(NT)]);
    % ae partial wave minus ps partial wave
    % i.e. r^2(|phi_i^a>-|\tilde{phi}_i^a>)on radial grids
    WAE_PS{NT}=fread(FID,[NMAX{NT},LLMAX{NT}],'double');
    fclose(FID);
    
    FID=fopen([ROOT_DIR,'/BEAD_1/WAE_1_',int2str(NT)]);
    % ae partial wave on radial grids
    WAE{NT}=fread(FID,[NMAX{NT},LLMAX{NT}],'double');
    fclose(FID);
    
    FID=fopen([ROOT_DIR,'/BEAD_1/WPS_1_',int2str(NT)]);
    % ps partial wave on radial grids
    WPS{NT}=fread(FID,[NMAX{NT},LLMAX{NT}],'double');
    fclose(FID);
end

%%
LMAX=zeros(1,NTYP);LMMAXC=zeros(1,NTYP);
CH0=cell(1,NTYP);CH1=cell(1,NTYP);CH2=cell(1,NTYP);CH3=cell(1,NTYP);
for NT=1:NTYP
    %% Calculate <\psi_\bm{k}|\phi_i^a> around each type of ion
    PROJJ=fopen([ROOT_DIR,'/BEAD_1/PROJ_1_1_1_',int2str(NT)]);
    LMAX(NT)=fread(PROJJ,1,'int');LMMAXC(NT)=fread(PROJJ,1,'int');
    CH0{NT}=fread(PROJJ,2,'int');CH1{NT}=fread(PROJJ,2,'int');
    CH2{NT}=fread(PROJJ,2,'int');CH3{NT}=fread(PROJJ,2,'int');
    fclose(PROJJ);
end

LPS = cell(1,NTYP);
for NT=1:NTYP
    %% Calculate <\psi_\bm{k}|\phi_i^a> around each type of ion
    LM=fopen([ROOT_DIR,'/BEAD_1/LM_INDEX_1_',int2str(NT)]);
    LPS{NT} = fread(LM,LMAX(NT),'int');
    fclose(PROJJ);
end

end