addpath ./0_Public
parameter;

% First generate the overlap core function
if L_OVERL_CORE
    cd ./1_GEN_CORE
    if NCL==1
	GEN_CORE_s;
    else
	while 1
            if ~exist('../CQIJ3.mat','file')
	        pause(600)
            else
		status=1;
	        try
	            load('../CQIJ3.mat');
	        catch
	            warning('Saving CQIJ3.mat. Wait for another 1 minute.')
	            pause(60)
	            status=0;
	        end
		if status==1
		    break;
	        end
	    end
        end
    end
    cd ../
end

% Then determine the new positions if necessary
clearvars -except NCL
parameter;
if L_CAC_POS_new
    cd ./2_Dens
    if NCL==1
	find_position_new;
    else
	while 1
            if ~exist('../POSION_new.mat','file')
	        pause(600)
            else
		status=1;
	        try
	            load('../POSION_new.mat');
	        catch
	            warning('Saving POSION_new.mat. Wait for another 1 minute.')
	            pause(60)
	            status=0;
	        end
		if status==1
		    break;
	        end
	    end
        end
    end
    cd ../
end	

% Then the Fermi energy of the system
clearvars -except NCL
parameter;

% Always first simply average 
cd ./3_EFERMI
CAC_EFERMI;
cd ../

if I_Ef==2 % From all eigenenergies
elseif I_Ef==3 % From all eigenenergies, recaulated on a denser k-mesh
    cd ./4_Spec
    LCAC_Ef=true;
    Spec_main;
    if NCL==1
	while 1
            existing=false(1,NKPTS_sym);
            for NK=1:NKPTS_sym
		existing(NK)=exist(['../Spectral_',int2str(NK),'.mat'],'file')
	    end
	    if all(existing)
		CAC_EFERMI_r;
		break;
	    end
	end
    else
        while 1
            if ~exist('../EFERMI.mat','file')
                pause(600)
            else
                status=1;
                try
                    load('../EFERMI.mat');
                catch
                    warning('Saving EFERMI.mat. Wait for another 1 minutes.')
                    pause(600)
                    status=0;
                end
                if status==1
                    break;
                end
            end
        end
    end
end

clearvars -except NCL
parameter;
if ~L_path
    cd ./5_SPIA
    Gbar_ND;
    
    clearvars -except NCL L_path L_SYM L_IBZ L_Liquid
    if ~L_Liquid
        if L_SYM && ~L_IBZ
            if NCL==1
                while 1
                    existing=false(1,NKPTS_sym);
                    for NK=1:NKPTS_sym
                	    existing(NK)=exist(['../Gbar_OK__',int2str(NK),'.mat'],'file')
                    end
                    if all(existing)
                        for NK=1:NKPTS_sym
            		delete(['../Gbar_OK__',int2str(NK),'.mat'])
            	    end
                	    break;
                    end
                end
                New_basis_sym;
            else
                while 1
                    if ~exist('../EFERMI_new.mat','file')
                        pause(600)
                    else
                        status=1;
                        try
                            load('../EFERMI_new.mat');
                        catch
                            warning('Saving new basis. Wait for another 1 minutes.')
                            pause(600)
                            status=0;
                        end
                        if status==1
                            break;
                        end
                    end
                end
            end
        else
            New_basis;
        end
    end

    clearvars -except NCL L_path
    Tbar_ND;
    cd ../
else
end

if NCL==1 && ~L_Path
    cd ./6_Collect
    Collectq_NDsc_lor_sym;
end
