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

clearvars -except NCL
Tbar_ND;
