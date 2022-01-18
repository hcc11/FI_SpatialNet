function CollectData(task)
    % Collect Fisher information vs. number of neuron (N) data from different parameter sets
    
    data_folder='';  % folder name to save data
    switch task
        case 'alpha_ffwd'
            Types={'L1_sigmaRF0d05_sigmaRR0d1','L1_sigmaRF0d0625_sigmaRR0d125','L1_sigmaRF0d075_sigmaRR0d15',...
                'L1_sigmaRF0d1_sigmaRR0d2','L1_sigmaRF10_sigmaRR10'};
            
        case 'alpha_rec'
            Types={'L1_sigmaRF0d05_sigmaRR0d05','L1_sigmaRF0d0625_sigmaRR0d15','L1_sigmaRF0d05_sigmaRR0d2',...
                'L1_sigmaRF0d075_sigmaRR0d05','L1_sigmaRF0d075_sigmaRR0d15','L1_sigmaRF0d075_sigmaRR0d2'};
            
        case 'sigmaI'
            Types={'L1_sigmaRF0d05_sigmaRR0d1','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d2','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3'};
            
        case 'muI'
            Types={'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d4',...
                'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d8','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI1d2'};
            
        case 'multiL'
            Types={'multiL_L1_muI0_','multiL_L1_muI1d2_','multiL_L2_muI0_','multiL_L2_muI1d2_',...
                'multiL_L3_muI0_','multiL_L3_muI1d2_'};
            
        case 'P_ts'
            Types={'L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3_Pts0d1',...
                'L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3_Pts0d2',...
                'L1_sigmaRF0d05_sigmaRR0d1_Pts0d1',...
                'L1_sigmaRF0d05_sigmaRR0d1_Pts0d2',...
                'multiL_L1_muI0_Pts0d2','multiL_L1_muI1d2_Pts0d2',...
                'multiL_L2_muI0_Pts0d2','multiL_L2_muI1d2_Pts0d2',...
                'multiL_L3_muI0_Pts0d2','multiL_L3_muI1d2_Pts0d2'};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    sigma_n=3.5;
    % sigma_n=0;
    
    Np=length(Types);
    Nn=9; % # of N
    Nrun=20;  % # of samplings per N 
    fname=@(Type,N,run) strrep(sprintf('%sGD_GaborTheta_sigma_n%.03g_%s_N%d_%d',...
        data_folder,sigma_n,Type,N,run),'.','d');  %filename of FI data as defined in FIdecoder_cluster_L1.m
    
    N_range=[50 100 200 400 800 1600 3200  6400  12800 25600];
    
    fnamesave=strrep(sprintf('%sGD_GaborTheta_sigma_n%.03g_%ssum',...
        data_folder,sigma_n,task),'.','d'), %filename to save summarized data 
    
    NR=5;
    
    FITR=NaN(Nn,NR,Nrun,Np);
    FIVAL=NaN(Nn,NR,Nrun,Np);
    FIBC=NaN(Nn,Nrun,Np);
    FInaive=NaN(Nn,Nrun,Np);
    for pid=1:Np
        Type=Types{pid},
        for run=1:Nrun
            for ipN=1:Nn
                N=N_range(ipN);
                if exist([fname(Type,N,run) '.mat'], 'file')
                    data=load(fname(Type,N,run));
                    listOfVariables = who('-file', fname(Type,N,run));
                    if ismember('FITR0', listOfVariables)==0
                        %  sprintf('file %s, FITR0 does not exist\n',fname(Type,N,run))
                    else
                        FITR(ipN,:,run,pid)=data.FITR0;
                        FIVAL(ipN,:,run,pid)=data.FIVAL0;
                    end
                    if ismember('FI_BC', listOfVariables)==0
                        sprintf('file %s, FI_BC does not exist\n',fname(Type,N,run))
                    else
                        FIBC(ipN,run,pid)=data.FI_BC;
                        FInaive(ipN,run,pid)=data.FInaive;
                    end
                else
                    sprintf('file %s does not exist\n',fname(Type,N,run))
                end
            end
            if exist([fname(Type,1600,run) '.mat'], 'file')
                data=load(fname(Type,1600,run),'covm','corr','Cd','daxis','Nfile','Ntr');
                res(pid).covm(run)=data.covm;
                res(pid).corr(run)=data.corr;
                res(pid).Cd(:,run)=data.Cd;
                res(pid).Nfile=data.Nfile;
                res(pid).Ntr=data.Ntr;
                res(pid).daxis=data.daxis;
            end
        end
    end
    
    save(fnamesave,'FITR','FIVAL','FIBC','N_range','res','Types')
 
end

