function CollectSpk(Task,RunPart)
    %%%%%%%%%%  parameter sets  %%%%%%%%%%%%%%%%%
    
%     %%%%  alpha_ffwd %%%%%e
%     Types={'L1_sigmaRF0d05_sigmaRR0d1','L1_sigmaRF0d0625_sigmaRR0d125','L1_sigmaRF0d075_sigmaRR0d15',...
%         'L1_sigmaRF0d1_sigmaRR0d2','L1_sigmaRF10_sigmaRR10'};
%     
%     %%%%%%% alpha_rec  %%%%%%%%
%     Types={'L1_sigmaRF0d05_sigmaRR0d05','L1_sigmaRF0d0625_sigmaRR0d15','L1_sigmaRF0d05_sigmaRR0d2',...
%         'L1_sigmaRF0d075_sigmaRR0d05','L1_sigmaRF0d075_sigmaRR0d15','L1_sigmaRF0d075_sigmaRR0d2'};
%     
%     %%%%%%% sigma_i %%%%%%%%%%%
%     Types={'L1_sigmaRF0d05_sigmaRR0d1','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d2','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3'};
%     
%     %%%%  mu_i %%%%%e
%     Types={'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d8',...
%         'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d4','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI1d2'};
%     
%     %%%%%%% multi layer %%%%%%%%%%%%%%
%     Types={'multiL_L1_muI0_','multiL_L1_muI1d2_','multiL_L2_muI0_','multiL_L2_muI1d2_',...
%         'multiL_L3_muI0_','multiL_L3_muI1d2_'};

%%%%%%%%%%%%%%%%%% sigI = 0.1, 0.3, P_ts = 0.1, 0.2 %%%%%%%%%%%%%%%%%%
%  Types={'L1_sigmaRF0d05_sigmaRR0d1_Pts0d1',...
%             'L1_sigmaRF0d05_sigmaRR0d1_Pts0d2',...
%             'L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3_Pts0d1',...
%             'L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3_Pts0d2'}; 
    
%%%%%%% multi layer w/ tuning dependent connections %%%%%%%%%%%%%%
    Types={'multiL_L1_muI0_Pts0d2','multiL_L1_muI1d2_Pts0d2',...
        'multiL_L2_muI0_Pts0d2','multiL_L2_muI1d2_Pts0d2',...
        'multiL_L3_muI0_Pts0d2','multiL_L3_muI1d2_Pts0d2'};

    data_folder='';
    sigma_n=3.5;
    
    if RunPart==1
        switch Task
            %% collect spkcounts from 500 trials into a file w/ name 'fnamesave' 
            case 'L1'
                AI = getenv('SLURM_ARRAY_TASK_ID');
                job_dex = str2num(AI);
                
                Pop='E1';
%                 Pop='I1';
                Np=length(Types);
                
                ntype=mod(job_dex-1,Np)+1;
                num=ceil(job_dex/Np);                
                dt=0.05;
                
                type=Types{ntype},
                
                if strcmp(type(1:6),'Tuning')
                    datafname=@(ID) strrep(sprintf('%sRF2D3layer_GaborTheta_sigma_n%.03g_%s_dt%.03g_ID%.0f',...
                        data_folder,sigma_n,type,dt,ID),'.','d');
                else
                    datafname=@(ID) strrep(sprintf('%sRF2D3layer_GaborTheta_Tseg_sigma_n%.03g_%s_dt%.03g_ID%.0f',...
                        data_folder,sigma_n,type,dt,ID),'.','d');
                end
                
                switch Pop
                    case 'E1'
                        fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_%ssum_%d',...
                            data_folder,sigma_n,type,num),'.','d'),
                        N=4e4;
                    case 'I1'
                        fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_%s_%ssum_%d',...
                            data_folder,sigma_n,type,Pop,num),'.','d'),
                        N=1e4;
                    case 'X'
                        fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_%s_%ssum_%d',...
                            data_folder,sigma_n,type,Pop,num),'.','d'),
                        N=2500;
                end
                
                Nstim=39;
                Ntrial=500;
                ns=Nstim*Ntrial;
                
                %%%%%%%%%%%%%% check if files are complete %%%%%%%%%%%%%%%%%%%%%%
                IDval=[];
                for k=1:Ntrial
                    ID=k+(num-1)*Ntrial;
                    if exist([datafname(ID) '.mat'], 'file')==0
                        sprintf('file %s, does not exist\n',datafname(ID))
                    else
                        listOfVariables = who('-file', datafname(ID));
                        if ismember('th_id', listOfVariables)==0||ismember('E1', listOfVariables)==0
                            sprintf('file ID %s, variables not complete\n',datafname(ID))
                        else
                            IDval=[IDval k];
                        end
                    end
                end
                setdiff(1:Ntrial, IDval),
                if nnz(IDval)<Ntrial
                    error('files not complete')
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                X=zeros(N,ns,'int8');
                th_id=zeros(ns,1,'int8');
                
                load(datafname(1),'Tw','T','p_stim')
                Nt=floor(T/Tw);
                Nt_on=p_stim.T_on/Tw;
                Nt_off=p_stim.T_off/Tw;
                Nseg=Nt_on+Nt_off;
                idt=1:Nt;
                idt(mod(idt-1,Nseg)+1<=Nt_off)=0;
                idt=idt(idt>0);
                nnz(idt),
                tic
                for k=1:Ntrial
                    ID=k+(num-1)*Ntrial;
                    data=load(datafname(ID));
                    th_id((1:Nstim)+(k-1)*Nstim)=data.th_id(2:end);
                    switch Pop
                        case 'E1'
                            sc=permute(sum(reshape(data.E1(:,idt),N,Nt_on,[]),2),[1 3 2]);
                        case 'I1'
                            sc=permute(sum(reshape(data.I1(:,idt),N,Nt_on,[]),2),[1 3 2]);
                        case 'X'
                            sc=data.Xstim;
                    end
                    
                    X(:,(1:Nstim)+(k-1)*Nstim)=sc(:,2:end);
                    
                end
                size(sc)
                toc
                save(fnamesave,'X','th_id','Ntrial','datafname')
                
                %   save parameters
                if num==1
                    load(datafname(1), 'T','param','p_stim','testp','Jx')
                    save(fnamesave,'T','param','p_stim','testp','Jx','-append')
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
            case 'multiL' 
                %% collect spkcounts from multi-layer network simulations  
                muI_range =[0    1.2];
                P_ts = 0.2; 
                
                AI = getenv('SLURM_ARRAY_TASK_ID');
                job_dex = str2num(AI);
                Np=2;
                pid=mod(job_dex-1,Np)+1;
                num=mod(ceil(job_dex/Np)-1,6)+1,
                nLayer=ceil(job_dex/12),  % nLayer=1, 2, 3
                
                muI=muI_range(pid);
%                 Pop='Inh';
                Pop='Exc';
                switch Pop
                    case 'Exc'
                        if P_ts
                            fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_multiL_L%d_muI%.03g_Pts%.03gsum_%d',...
                                data_folder,sigma_n,nLayer,muI,P_ts,num),'.','d'),
                        else
                            fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_multiL_L%d_muI%.03g_sum_%d',...
                                data_folder,sigma_n,nLayer,muI,num),'.','d'),
                        end
                        N=4e4;
                    case 'Inh'
                        if P_ts
                            fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_multiL_L%d_muI%.03g_Pts%.03g_%ssum_%d',...
                                data_folder,sigma_n,nLayer,muI,P_ts,Pop,num),'.','d'),
                        else
                            fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_multiL_L%d_muI%.03g_%ssum_%d',...
                                data_folder,sigma_n,nLayer,muI,Pop,num),'.','d'),
                        end
                        N=1e4;
                end
                dt=0.05;
                Ntshift=0;
                
                if nLayer==1
                    nLayer2=2;
                else
                    nLayer2=nLayer;
                end
                if P_ts 
                    datafname=@(ID) strrep(sprintf('%sRF2D3layer_L%d_GaborTheta_Tseg_sigma_n%.03g_muI%.03g_Pts%.03g_dt%.03g_ID%.0f',...
                        data_folder2,nLayer2,sigma_n,muI,P_ts,dt,ID),'.','d'),
                else
                    datafname=@(ID) strrep(sprintf('%sRF2D3layer_L%d_GaborTheta_Tseg_sigma_n%.03g_muI%.03g_dt%.03g_ID%.0f',...
                        data_folder2,nLayer2,sigma_n,muI,dt,ID),'.','d'),
                end
                
                Ntrial=500;
                Nstim=39;
                ns=Nstim*Ntrial;
                %%%%%%%%%%%%%% check files %%%%%%%%%%%%%%%%%%%%%%
                IDval=[];
                for k=1:Ntrial
                    ID=k+(num-1)*Ntrial;
                    if exist([datafname(ID) '.mat'], 'file')==0
                        sprintf('file %s, does not exist\n',datafname(ID))
                    else
                        listOfVariables = who('-file', datafname(ID));
                        if ismember('th_id', listOfVariables)==0||ismember('E2', listOfVariables)==0
                            sprintf('file ID %s, variables not complete\n',datafname(ID))
                        else
                            IDval=[IDval k];
                        end
                    end
                end
                setdiff(1:Ntrial, IDval),
              
                if nnz(IDval)<Ntrial
                    error('files not complete')
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                X=zeros(N,ns,'int8');
                th_id=zeros(ns,1,'int8');
                
                load(datafname(1),'T')
                Tw=50;
                p_stim.T_on=200;
                p_stim.T_off=300;
                Nt=floor(T/Tw);
                Nt_on=p_stim.T_on/Tw;
                Nt_off=p_stim.T_off/Tw;
                Nseg=Nt_on+Nt_off;
                idt=1:Nt;
                idt(mod(idt-1,Nseg)+1<=Nt_off)=0;
                if Ntshift
                    idt=idt(Nt_off+1:end-Nt_on);
                end
                idt=idt(idt>0);
                idt=idt+Ntshift;
                
                nnz(idt)
                
                tic
                for k=1:Ntrial
                    ID=k+(num-1)*Ntrial;
                    data=load(datafname(ID));
                    if Ntshift
                        th_id((1:Nstim)+(k-1)*Nstim)=data.th_id(2:end-1);
                    else
                        th_id((1:Nstim)+(k-1)*Nstim)=data.th_id(2:end);
                    end
                    switch Pop
                        case 'Exc'
                            if nLayer==1
                                sc=permute(sum(reshape(data.E1(:,idt),N,Nt_on,[]),2),[1 3 2]);
                            else
                                sc=permute(sum(reshape(data.E2{1}(:,idt),N,Nt_on,[]),2),[1 3 2]);
                            end
                        case 'Inh'
                            if nLayer==1
                                sc=permute(sum(reshape(data.I1(:,idt),N,Nt_on,[]),2),[1 3 2]);
                            else
                                sc=permute(sum(reshape(data.I2{1}(:,idt),N,Nt_on,[]),2),[1 3 2]);
                            end
                    end
                    
                    X(:,(1:Nstim)+(k-1)*Nstim)=sc(:,2:end);
                end
                size(sc)
                toc
                save(fnamesave,'X','th_id','Ntrial','Ntshift','datafname')
                if num==1
                    load(datafname(1), 'T','param','p_stim','muI','W_fname')
                    save(fnamesave,'T','param','p_stim','muI','W_fname','-append')
                end
                
        end
        
    elseif RunPart==2
        %% compute mean and variance of spike counts 
        rng('shuffle');
        AI = getenv('SLURM_ARRAY_TASK_ID');
        job_dex = str2num(AI);
        seed_offset = randi(floor(intmax/10));
        rng(job_dex + seed_offset);
        
        type=Types{job_dex},
        fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_%s_fm',...
            data_folder,sigma_n,type),'.','d'),
        datafname=@(ID) strrep(sprintf('%sGaborTheta_sigma_n%.03g_%ssum_%d',...
            data_folder,sigma_n,type,ID),'.','d');
        
        Ntrial=1;  % one file to compute Fm and Var
        
        data=load(datafname(1));
        Nstim=size(data.X,2);
        N=size(data.X,1);
        
        ns=Nstim*Ntrial;
        Fm=zeros(N,2); % mean of spike counts, each column is for one orientation  
        Var=zeros(N,2); % variance of spike counts, each column is for one orientation 
        
        X=zeros(ns,N);
        th=zeros(ns,1);
        
        for ID=1:Ntrial
            data=load(datafname(ID));
            th((1:Nstim)+(ID-1)*Nstim)=data.th_id;
            X((1:Nstim)+(ID-1)*Nstim,:)=data.X';
        end
        
        theta0=unique(th),
        
        Fm(:,1)=mean(X(th==theta0(1),:),1);
        Fm(:,2)=mean(X(th==theta0(2),:),1);
        
        Var(:,1)=var(X(th==theta0(1),:),[],1);
        Var(:,2)=var(X(th==theta0(2),:),[],1);
        
        save(fnamesave,'Fm','Var','type','datafname','Ntrial')
    elseif RunPart==3 % collect Fm & Var from different parameter sets
        
        Np=length(Types);       
        fnamesave=strrep(sprintf('%sGaborTheta_sigma_n%.03g_fm_sum_Pts',...
            data_folder,sigma_n),'.','d'),
        
        for pid=1:Np
            type=Types{pid};
            fname=strrep(sprintf('%sGaborTheta_sigma_n%.03g_%s_fm',...
                data_folder,sigma_n,type),'.','d'),
            data=load(fname);
            res(pid).Fm=data.Fm;
            res(pid).Var=data.Var;
            res(pid).Ntrial=data.Ntrial;
            res(pid).datafname=data.datafname;
        end
        save(fnamesave,'Types','res')       
        
    end
