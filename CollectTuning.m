function CollectTuning(Task,RunPart)
% compute tuning curves
% First, run Sim_Ori_gabor_L1.m to simulate spike counts
% with
% testp.theta0=0.02:.02:1;
% Nth=length(testp.theta0);
% Ntrial=250;

switch Task
    case 'L1'
        %%%%%%%%%%  parameter sets  %%%%%%%%%%%%%%%%%
        
        %%%%  alpha_ffwd %%%%%e
        Types={'L1_sigmaRF0d05_sigmaRR0d1','L1_sigmaRF0d0625_sigmaRR0d125','L1_sigmaRF0d075_sigmaRR0d15',...
            'L1_sigmaRF0d1_sigmaRR0d2','L1_sigmaRF10_sigmaRR10'};
        
        %%%%%%% alpha_rec  %%%%%%%%
        Types={'L1_sigmaRF0d05_sigmaRR0d05','L1_sigmaRF0d0625_sigmaRR0d15','L1_sigmaRF0d05_sigmaRR0d2',...
            'L1_sigmaRF0d075_sigmaRR0d05','L1_sigmaRF0d075_sigmaRR0d15','L1_sigmaRF0d075_sigmaRR0d2'};
        
        %%%%%%% . sigma_i %%%%%%%%%%%
        Types={'L1_sigmaRF0d05_sigmaRR0d1','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d2','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3'};
        
        %%%%  mu_i %%%%%e
        Types={'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d8',...
            'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d4','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI1d2'};
        
        %%%%%%% multi layer %%%%%%%%%%%%%%
        Types={'multiL_L1_muI0_','multiL_L1_muI1d2_','multiL_L2_muI0_','multiL_L2_muI1d2_',...
            'multiL_L3_muI0_','multiL_L3_muI1d2_'};
        
        data_folder=''; % folder name to save data
        
        sigma_n=3.5;
        Np=length(Types);
        
        dt=0.05;
        fnamesave0=@(type) strrep(sprintf('%sRF2D3layer_GaborTheta_sigma_n%.03g_%s_Tuning',...
            data_folder,sigma_n,type),'.','d'),  % filename to save each tuning curves
        
        fnamesave=strrep(sprintf('%sRF2D3layer_GaborTheta_sigma_n%.03g_L1_Tuning_sum',...
            data_folder,sigma_n),'.','d'), % filename to save a collection of tuning curves
        
        filename=@(type,ID) strrep(sprintf('%sRF2D3layer_GaborTheta_sigma_n%.03g_Tuning%s_dt%.03g_ID%.0f',...
            data_folder,sigma_n,type,dt,ID),'.','d');  % filename of simulations
        
        %% Tuning L1
        if RunPart==1
            rng('shuffle');
            %%%%%%%%%%%%%%%%%%%%%
            AI = getenv('SLURM_ARRAY_TASK_ID');
            job_dex = str2num(AI);
            seed_offset = randi(floor(intmax/10));
            rng(job_dex + seed_offset);
            %%%%%%%%%%%%%%%%%%%%%
            nType=job_dex,
            type=Types{nType};
            Ntrial=250;
            
            Ne=4e4;
            Ni=1e4;
            load(filename(type,1),'Tw','T','p_stim','testp')
            Nt=floor(T/Tw);
            Nt_on=p_stim.T_on/Tw;
            Nt_off=p_stim.T_off/Tw;
            Nseg=Nt_on+Nt_off;
            idt=1:Nt;
            idt(mod(idt-1,Nseg)+1<=Nt_off)=0;
            Nstim=ceil(T/(p_stim.T_on+p_stim.T_off))-1,
            sc1=zeros(Ne,Nstim*Ntrial);
            sc2=zeros(Ni,Nstim*Ntrial);
            th=zeros(Nstim*Ntrial,1);
            %%%%%%%%%%%%%% check files %%%%%%%%%%%%%%%%%%%%%%
            %                 IDval=[];
            %                 for k=1:Ntrial
            %                     if exist([filename(type,k) '.mat'], 'file')==0
            %                         sprintf('file %s, does not exist\n',filename(type,k))
            %                     else
            %                         listOfVariables = who('-file', filename(type,k));
            %                         if ismember('th_id', listOfVariables)==0||ismember('E1', listOfVariables)==0
            %                             sprintf('file ID %s, variables not complete\n',filename(type,k))
            %                         else
            %                             IDval=[IDval k];
            %                         end
            %                     end
            %                 end
            %                 if nnz(IDval)<Ntrial
            %                     error('files not complete')
            %                 end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k=1:Ntrial
                data=load(filename(type,k),'th_id','E1','I1');
                th((1:Nstim)+(k-1)*Nstim)=data.th_id(2:end);
                Y=permute(sum(reshape(data.E1(:,idt>.5),Ne,Nt_on,[]),2),[1 3 2]);
                sc1(:,(1:Nstim)+(k-1)*Nstim)=Y(:,2:end);
                Y=permute(sum(reshape(data.I1(:,idt>.5),Ne,Nt_on,[]),2),[1 3 2]);
                sc2(:,(1:Nstim)+(k-1)*Nstim)=Y(:,2:end);
            end
            
            disp('collected spk counts')
            
            theta=unique(th);
            Nth=numel(theta);
            E1_th=zeros(Ne,Nth);
            I1_th=zeros(Ni,Nth);
            
            for k=1:Nth
                idx=th<k+.1 & th>k-.1;
                E1_th(:,k)=mean(sc1(:,idx),2);
                I1_th(:,k)=mean(sc2(:,idx),2);
            end
            
            save(fnamesave0(type),'theta','E1_th','I1_th','filename','testp')
            
        elseif RunPart==2
            for pid=1:Np
                type=Types{pid};
                data=load(fnamesave0(type));
                res(pid).E1_th=data.E1_th;
                res(pid).I1_th=data.I1_th;
                res(pid).theta=data.theta;
                res(pid).testp=data.testp;
            end
            save(fnamesave,'res','Types','fnamesave0')
        end
end

end
