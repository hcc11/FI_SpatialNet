% To compute the linear Fisher information and correlations, modified codes
% from 
% I. Kanitscheider, R. Coen-Cagli, A. Kohn, A. Pouget, 
% Measuring fisher information accurately in correlated neural populations. 
% PLoS computational biology 11, e1004218 (2015).

% outputs: 
% FI_BC: the bias-corrected Fisher information 
% FIVAL: the Fisher information on the Validation set using early
% stopping method 
% FITR: the Fisher information on the Training set using early
% stopping method  
% (FI_BC is typically between FITR and FIVAL, and converge as # of trials
% increase.  Can set NR=0 to not run early stopping method, which can be slow for large N. )
% Nid: sampled neuron index for computing FI  
% compute and save correlation vs distance when N=1600

% need to define: 
% datafname: data filename for spike count matrix (#neurons x #trials) 
% fnamesave: filename to save computed FI 
% data_folder: folder name to read spike count data and to save FI 


data_folder='data/'; % folder name to read spike count data and to save FI 

%%%%%%%%% to use on cluster %%%%%%%%%
rng('shuffle');
AI = getenv('SLURM_ARRAY_TASK_ID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% job_dex=1; % for indexing the parameter set, # of neurons and # of repeats

%%%%%%%%%%%%%%% param sets %%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%  alpha_ffwd %%%%%e
% Types={'L1_sigmaRF0d05_sigmaRR0d1','L1_sigmaRF0d0625_sigmaRR0d125','L1_sigmaRF0d075_sigmaRR0d15',...
%     'L1_sigmaRF0d1_sigmaRR0d2','L1_sigmaRF10_sigmaRR10'};
% 
% %%%%%%% alpha_rec  %%%%%%%%
% Types={'L1_sigmaRF0d05_sigmaRR0d05','L1_sigmaRF0d0625_sigmaRR0d15','L1_sigmaRF0d05_sigmaRR0d2',...
%     'L1_sigmaRF0d075_sigmaRR0d05','L1_sigmaRF0d075_sigmaRR0d15','L1_sigmaRF0d075_sigmaRR0d2'};
% 
% %%%%%%%  sigma_i %%%%%%%%%%%
% Types={'L1_sigmaRF0d05_sigmaRR0d1','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d2','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3'};
% 
% %%%%  mu_i %%%%%e
% Types={'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d4',...
%     'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d8','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI1d2'};
% 
% %%%%%%% multi layer %%%%%%%%%%%%%%
% Types={'multiL_L1_muI0_','multiL_L1_muI1d2_','multiL_L2_muI0_','multiL_L2_muI1d2_',...
%     'multiL_L3_muI0_','multiL_L3_muI1d2_'};
%
% %%%%%%%%% tuning dependent connection %%%%%%%%
% Types={'L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3_Pts0d1',...
%     'L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3_Pts0d2',...
%     'L1_sigmaRF0d05_sigmaRR0d1_Pts0d1',...
%     'L1_sigmaRF0d05_sigmaRR0d1_Pts0d2'};
% 
% Types={'multiL_L1_muI0_Pts0d2','multiL_L1_muI1d2_Pts0d2',...
%     'multiL_L2_muI0_Pts0d2','multiL_L2_muI1d2_Pts0d2',...
%     'multiL_L3_muI0_Pts0d2','multiL_L3_muI1d2_Pts0d2'};

Nfile=6; % number of spike count files per parameter set 

ds=0.01;
sigma_n=3.5; % noise intensity 
Np=length(Types),  % # of parameter sets to run 
task='N'; % compute FI for different number of neurons 
% task='Ntr'; % compute FI for different number of trials
switch task
    case 'N'
        Nn=5; % # of N to compute
        Nrep=20; % # of repetition of neuron sampling  
        ntype=ceil(job_dex/Nrep/Nn); % index for parameter set in 'Types'
        Nrun=mod(job_dex-1,Nrep)+1; % index for the repetition number 
        ipN= mod(ceil(job_dex/Nrep)-1,Nn)+1; % index for N in 'N_range'
    case 'Ntr'
        Nn=10; % # of Ntr to compute 
        Nrep=10; % # of repetition of neuron sampling 
        ntype=ceil(job_dex/Nrep/Nn); % index for parameter set in 'Types'
        Nrun=mod(job_dex-1,Nrep)+1; % index for the repetition number 
        ipN=mod(ceil(job_dex/Nrep)-1,Nn)+1; % index for Ntr in 'Ntr_range'
end

Type=Types{ntype}, 

datafname=@(ID) strrep(sprintf('%sGaborTheta_sigma_n%.03g_%ssum_%d',...
    data_folder,sigma_n,Type,ID),'.','d'); % data filename for spike count matrix (#neurons x #trials) 

data=load(datafname(1));
Nstim=size(data.X,2);
ns=Nstim*Nfile;

%%%%%%%%%%%% select neurons with FR > 1 Hz %%%%%%%%%%%%%%
Tw=0.2;
FR_th=1;
Fm = (mean(data.X(:,data.th_id==1),2)+mean(data.X(:,data.th_id==2),2))/2; % average spike count per neuron
ind_FR=find(Fm>FR_th*Tw);
sprintf('number of neurons w/ rate larger than %d Hz: ',FR_th, nnz(ind_FR)) 
%%%%%%%%% all neurons %%%%%%%%%%%%%%%%%%
% ind_FR=(1:4e4)';

switch task
    case 'N'
        %%%%%%%%%%%% vary N %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N_range=[50 100 200 400 800 1600 3200 6400  12800  4e4];
        N=N_range(ipN),
        fnamesave=strrep(sprintf('%sGD_GaborTheta_sigma_n%.03g_%s_N%d_%d',...
            data_folder2,sigma_n,Type,N,Nrun),'.','d'),
        Ntr=ns/2;
        if N==4e4
            Nid=(1:4e4)';
        else
            Nid=randsample(ind_FR, min(N,nnz(ind_FR))); % sampled neuron index for computing FI 
        end
    case 'Ntr'
        %%%%%%%%%%%% vary Ntr %%%%%%%%%%%%%%%%%%%%%%%%%%%
        N=6400; % # of neuron 
        Ntr_range=round(exp(linspace(log(N*2),log(ns/2),10))); % range for trial #   
        Ntr=Ntr_range(ipN),
        fnamesave=strrep(sprintf('%sGD_GaborTheta_sigma_n%.03g_%s_N%d_Ntr%d_%d',...
            data_folder2,sigma_n,Type,N,Ntr,Nrun),'.','d'),
        if ipN==1
            Nid=randsample(ind_FR, N);
        else
            load(strrep(sprintf('%sGD_GaborTheta_sigma_n%.03g_%s_N%d_Ntr%d_%d',...
                data_folder,sigma_n,Type,N,Ntr_range(1),Nrun),'.','d'),'Nid')
        end     
end
NR=5;  % how many cross-validation splits for early stopping
if N>5e3
    NR=0; % do not run early stopping for N>5000 
end

if exist([fnamesave '.mat'], 'file')
    delete([fnamesave '.mat'])
end

if N==4e4
    FIdecoder(datafname,fnamesave,Nid,NR,Nfile,Ntr,ds);
else
    COV=FIdecoder(datafname,fnamesave,Nid,NR,Nfile,Ntr,ds);
end

if N==1600 % compute correlation vs distance 
    U=triu(ones(size(COV)),1);
    covm=mean(COV(U==1)) % average covariance 
    R=corrcov(COV);
    corr=mean(R(U==1)) % average correlations 
    
    Ne1=200;
    Ix2=(ceil(Nid/Ne1))/Ne1;
    Iy2=(mod((Nid-1),Ne1)+1)/Ne1;
    D = pdist2([Ix2,Iy2],[Ix2,Iy2],'euclidean');
    D=D(U==1);
    R=R(U==1);
    dmax=0.5;
    dd=0.025;
    daxis=0:dd:dmax;
    [n, ind] = histc(D,daxis);
    n=n(1:end-1); % discard last bin (d>dmax)
    Cd=zeros(length(n),1); % correlatio vs distance 
    cov_m=zeros(length(n),1);
    for k=1:length(n)
        Cd(k)=mean(R(ind==k));
    end
    daxis=daxis(1:end-1)+dd/2; % distance axis 
    save(fnamesave,'covm','corr','Cd','daxis','-append')
end
