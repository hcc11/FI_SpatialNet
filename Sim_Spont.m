% Simulate spontaneous dynamics with uncorrelated and untuned inputs from
% L4 (related to Figs. 5C, 7C)
% 
% mex EIF1DRFfastslowSyn2.c
% mex spktime2count.c

data_folder=''; % folder name to save data 

Wseed1_range=[48395]; %random seed for generating connectivity matrix   
Wseed2_range=[421762]; 
nws=1; 
%%%%%%%%%%%%% this part is for running as a job array on cluster %%%%%%%%%%%
rng('shuffle');
AI = getenv('SLURM_ARRAY_TASK_ID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);
% job_dex=1;
%%%%%%%%%%%%%%%%%%%%%
Ntrial=5; % 5 trials for each parameter set
Np=9;
nt=mod(job_dex-1,Ntrial)+1+Ntrial*(ceil(job_dex/Ntrial/Np)-1),
pid= mod(ceil(job_dex/Ntrial)-1,Np)+1,

if pid<=3
    nType=pid;
    Jx=[240; 400];
    Sigs=[0.1; 0.2; 0.3];
    sigmaRF=0.05;sigmaE=0.1; sigmaI=Sigs(nType);
    Iapp=[0;0]; 
    Type=strrep(sprintf('L1_sigmaRF%.03g_sigmaE%.03g_sigmaI%.03g',sigmaRF,sigmaE,sigmaI),'.','d')
    param(1).sigmaRX=sigmaRF*ones(2,1);
    param(1).sigmaRR=[sigmaE sigmaI; sigmaE sigmaI];
    opt.useWfile=1;
    W_fname=sprintf('%sweight%s_v2',data_folder,Type),
    
elseif pid<=9
    nType=pid-3;
    Jx=[140; 0];
    muI_range=[0  .4  .8  1.2  1.6  2];
    Iapp=[0; muI_range(nType)];
    sigmaRF=0.05;sigmaRR=0.1;
    Type=strrep(sprintf('L1_sigmaRF%.03g_sigmaRR%.03g_Jix%.03g_muI%.03g',sigmaRF,sigmaRR,Jx(2),Iapp(2)),'.','d')
    param(1).sigmaRX=sigmaRF*ones(2,1);
    param(1).sigmaRR=sigmaRR*ones(2,2);
    opt.useWfile=1;
    W_fname=sprintf('%sweight%s_v2',data_folder,'L1_sigmaRF0d05_sigmaRR0d1'),
end


for trial=nt
rng(trial + seed_offset);

dt=.05;

opt.save=1; % save data 
opt.saveS1=1; 
opt.CompCorr=0; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 2 & 1  
opt.loadS1=0; 
opt.plotPopR=0; % plot population rate
opt.fixW=1; 
    Wseed1=Wseed1_range(nws);
    Wseed2=Wseed2_range(nws); 
opt.saveW=0;
opt.savecurrent=0; 
opt.saveRm=0;
opt.Layer1only=1;
if trial==1
    opt.saveParam=1; 
else
    opt.saveParam=0; 
end

T=20000;
p_stim.stim_type='Uncorr';
p_stim.Nstim=1;
p_stim.Nsource=1;
p_stim.rX=.01; 
p_stim.Nsource=1;

filename=strrep(sprintf('%sRF2D3layer_Spont_%s_dt%.03g_ID%.0f',...
    data_folder,Type,dt,trial),'.','d'),


ParamChange={ 'filename', filename;'dt', dt; 'T', T; 'Nc',Nc; 'p_stim',p_stim;...
    'param(1).sigmaRX',param(1).sigmaRX; 'param(1).sigmaRR',param(1).sigmaRR };
if exist('s1_fname','var')
    ParamChange=[ParamChange;{'s1_fname',s1_fname}];
end

ParamChange=[ParamChange;{'param(1).Jx', Jx; 'param(1).Iapp', Iapp}]; 

if opt.fixW && opt.useWfile==0
    ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
end

if opt.useWfile==1
    if exist([W_fname '.mat'], 'file')==0  % if W_fname does not exist, save W in simulation 
        opt.saveW=1;
    end
end
if opt.saveW==1||opt.useWfile==1
    ParamChange=[ParamChange;{'W_fname',W_fname}];
end
   

clear F Filter theta_map X Y Imag;

RF2D3layer(opt, ParamChange)

save(filename,'Jx','Iapp','ParamChange','-append')

end

