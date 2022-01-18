% simulation code for Figs. 1-4, 6, 8
% Run the following the first-time use: 
% mex EIF1DRFfastslowSyn.c
% mex spktime2count.c

% Run the following to generate orientation preference map before simulation: 
% Nx=50;
% Kc=5*2*pi; % spatial frequency
% theta_map=ori_map(Nx,Kc);
% save([data_folder 'theta_map'],'theta_map');

% all updated parameters need to be saved in struct 'ParamChange' 

data_folder=''; % folder name to save data 

Wseed1_range=[48395]; %random seed for generating connectivity matrix   
Wseed2_range=[421762]; 
nws=1; 
%%%%%%%%%%%%% this part is for running as a job array on cluster %%%%%%%%%%%
rng('shuffle');
% AI = getenv('SLURM_ARRAY_TASK_ID');
% job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);

job_dex = 1; 
%%%%%%%%% two orientations (for computing Fisher Info) %%%%%%%%%%%%
Ntrial=100; % per param set
Np=16;  % number of parameter sets to test 
nt=mod(job_dex-1,Ntrial)+1+Ntrial*(ceil(job_dex/Ntrial/Np)-1), % trial number = (nt-1)*5+(1:5)
pid=mod(ceil(job_dex/Ntrial)-1,Np)+1,  % choose param set 

dtheta=0.01;
testp.theta0=[ 0 ];  % orientations to test
% testp.theta0=[.5-dtheta/2, .5+dtheta/2];  % orientations to test
Nth=length(testp.theta0); 

%%%%%% check tuning curve %%%%%%%%%%
% testp.theta0=0.02:.02:1; 
% Nth=length(testp.theta0); 
% Ntrial=50;  % job_dex: 1:50*Np
% Np=16; 
% nt=mod(job_dex-1,Ntrial)+1+Ntrial*(ceil(job_dex/Ntrial/Np)-1), % trial number = (nt-1)*5+(1:5)
% pid=mod(ceil(job_dex/Ntrial)-1,Np)+1,


%%%%%% define parameters %%%%%%%%%%

if pid<=10
    task=1
    nType=pid;
elseif pid-10<=2
    task=2
    nType=pid-10;
elseif pid-12<=3
    task=3,
    nType=pid-12;
end

switch task
    case 1 % Figs 1-4
        Sigs=[0.05   .1;   .0625  .125;   0.075  .15;    .1    .2;    10  10;... % Fig 3
               .05  .05;     .05   .15;     .05   .2;   .075  .05;  .075  .2 ]; % Fig 4
        sigmaRF=Sigs(nType,1); % ffwd width 
        sigmaRR=Sigs(nType,2); % recurrent width 
        Type=strrep(sprintf('L1_sigmaRF%.03g_sigmaRR%.03g',sigmaRF,sigmaRR),'.','d'), 
        param(1).sigmaRX=sigmaRF*ones(2,1);
        param(1).sigmaRR=sigmaRR*ones(2,2);
        Jx=[240; 400]; % feedforward connection strengths 
        Iapp=[0; 0]; % static current inputs 
        W_fname=sprintf('%sweight%s_v2',data_folder,Type),  %filename to save connectivity matrix 

    case 2 % Fig. 6
        Sigs=[0.2; 0.3];   
        sigmaRF=0.05;sigmaE=0.1; sigmaI=Sigs(nType);
        Type=strrep(sprintf('L1_sigmaRF%.03g_sigmaE%.03g_sigmaI%.03g',sigmaRF,sigmaE,sigmaI),'.','d')
        param(1).sigmaRX=sigmaRF*ones(2,1);
        param(1).sigmaRR=[sigmaE sigmaI; sigmaE sigmaI];
        Jx=[240; 400]; 
        Iapp=[0; 0]; 
        W_fname=sprintf('%sweight%s_v2',data_folder,Type),  

    case 3 % Fig. 8
        sigmaRF=0.05;sigmaRR=0.1;
        Jx=[140; 0];
        muI_range=[0  .4   0.8  1.2]; 
        Iapp=[0; muI_range(nType)];
        Type=strrep(sprintf('L1_sigmaRF%.03g_sigmaRR%.03g_Jix%.03g_muI%.03g',sigmaRF,sigmaRR,Jx(2),Iapp(2)),'.','d')       
        param(1).sigmaRX=sigmaRF*ones(2,1);
        param(1).sigmaRR=sigmaRR*ones(2,2);
        W_fname=sprintf('%sweight%s_v2',data_folder,'L1_sigmaRF0d05_sigmaRR0d1'), 
        
end

%%%%%%%%%%%%% simulations %%%%%%%%%%%%%%
for trial= 1 %(nt-1)*5+(1:5)  % use trial =1 for testing 
rng(trial + seed_offset);

sigma_n=3.5; % input noise intensity 
% sigma_n=0;
tau_n=40;  % tau of noise 

dt=.05;  % time step (ms) 

opt.save=1; % save data 
opt.saveS1=0; % don't save spike times from Layer 1 (default is 1, save to filename, or s1_fname if specified)  
opt.saveS2=0; % don't save spike times from Layer 2 (default is 1, save to filename)
opt.CompCorr=0; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 2 & 1  
opt.loadS1=0; 
%     s1_fname=sprintf('%sRF2D3layer_GaborTheta%.04g_sigma_n%.03g_Jex%.03g_Jix%.03g_inE%.03g_inI%.03g_ID%.0f',...
%     data_folder,theta0,sigma_n,Jx(1),Jx(2),inE,0,trial);
%     s1_fname=strrep(s1_fname,'.','d');
opt.plotPopR=0; % plot population rate
opt.fixW=1; 
    Wseed1=Wseed1_range(nws);
    Wseed2=Wseed2_range(nws); 
opt.saveW=0;  % 1:save weight matrix to W_fname 
opt.savecurrent=0; 
opt.saveRm=0;
opt.Layer1only=1;  % only simulate one layer 
if trial==1
    opt.saveParam=1; 
else
    opt.saveParam=0; 
end
opt.useWfile=1; % use weight matrices from W_fname 

T=20000;  % total time of simulation (ms) 
p_stim.stim_type='OriMap_gabor_Tseg';

%%%%%%%% define firing rate of L4 neurons %%%%%%%%
rX=.01; % mean input rate 
Kc=5*2*pi; % spatial frequency 
% Nx=50; 
% theta_map=ori_map(Nx,Kc);  % generate orientation map (size Nx by Nx) 
% save([data_folder 'theta_map'],'theta_map')
load([data_folder 'theta_map4'],'theta_map')
dx=0.04;
x=-.48:dx:.48;
[X,Y]=meshgrid(x,x);
X=X(:);Y=Y(:); 
sigma=0.2;
lambda=.6;
Imag=@(theta) (exp(-(X.^2+Y.^2)/(2*sigma^2))*ones(size(theta)))...
    .*cos(2*pi/lambda*(X*cos(theta)+Y*(sin(theta))));  % Gabor image (size 625x1) 
NI=numel(Imag(0));

Filter=@(theta) exp(-(X.^2+Y.^2)/(2*sigma^2))*ones(size(theta))...
    .*cos(2*pi/lambda*(X*cos(theta)+Y*(sin(theta))))...
    /(sum(Imag(0).^2));  % Gabor filters for each L4 neuron (size(Nx^2 x NI) 
F=Filter(theta_map(:)'*pi)';  % Gabor filters for each 
fr=F*Imag(mean(testp.theta0)*pi);
F=F/mean(fr)*rX;

fr=F*Imag(testp.theta0*pi);  % firing rate for each theta  (size Nx^2 x Nth) 

p_stim.F=F; % firing rate (Nx by Nx) 
p_stim.theta_map=theta_map;
p_stim.NI=NI;
p_stim.rX=rX;
p_stim.fr=fr;
p_stim.sigma_n=sigma_n;
p_stim.sigma=sigma;
p_stim.lambda=lambda;
p_stim.T=T;
p_stim.tau_n=tau_n; 
p_stim.T_on=200;   % stim. is ON for 200 ms, then OFF for 300 ms  
p_stim.T_off=300;  
p_stim.rX_off=.005; % firing rate during OFF intervals (kHz) 
Nseg=ceil(T/(p_stim.T_on+p_stim.T_off)); % number of stim. presentations per trial 
p_stim.th_id=int8(randsample(Nth,Nseg,1));  % randomly select orientation id for each stim. presentation 
p_stim.theta0=testp.theta0;

if Nth>=3
    filename=strrep(sprintf('%sRF2D3layer_GaborTheta_sigma_n%.03g_Tuning%s_dt%.03g_ID%.0f',...
    data_folder,sigma_n,Type,dt,trial),'.','d'),
else
filename=strrep(sprintf('%sRF2D3layer_GaborTheta_Tseg_sigma_n%.03g_%s_dt%.03g_ID%.0f',...
    data_folder,sigma_n,Type,dt,trial),'.','d'),
end

ParamChange={ 'filename', filename;'dt', dt; 'T', T; 'Nc',Nc; 'p_stim',p_stim;...
    'param(1).sigmaRX',param(1).sigmaRX; 'param(1).sigmaRR',param(1).sigmaRR; 'param(1).Iapp', Iapp };
if opt.loadS1
    ParamChange=[ParamChange;{'s1_fname',s1_fname}];
end

ParamChange=[ParamChange;{'param(1).Jx', Jx}]; 

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

RF2D3layer(opt, ParamChange) % main simulation 

th_id=p_stim.th_id;
save(filename,'Jx','testp','th_id','ParamChange','-append')

end

