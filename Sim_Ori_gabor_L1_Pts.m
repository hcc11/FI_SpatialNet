% simulation code for Figs. 10C,D, with long-range tuning similarity dependent
% connections. 

% Run the following the first-time use: 
% mex EIF1DRFfastslowSyn.c
% mex spktime2count.c

% Run the following to generate orientation preference map before simulation: 
% Nx=50;
% Kc=5*2*pi; % spatial frequency
% theta_map=ori_map(Nx,Kc);
% save([data_folder 'theta_map'],'theta_map');

% all updated parameters need to be saved in struct 'ParamChange' 

data_folder=''; 

Wseed1_range=[48395]; %random seed for generating connectivity matrix   
Wseed2_range=[421762]; 
nws=1; 
%%%%%%%%%%%%% this part is for running as a job array on cluster %%%%%%%%%%%
rng('shuffle');
AI = getenv('SLURM_ARRAY_TASK_ID');
job_dex = str2num(AI);
seed_offset = randi(floor(intmax/10));
rng(job_dex + seed_offset);

%%%%%%%%% two orientations (for computing Fisher Info) %%%%%%%%%%%%
Ntrial=100; % per param set
Np=2;  % number of parameter sets to test 
nt=mod(job_dex-1,Ntrial)+1+Ntrial*(ceil(job_dex/Ntrial/Np)-1), % trial number = (nt-1)*5+(1:5)
pid=mod(ceil(job_dex/Ntrial)-1,Np)+1,  % choose param set 

dtheta=0.01;
testp.theta0=[.5-dtheta/2, .5+dtheta/2];  % orientations to test
Nth=length(testp.theta0); 

%%%%%% check tuning curve %%%%%%%%%%
% testp.theta0=0.02:.02:1; 
% Nth=length(testp.theta0); 
% Ntrial=50;  % job_dex: 1:50*Np
% Np=16; 
% nt=mod(job_dex-1,Ntrial)+1+Ntrial*(ceil(job_dex/Ntrial/Np)-1), % trial number = (nt-1)*5+(1:5)
% pid=mod(ceil(job_dex/Ntrial)-1,Np)+1,


%%%%%% define parameters %%%%%%%%%%
nType=pid;

%%%%%%  panel C  %%%%%%%%
Jx=[240; 400];
Iapp=[0; 0];
sigmaRF=0.05;sigmaRR=0.1; 
param_range=[.1  .2];
P_ts=[param_range(nType); 0]; % [P_EE; P_EX], percentage of tuning-similarity dependent connections 
Type=strrep(sprintf('L1_sigmaRF%.03g_sigmaRR%.03g_Pts%.03g',sigmaRF,sigmaRR,P_ts(1)),'.','d')
param(1).sigmaRX=sigmaRF*ones(2,1);
param(1).sigmaRR=sigmaRR*ones(2,2);
W_fname=sprintf('%sweight%s_v2',data_folder,Type),

%%%%%%  panel D  %%%%%%%%
% Jx=[240; 400];
% Iapp=[0; 0];
% sigmaRF=0.05;sigmaE=0.1; sigmaI=0.3;
% param_range=[.1  .2];
% P_ts=[param_range(nType); 0]; % [P_EE; P_EX]
% Type=strrep(sprintf('L1_sigmaRF%.03g_sigmaE%.03g_sigmaI%.03g_Pts%.03g',sigmaRF,sigmaE,sigmaI,P_ts(1)),'.','d')
% param(1).sigmaRX=sigmaRF*ones(2,1);
% param(1).sigmaRR=[sigmaE sigmaI; sigmaE sigmaI];
% W_fname=sprintf('%sweight%s_v2',data_folder,Type),

%%%%%%%%%%%%% simulations %%%%%%%%%%%%%%
for trial= (nt-1)*5+(1:5)  % use trial =1 for testing 
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
opt.saveW=0;
opt.savecurrent=0; 
opt.saveRm=0;
opt.Layer1only=1;  % only simulate one layer 
if trial==1
    opt.saveParam=1; 
else
    opt.saveParam=0; 
end
opt.useWfile=1; % save weight matrices to W_fname 

T=20000;  % total time of simulation (ms) 
p_stim.stim_type='OriMap_gabor_Tseg';

%%%%%%%% firing rate of L4 neurons %%%%%%%%
rX=.01; % mean input rate 
Kc=5*2*pi; % spatial frequency 
% Nx=50; 
% theta_map=ori_map(Nx,Kc);  % generate orientation map (size Nx by Nx) 
% save([data_folder 'theta_map'],'theta_map')
load([data_folder 'theta_map'],'theta_map')
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


if opt.useWfile==1  % generate connectivity matrix  
    ParamChange=[ParamChange;{'W_fname',W_fname}];
    if exist([W_fname '.mat'], 'file')==0
        Wseed1=Wseed1_range(nws);
        param(1).Ne=200^2;
        param(1).Ni=100^2;
        param(1).Nx=50^2;
        param(1).Prr=[.01, .04; .03, .04];
        param(1).Prx=[ .1; .05];
        scurr = rng;
        if P_ts(1)>0
            th=0:.01:.99;
            Ne1=200;
            nth=numel(th);
            FR0=zeros(Ne1,Ne1,nth);
            FR1=zeros(Ne1,Ne1,nth);
            fr_th=F*Imag(th*pi);  % Nx^2 x nth
            for k=1:nth
                FR0(:,:,k)=circshift(kron(reshape(fr_th(:,k),50,50),ones(4)),[2 2]);
                FR1(:,:,k)=imgaussfilt(FR0(:,:,k),sigmaRF*Ne1,'Padding','circular');
            end
            [m1,I1]=max(FR1,[],3);
            I1=I1/nth;
            rng(Wseed1,'twister')
            [Wrr1,Wrf1]=gen_weights_tuning(param(1).Ne,param(1).Ni,param(1).Nx,param(1).sigmaRX,...
                param(1).sigmaRR,param(1).Prr, param(1).Prx,P_ts,theta_map,I1);
            save(W_fname,'Wrr1','Wrf1','P_ts','I1','Wseed1','param')
            clear FR1 FR0 m1 I1 Wrr1 Wrf1;
        else
            rng(Wseed1,'twister')
            [Wrr1,Wrf1]=gen_weights(param(1).Ne,param(1).Ni,param(1).Nx,param(1).sigmaRX,...
                param(1).sigmaRR,param(1).Prr, param(1).Prx,'2D');
            save(W_fname,'Wrr1','Wrf1','Wseed1','param')
            clear Wrr1 Wrf1;
        end
        rng(scurr);
    end
end


clear F Filter theta_map X Y Imag;

RF2D3layer(opt, ParamChange) % main simulation 

th_id=p_stim.th_id;
save(filename,'Jx','testp','th_id','ParamChange','-append')

end

