% Simulation code for multi-layer networks (Fig. 9) 

% Run the following the first-time use: 
% mex EIF1DRFfastslowSyn.c
% mex spktime2count.c

% Run the following to generate orientation preference map before simulation: 
% Nx=50;
% Kc=5*2*pi; % spatial frequency
% theta_map=ori_map(Nx,Kc);
% save([data_folder 'theta_map'],'theta_map');

%  all updated parameters need to be saved in struct 'ParamChange' 

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
    

Ntrial=100; 
Np=2; % # of parameter sets 
nt=mod(job_dex-1,Ntrial)+1+Ntrial*(ceil(job_dex/Ntrial/Np)-1),  % trial number 
pid=mod(ceil(job_dex/Ntrial)-1,Np)+1,
nType=pid;

nLayer=2; % which layer to simulate 
% nLayer=2, 3, 4

%%%%%%%%  parameters  %%%%%%%%%%%%%%
sigmaRF=0.05;sigmaRR=0.1;
muI_range =[0    1.2];
Jex1_range=[140  180]; 
Jex2_range=[35   54]; 

param(1).Jx=[Jex1_range(nType); 0]; 
muI=muI_range(nType);
param(1).Iapp=[0; muI];
param(2).Iapp=[0; muI];

param(1).sigmaRX=sigmaRF*ones(2,1);
param(1).sigmaRR=sigmaRR*ones(2,2);
param(2).sigmaRX=sigmaRF*ones(2,1);
param(2).sigmaRR=sigmaRR*ones(2,2);
param(2).Prx=[.01; 0];
param(2).Jx=[Jex2_range(nType); 0];  

% param(1): parameters from Layer 1 to Layer 2 
% param(2): parameters for higher-order layers 

%%%%%%%%% two orientations (for computing Fisher Info) %%%%%%%%%%%%
dtheta=0.01;
testp.theta0=[.5-dtheta/2, .5+dtheta/2]; 
Nth=length(testp.theta0); 

%%%%%% check tuning curve %%%%%%%%%%
% testp.theta0=0.02:.02:1; 
% Nth=length(testp.theta0); 
% Ntrial=250; 

T=20000; 
for trial=nt
% rng(trial + seed_offset);

sigma_n=3.5;
tau_n=40; 
dt=.05;

opt.save=1; % save data 
opt.saveS1=0; % save s1 to s1_fname
opt.saveS2=1; % save s2 to filename
opt.CompCorr=0; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 2 & 1  
opt.plotPopR=0; % plot population rate
opt.fixW=1; 
    Wseed1=Wseed1_range(nws);
    Wseed2=Wseed2_range(nws); 
opt.useWfile=1;
    W_fname=[data_folder 'weight_L2_sigmaRF0d05_sigmaRR0d1'];
    W_fname_L1=sprintf('%sweight%s_v2',data_folder,'L1_sigmaRF0d05_sigmaRR0d1');
opt.savecurrent=0; 
opt.saveRm=0;
opt.saveW=0; 
if trial==1
    opt.saveParam=1; 
else
    opt.saveParam=0; 
end

if Nth>=3
    filename=strrep(sprintf('%sRF2D3layer_L%d_GaborTheta_th_sigma_n%.03g_muI%.03g_dt%.03g_ID%.0f',...
        data_folder,nLayer,sigma_n,muI,dt,trial),'.','d'),
else
    filename=strrep(sprintf('%sRF2D3layer_L%d_GaborTheta_Tseg_sigma_n%.03g_muI%.03g_dt%.03g_ID%.0f',...
        data_folder,nLayer,sigma_n,muI,dt,trial),'.','d'), 
end

if nLayer>2  % for higher-order layers, load input spike trains from lower layer
    opt.loadS1=1;
    if Nth>=3
        fname_prevL=strrep(sprintf('%sRF2D3layer_L%d_GaborTheta_th_sigma_n%.03g_muI%.03g_dt%.03g_ID%.0f',...
            data_folder,nLayer-1,sigma_n,muI,dt,trial),'.','d'),
    else
        fname_prevL=strrep(sprintf('%sRF2D3layer_L%d_GaborTheta_Tseg_sigma_n%.03g_muI%.03g_dt%.03g_ID%.0f',...
            data_folder,nLayer-1,sigma_n,muI,dt,trial),'.','d'),
    end
    load(fname_prevL,'s2','th_id');
    s1=s2;
    s1_fname=[fname_prevL '_s1'];
    save(s1_fname,'s1','th_id')
    clear s1 s2;
else
    opt.loadS1=0;
end

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
p_stim.T_on=200; 
p_stim.T_off=300; 
p_stim.rX_off=.005; 
Nseg=ceil(T/(p_stim.T_on+p_stim.T_off));
if opt.loadS1
    load(s1_fname,'th_id')
    p_stim.th_id=th_id;
else
    p_stim.th_id=int8(randsample(Nth,Nseg,1));
end
p_stim.theta0=testp.theta0;

F=[]; Filter=[]; theta_map=[]; X=[]; Y=[]; Imag=[];

ParamChange={'param(2).Iapp', param(2).Iapp;'param(1).Iapp', param(1).Iapp; 'filename', filename;...
    'dt', dt; 'T', T; 'Nc',Nc;'param(2).Jx',param(2).Jx;'param(1).Jx',param(1).Jx;...
    'param(1).sigmaRX',param(1).sigmaRX;'param(1).sigmaRR',param(1).sigmaRR;...
    'param(2).sigmaRX',param(2).sigmaRX;'param(2).sigmaRR',param(2).sigmaRR;...
    'param(2).Prx',param(2).Prx; 'p_stim',p_stim};
if exist('s1_fname','var')
    ParamChange=[ParamChange;{'s1_fname',s1_fname}];
end
ParamChange=[ParamChange;{'param(2).Psyn', [1, 0; 1, 0; 1, 0] }]; % only fast feedforward synapses 

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

if opt.fixW
    if opt.useWfile==1
        ParamChange=[ParamChange;{'W_fname',W_fname}];
    else
        ParamChange=[ParamChange;{'Wseed1',Wseed1; 'Wseed2',Wseed2}];
    end
end

RF2D3layer(opt, ParamChange)

th_id=p_stim.th_id;
save(filename,'muI','th_id','W_fname','-append')

if nLayer>2
    delete([s1_fname '.mat'])
end

end