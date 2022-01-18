%  demo code for simulating a two-layer network with Gabor images 

mex EIF1DRFfastslowSyn.c
mex spktime2count.c

% all updated parameters need to be saved in struct 'ParamChange' 

clear 
rng('shuffle');

data_folder=''; 

%%%%%%%%% randomly select from two orientations (for computing Fisher Info) %%%%%%%%%%%%
dtheta=0.01;
testp.theta0=[.5-dtheta/2, .5+dtheta/2];  % orientations to test
Nth=length(testp.theta0); 

%%%%%% define parameters %%%%%%%%%%
sigmaRF=0.05; % ffwd width
sigmaRR=0.1; % recurrent width
Type=strrep(sprintf('L1_sigmaRF%.03g_sigmaRR%.03g',sigmaRF,sigmaRR),'.','d'),
param(1).sigmaRX=sigmaRF*ones(2,1);
param(1).sigmaRR=sigmaRR*ones(2,2);
Jx=[240; 400];
Iapp=[0; 0];
W_fname=sprintf('%sweight%s_v2',data_folder,Type),
trial=1; 

sigma_n=3.5; % input noise intensity 
tau_n=40;  % tau of noise 

dt=.05;  % time step (ms) 

opt.save=1; % save data 
opt.saveS1=1; % don't save spike times from Layer 1 (default is 1, save to filename, or s1_fname if specified)  
opt.saveS2=0; % don't save spike times from Layer 2 (default is 1, save to filename)
opt.CompCorr=0; % compute correlations 
    Nc=[500 500];  % # of E neurons sampled from Layer 2 & 1  
opt.loadS1=0; 
opt.plotPopR=0; % plot population rate
opt.fixW=0; 
opt.saveW=0; %save weight matrices to W_fname
opt.savecurrent=0; 
opt.saveRm=0;
opt.Layer1only=1;  % only simulate one layer 
if trial==1
    opt.saveParam=1; 
else
    opt.saveParam=0; 
end
opt.useWfile=0; % load weight matrices from W_fname 

T=2000;  % total time of simulation (ms) 
p_stim.stim_type='OriMap_gabor_Tseg';


%%%%%%%%%%%% define image %%%%%%%%%%
dx=0.04;
x=-.48:dx:.48;
[X,Y]=meshgrid(x,x);
X=X(:);Y=Y(:); 
sigma=0.2;
lambda=.6;
Imag=@(theta) (exp(-(X.^2+Y.^2)/(2*sigma^2))*ones(size(theta)))...
    .*cos(2*pi/lambda*(X*cos(theta)+Y*(sin(theta))));  % Gabor image (size 625x1) 
NI=numel(Imag(0));

%%%%%%%% define firing rate of L4 neurons %%%%%%%%
rX=.01; % mean input rate 
if exist('theta_map.mat','file')==0
    Nx=50;
    Kc=5*2*pi; % spatial frequency
    theta_map=ori_map(Nx,Kc);  % generate orientation map (size Nx by Nx)
    save([data_folder 'theta_map'],'theta_map')
else
    load([data_folder 'theta_map'],'theta_map')
end

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

%%%%%%%%%%%%% simulations %%%%%%%%%%%%%%
RF2D3layer(opt, ParamChange) % main simulation 

th_id=p_stim.th_id;
save(filename,'Jx','testp','th_id','ParamChange','-append')

%% %%%%%%%%%%%% plot population rate %%%%%%%%%%%%%
load(filename,'s1','T','param')
time=0:1:T;
Tw=200;
Tburn=200;
Ne1=param(1).Ne;
re1=hist(s1(1,s1(2,:)<=Ne1&s1(2,:)>0),[time T+1])/Ne1*1e3;re1=re1(1:end-1);
re1_smoothed=imfilter(re1(Tburn+1:end),ones(1,Tw)/Tw);
re1_smoothed=re1_smoothed(Tw/2-1:end-Tw/2);

figure
subplot(2,1,1)
plot(time,re1)
box off
xlim([0, T])
ylabel('FR (Hz)')
xlabel('time (ms)')
h_l=legend('E');
set(h_l,'box','off')
title('Population rate')
subplot(2,1,2)
plot(time(Tburn+Tw/2-1:end-Tw/2),re1_smoothed)
xlim([0, T])
box off
ylabel('FR (Hz)')
xlabel('time (ms)')
title('smoothed')


