%  computes correlation as a function of preferred orientations of the pair of neurons (Fig 1D)  
data_folder='';

Nc=2000;  % # of neurons to sample each time 

sigma_n=3.5;
% sigma_n=0;

type='L1_sigmaRF0d05_sigmaRR0d1',

Ntrial=1;

datafname=@(ID) strrep(sprintf('%sGaborTheta_sigma_n%.03g_%ssum_%d',...
    data_folder,sigma_n,type,ID),'.','d'); % filename of the spike count matrix 

fname_tuning=strrep(sprintf('%sRF2D3layer_GaborTheta_sigma_n%.03g_%s_Tuning',data_folder,sigma_n,type),'.','d'),
% filename of the tuning curve function 

%%%%%%%%%%%% neurons with FR > 1 Hz %%%%%%%%%%%%%%
Tw=0.2;
FR_th=1;
data=load(datafname(1));
FR = mean(data.X,2)/Tw; 
ind_FR=find(FR>FR_th);
nnz(ind_FR),
load(fname_tuning,'E1_th','testp')
Nrep = 10; 
for rep=1:Nrep  % sample Nrep times 
    Ic=randsample(ind_FR, Nc);
    Nth=length(testp.theta0);   
    Tuning=double(E1_th(Ic,:));

    data=load(datafname(1));
    Nstim=size(data.X,2);
    ns=Nstim*Ntrial;
    
    X=zeros(Nc,ns);
    th_id=zeros(ns,1);
    
    for ID=1:Ntrial
        data=load(datafname(ID));
        th_id((1:Nstim)+(ID-1)*Nstim)=data.th_id;
        X(:,(1:Nstim)+(ID-1)*Nstim)=data.X(Ic,:);
    end
    
    COV1=cov(X(:,th_id==1)');
    COV2=cov(X(:,th_id==2)');
    
    COV=(COV1+COV2)/2;
    R=corrcov(COV);
    
    tauf=5;
    t=-tauf*5:tauf*5;
    hf=exp(-t.^2/(2*tauf^2));
    hf=hf/sum(hf);
    Tuning=imfilter(Tuning,hf,'circular');
    [alpha,Ipref]=max(Tuning,[],2);  % Ipref: preferred orientation (w/ max rate) 
    
    Cov_Ipref=zeros(Nth,Nth);
    Corr_Ipref=zeros(Nth,Nth);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [n,ind_pref]=histc(Ipref,1:Nth);
    
    for i=1:Nth
        for j=1:Nth
            if i==j
                temp=COV(ind_pref==i,ind_pref==j);
                Cov_Ipref(i,j)=mean(temp(triu(ones(size(temp)),1)==1));
                temp=R(ind_pref==i,ind_pref==j);
                Corr_Ipref(i,j)=mean(temp(triu(ones(size(temp)),1)==1));
            else
                temp=COV(ind_pref==i,ind_pref==j);
                Cov_Ipref(i,j)=mean(temp(:));
                temp=R(ind_pref==i,ind_pref==j);
                Corr_Ipref(i,j)=mean(temp(:));
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    res(rep).Corr_Ipref=Corr_Ipref;
end

%% plot 
Corr_Ipref=zeros(Nth,Nth); 
for rep=1:Nrep 
Corr_Ipref=Corr_Ipref + res(rep).Corr_Ipref; 
end 
Corr_Ipref=Corr_Ipref/Nrep; 
figure
th = 25; % index for stimulus, at 90 deg 
kshift=round(Nth/2)-th; 
Corr_Ipref=circshift(circshift(Corr_Ipref,kshift),kshift,2);  
imagesc((1:Nth)/Nth*180-90,(1:Nth)/Nth*180-90, Corr_Ipref) 
xlabel('pref. ori. -stim. ori (deg)') 
ylabel('pref. ori. -stim. ori (deg)') 

