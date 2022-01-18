%% compute FFT over space and time 
% run Sim_Spont.m first 

rng('shuffle');

Types={'L1_sigmaRF0d05_sigmaE0d1_sigmaI0d1','L1_sigmaRF0d05_sigmaE0d1_sigmaI0d2',...
    'L1_sigmaRF0d05_sigmaE0d1_sigmaI0d3',...
    'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d4',...
    'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI0d8','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI1d2',...
    'L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI1d6','L1_sigmaRF0d05_sigmaRR0d1_Jix0_muI2'};

data_folder='data/';

Np=length(Types);
Ntrial=5;

Ne=4e4;
T=2e4;
Tw=200;
dt=1;
Tburn=1000;
Tshift=100;
Nt=floor((T-Tburn-Tw)/Tshift)+1;

fnamesave=[data_folder 'SpatioTempFreq_Spont_sigI_muI'],

for pid=1:Np
    Type=Types{pid};
    
    Y=zeros(200,200,floor(Tw/dt));
    rate=zeros(1,Ne);
    for trial=1:Ntrial
        filename=strrep(sprintf('%sRF2D3layer_Spont_%s_dt0d05_ID%.0f',...
            data_folder,Type,trial),'.','d'),
        load(filename,'s1')
        rate=rate+hist(s1(2,s1(1,:)<=T&s1(1,:)>Tburn),1:Ne)/(T-Tburn);
    end
    rate=rate/Ntrial;
    
    for trial=1:Ntrial
        filename=strrep(sprintf('%sRF2D3layer_Spont_%s_dt0d05_ID%.0f',...
            data_folder,Type,trial),'.','d'),
        load(filename,'s1')
        
        for k=1:Nt
            s1_seg=s1(:,s1(1,:)<(Tburn+(k-1)*Tshift+Tw)&s1(1,:)>=(Tburn+(k-1)*Tshift));
            s1_seg(1,:)=s1_seg(1,:)-(Tburn+(k-1)*Tshift);
            E1=(spktime2count(s1_seg,1:Ne,dt,floor(Tw/dt),1))-repmat(rate',[1 Tw]);  % subtract mean
            Y = Y+abs(fftshift(fftn(reshape(E1,Ne1,Ne1,[])))).^2;
        end
    end
    Y=Y/Nt/Ntrial;
    res(pid).Y=Y;
    res(pid).rate=rate;
    
end

save(fnamesave,'res','Ntrial','Types','Tw','dt','Tburn','Tshift','Nt','T')

%% plot 
data_folder='data/';
load([data_folder 'SpatioTempFreq_Spont_sigI_muI'])

Np=length(Types);
Fs=1e3/dt;
f_range=Fs*(0:(Tw/2-1))/Tw;

for pid=1:Np
    Y=res(pid).Y;
    Y=Y/Ne^2/(Tw)*1e3*dt^2;
    L = size(Y,3);
    
    Y0=Y(:,:,round(L/2)+1);
    Y=reshape(Y,[],L);
    
    xn=-100:99;
    [Xn,Yn]=meshgrid(xn,xn);
    Wavn=(Xn.^2+Yn.^2);  % wave number 
    wn=unique(Wavn(:));
    wn=wn(wn<=100);
    [n, edges, bin]=histcounts(Wavn,[wn-.5; wn(end)+.5]);
    
    res(pid).SF=zeros(length(wn),L/2);
    res(pid).F0=zeros(length(wn),1);
    for i=1:length(wn)
        res(pid).SF(i,:)=mean(Y(bin(:)==i,(L/2+1):end),1);
        res(pid).F0(i)=mean(Y0(bin(:)==i));
    end
    res(pid).Y0 = Y0; % power at omega=0;

end

figure % sigma_i
colororder1=copper(3+1);
colororder1=colororder1((3:-1:1)+1,:);
for k=1:3
    subplot(1,4,k)
    pid=k;
    imagesc(xn,xn,res(pid).Y0)
    axis([-10 10 -10 10])
    subplot(1,4,4)
    hold on
    plot(sqrt(wn),res(pid).F0,'color',colororder1(k,:))
    
end

figure  % mu_i 
colororder1=copper(5+1);
colororder1=colororder1((5:-1:1)+1,:);
for k=1:5
    subplot(1,6,k)
    pid=k+3;
    imagesc(xn,xn,res(pid).Y0)
    axis([-10 10 -10 10])
    subplot(1,6,6)
    hold on
    plot(sqrt(wn),res(pid).F0,'color',colororder1(k,:))
    
end
