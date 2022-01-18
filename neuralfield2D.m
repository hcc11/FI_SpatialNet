% codes for computing eigenvalues per wavenumber of the neural field model 
% (Figs. 5A,B, 7A,B) 

clear 

data_folder='data/';
fnamesave=[data_folder 'neuralfield_stability'];

k_e=1; 
k_i=1;
phi=@(x) [k_e * x(1,:).^2.*(x(1,:)>0); k_i * x(2,:).^2.*(x(2,:)>0)];

wee0 = 80;
wei0 = -72;
wie0 = 120;
wii0 = -90;
W0 = [wee0 wei0;wie0 wii0];

taue=5; % (ms)
taui=8;

sigmae=.1; % spatial connection width of E pop 
mue=0.5; % input to E 

mui_use = 0.1:.002:.7; % input to I 
sigmai_use = .05:(.001/2):.3;  % spatial connection width of I pop 

maxreallamplot=nan(numel(sigmai_use),numel(mui_use));  % save the maximum real part of eigenvalues 
maxreallamindex=nan(numel(sigmai_use),numel(mui_use));  % the wavenumber index of the maximum real part of eigenvalues 
imag_lam = nan(numel(sigmai_use),numel(mui_use));  % save the imaginary part of eigenvalues 

Fmodes=0:1:5;
[f1,f2]=meshgrid(Fmodes,Fmodes);
wavenums=sqrt(sort(unique(f1.^2+f2.^2))); % wave numbers
Soln=zeros(2,numel(mui_use)); % save steady state firing rates 
fsolve_flag=zeros(1,numel(mui_use));

for j = 1:numel(mui_use)
    Fin = [mue; mui_use(j)];
    % find steady state for rates
    r0 = [0.1; 0.];  % initial condition  
    fun=@(r) r - phi(W0*r+Fin);
    [r,fval,exitflag,output]=fsolve(fun, r0);
    Soln(:,j)=r;
    fsolve_flag(j)=exitflag; 
    
    % Calculate eigenvalues
    for i = 1:numel(sigmai_use)
        maxreallam = zeros(1,numel(wavenums)); % eigenvalue w/ larger real part for each wave number
        imag_lam_wn = zeros(1,numel(wavenums)); 
        % Fourier coefficients of W's, wn: wave number
        weet=@(wn)(wee0*exp(-2*wn.^2*pi^2*sigmae^2));
        weit=@(wn)(wei0*exp(-2*wn.^2*pi^2*sigmai_use(i)^2));
        wiet=@(wn)(wie0*exp(-2*wn.^2*pi^2*sigmae^2));
        wiit=@(wn)(wii0*exp(-2*wn.^2*pi^2*sigmai_use(i)^2));
        Wt=@(wn)([weet(wn) weit(wn); wiet(wn) wiit(wn)]);
        
        Tau_e = taue;
        Tau_i = taui;
        
        Tau = [-1/Tau_e 0; 0 -1/Tau_i];
        
        u= W0*r+Fin;   
        ge = 2*k_e*u(1).*(u(1)>0);
        gi = 2*k_i*u(2).*(u(2)>0);
        
        GTau=[ge/Tau_e 0; 0 gi/Tau_i];
        
        for k=1:numel(wavenums)
            Matrix = Tau + GTau*Wt(wavenums(k));
            tmp = eig(Matrix);
            maxreallam(k) = max(real(tmp)); 
            imag_lam_wn(k) = abs(imag(tmp(1)));
        end
        [maxreallamplot(i,j),maxreallamindex(i,j)] = max(real(maxreallam));
        imag_lam(i,j) = imag_lam_wn(maxreallamindex(i,j));
    end
end

maxReWavNum = maxreallamindex.*(maxreallamplot>=0); % 0 for negative eigenvalues
maxReWavNum(maxReWavNum==0) = maxReWavNum(maxReWavNum==0) -1;  % convert from index to wn, -1 for negative eigenvalues
maxReWavNum(maxReWavNum>0) = wavenums(maxReWavNum(maxReWavNum>0)); % convert from index to wn

%%
figure
hold on
imagesc(mui_use,sigmai_use,maxReWavNum)
xlabel('mu_I/mu_e')
ylabel('Sigma_I/Sigma_E')

wavenums=unique(maxReWavNum); 
maxval = max(wavenums); %find maximum intensity
map = colormap; %get current colormap (usually this will be the default one)
imdata=maxReWavNum;
ind=find(imdata<-0.1);
imdata = floor((imdata./maxval)*(length(map)-1))+1; 
imrgb=ind2rgb(imdata, map);
[x, y] = ind2sub(size(imrgb),ind);
% now set colors to grey
for indx = 1:length(x)
      imrgb(x(indx),y(indx),:) = [.5 .5 .5];
end
caxis([0 maxval])
h = narrow_colorbar('vert');
Pos = get(h,'Position');
set(h,'Position',Pos+[0.01,0,0,0],'YAxisLocation','right');
set(h,'ytick',wavenums)
image(mui_use/Fin(1),sigmai_use/sigmae,imrgb) 
axis xy
dx=mui_use(2)-mui_use(1); 
dy=sigmai_use(2)-sigmai_use(1);
xlim([mui_use(1)-dx/2 mui_use(end)+dx/2]/Fin(1))
ylim([sigmai_use(1)-dy/2 sigmai_use(end)+dy/2]/sigmae)

set(gcf,'color','w')
width=7;height=6;
set(gcf,'PaperUnits','centimeters')
set(gcf,'Papersize',[width height])
set(gcf,'Paperposition',[0 0 width height])
% print('BifDiag','-dpdf','-r300','-painters')

% save(fnamesave,'res','sigmai_use','taui_use','wavenums','taue','Fin') 

%% plot imaginary part of eigenvalues  
figure
imagesc(mui_use/Fin(1),sigmai_use/sigmae,imag_lam)
hold on 
contour(mui_use/Fin(1),sigmai_use/sigmae,imag_lam,[1e-6 1e-6],'color',[1 1 1],'linewidth',2)

a=zeros(size(mui_use));
for k=1:length(a)
    tmp=find(maxReWavNum(:,k)<0,1,'last');
    if isempty(tmp)
        a(k) = 0;
    else
        a(k) = sigmai_use(tmp);
    end
end
plot(mui_use/Fin(1), a/sigmae,'color',[0.5 0.5 0.5],'linewidth',2)


%% eigenval (max real) vs wavenumber, for various sigma_i

k_e=1; 
k_i=1;
phi=@(x) [k_e * x(1,:).^2.*(x(1,:)>0); k_i * x(2,:).^2.*(x(2,:)>0)];

wee0 = 80;
wei0 = -72;
wie0 = 120;
wii0 = -90;
W0 = [wee0 wei0;wie0 wii0];

sigmai_use = .05:(.02):.2;

Fin = [0.5;  0.5]; % external inputs [mue; mui]

r0 = [0.1; 0.];

taue = 5;
taui = 8;
sigmae=.1;

fun=@(r) r - phi(W0*r+Fin);
[r,fval,exitflag,output]=fsolve(fun, r0);
r, 
exitflag, 

Fmodes=0:1:5;
[f1,f2]=meshgrid(Fmodes,Fmodes);
wavenums=sqrt(sort(unique(f1.^2+f2.^2)));
reallam=zeros(numel(sigmai_use),numel(wavenums),2);  %real parts of eigenvalues 
imaglam=zeros(numel(sigmai_use),numel(wavenums));  %imaginary parts

ss=0;
for sigmai=sigmai_use
    ss=ss+1;
    weet=@(wn)(wee0*exp(-2*wn.^2*pi^2*sigmae^2));
    weit=@(wn)(wei0*exp(-2*wn.^2*pi^2*sigmai^2));
    wiet=@(wn)(wie0*exp(-2*wn.^2*pi^2*sigmae^2));
    wiit=@(wn)(wii0*exp(-2*wn.^2*pi^2*sigmai^2));
    Wt=@(wn)([weet(wn) weit(wn); wiet(wn) wiit(wn)]);
   
    % Calculate eigenvalues
    
        Tau_e = taue;
        Tau_i = taui;
        
        Tau = [-1/Tau_e 0; 0 -1/Tau_i];
        
        u= W0*r+Fin;
        ge = 2*k_e*u(1).*(u(1)>0);
        gi = 2*k_i*u(2).*(u(2)>0);
        
        GTau=[ge/Tau_e 0; 0 gi/Tau_i];
        
        for k=1:numel(wavenums)
            Matrix = Tau + GTau*Wt(wavenums(k));
            reallam(ss,k,:)=sort((real(eig(Matrix))));
            tmp=eig(Matrix);
            imaglam(ss,k)=imag(tmp(1));
        end 
end

figure
Np=length(sigmai_use); 
colororder=copper(Np);
colororder=colororder(Np:-1:1,:); 

for ss=1:Np
    plot(wavenums,max(reallam(ss,:,:),[],3),'color',colororder(ss,:),'linewidth',1)
    hold on
end
box off 
xlim([0 7.5])
plot(xlim,[0 0],'k--','linewidth',1)
set(gcf,'color','w')
xlabel('spatial frequency')
ylabel('Max Re(\lambda)')
set(gca,'linewidth',1)
colormap(colororder)
h = narrow_colorbar('vert');
set(h,'ytick',1:Np)
set(h,'yticklabel',sigmai_use/sigmae)

width=6;height=4;
set(gcf,'PaperUnits','centimeters')
set(gcf,'Papersize',[width height])
set(gcf,'Paperposition',[0 0 width height])
% print('Lambda_sigmai','-dpdf','-r300','-painters')

%% eigenval (max real) vs wavenumber, for various mu_i
k_e=1; 
k_i=1;
phi=@(x) [k_e * x(1,:).^2.*(x(1,:)>0); k_i * x(2,:).^2.*(x(2,:)>0)];

wee0 = 80;
wei0 = -72;
wie0 = 120;
wii0 = -90;
W0 = [wee0 wei0;wie0 wii0];

taue = 5;
taui = 8;
sigmae=.1;
sigmai=.1;

mue=0.5; 
mui_use = 0.1:.1:.7;
fsolve_flag=zeros(1,numel(mui_use));

Fmodes=0:1:5;
[f1,f2]=meshgrid(Fmodes,Fmodes);
wavenums=sqrt(sort(unique(f1.^2+f2.^2)));
reallam=zeros(numel(mui_use),numel(wavenums),2);

weet=@(wn)(wee0*exp(-2*wn.^2*pi^2*sigmae^2));
weit=@(wn)(wei0*exp(-2*wn.^2*pi^2*sigmai^2));
wiet=@(wn)(wie0*exp(-2*wn.^2*pi^2*sigmae^2));
wiit=@(wn)(wii0*exp(-2*wn.^2*pi^2*sigmai^2));
Wt=@(wn)([weet(wn) weit(wn); wiet(wn) wiit(wn)]);

ss=0;
for mui=mui_use
    ss=ss+1;
    Fin=[mue; mui];
    r0 = [0.1; 0.];
    
    fun=@(r) r - phi(W0*r+Fin);
    [r,fval,exitflag,output]=fsolve(fun, r0);
    fsolve_flag(ss)=exitflag; 
    
    % Calculate eigenvalues
    
        Tau_e = taue;
        Tau_i = taui;
        
        Tau = [-1/Tau_e 0; 0 -1/Tau_i];
        
        u= W0*r+Fin;
        ge = 2*k_e*u(1).*(u(1)>0);
        gi = 2*k_i*u(2).*(u(2)>0);

        GTau=[ge/Tau_e 0; 0 gi/Tau_i];
        
        for k=1:numel(wavenums)
            Matrix = Tau + GTau*Wt(wavenums(k));
            reallam(ss,k,:)=sort((real(eig(Matrix))));
        end 
end

figure
Np=length(mui_use); 
colororder=copper(Np);
colororder=colororder(end:-1:1,:); 
for ss=1:Np
    plot(wavenums,max(reallam(ss,:,:),[],3),'color',colororder(ss,:),'linewidth',1)
    hold on
end
box off 
xlim([0 7.5])
plot(xlim,[0 0],'k--','linewidth',1)
set(gcf,'color','w')
xlabel('spatial frequency')
ylabel('Max Re(\lambda)')
set(gca,'linewidth',1)
colormap(colororder)
h = narrow_colorbar('vert');
set(h,'ytick',1:Np)
set(h,'yticklabel',mui_use/Fin(1))

width=6;height=4;
set(gcf,'PaperUnits','centimeters')
set(gcf,'Papersize',[width height])
set(gcf,'Paperposition',[0 0 width height])
% print('Lambda_mui','-dpdf','-r300','-painters')
