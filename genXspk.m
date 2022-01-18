function sx=genXspk(p,Nx,T)
% p: param struct 
% sx=genXspk(stim_type,rX, Nx,cx,sigmac,Nsource,T,dt,taucorr,varargin)
switch p.stim_type
    case 'LocalCorr'
        cx=p.cx; 
        rX=p.rX;
        Nx1=sqrt(Nx);
        kc=round(Nx*cx);
        nspikes=round(T*1.1*rX/cx);
        tempspikes=sort(rand(nspikes,1)*(T));  % uniform distribution
        sx=zeros(2,kc*numel(tempspikes));
        for j=1:numel(tempspikes)
            sx(1,(j-1)*kc+1:j*kc)=tempspikes(j)+randn(kc,1)*p.taucorr;
            
            %%% Localized correlations
            spikeloc=rand(1,2);
            spikeindXs=mod(round((randn(kc,1)*p.sigmac+spikeloc(1))*Nx1)-1,Nx1);
            spikeindYs=mod(round((randn(kc,1)*p.sigmac+spikeloc(2))*Nx1)-1,Nx1);
            sx(2,(j-1)*kc+1:j*kc)=spikeindXs*Nx1+spikeindYs+1;
        end
        
    case 'Uncorr'
        % Generate uncorrelated spike trains for the feedforward layer
        sx=[];
        Nsource=p.Nsource; 
        rX=p.rX;
        for ns=1:Nsource
        tempspikes=cumsum(-log(rand(1,round(rX(ns)*(Nx/Nsource)*T*1.1)))/(rX(ns)*(Nx/Nsource)));
        tempspikes=tempspikes(tempspikes<T&tempspikes>0);
        sx_temp=zeros(2,numel(tempspikes));
        sx_temp(1,:)=tempspikes;
        sx_temp(2,:)=ceil(rand(1,size(sx_temp,2))*Nx/Nsource)+(ns-1)*Nx/Nsource;
        sx=[sx,sx_temp];
        end

    case 'GlobalCorr'
        sx=[];
        Nsource=p.Nsource; 
        rX=p.rX;cx=p.cx;
        for ns=1:Nsource
            % divid Nx for ns sources of correlation
            %%% globally correlated input
            kc_global=round(ceil(Nx/Nsource)*cx(ns));
            nspikes=round(T*1.1*rX(ns)/cx(ns));
%             tempspikes=sort(rand(nspikes,1)*(T-4*dt)+2*dt);  % uniform distribution
            tempspikes=cumsum(-log(rand(1,nspikes))/(rX(ns)/cx(ns)));
            sx_temp=zeros(2,kc_global*numel(tempspikes));
            for j=1:numel(tempspikes)
                sx_temp(1,(j-1)*kc_global+1:j*kc_global)=tempspikes(j)+randn(kc_global,1)*p.taucorr;
                
                %%% Global correlations
                spikeinds=randi(ceil(Nx/Nsource),1,kc_global)+ceil(Nx/Nsource)*(ns-1);
                sx_temp(2,(j-1)*kc_global+1:j*kc_global)=spikeinds;
            end
            sx=[sx, sx_temp];
        end
        
    case 'MultiSource' % when sigmaRx is not global 
        sx=[];
        Nsource=p.Nsource; 
        rX=p.rX;cx=p.cx;
        for ns=1:Nsource
            % divid Nx for ns sources of correlation
            %%% globally correlated input
            kc_global=round(ceil(Nx/Nsource)*cx(ns));
            nspikes=round(T*1.1*rX(ns)/cx(ns));
            tempspikes=cumsum(-log(rand(1,nspikes))/(rX(ns)/cx(ns)));
%             tempspikes=sort(rand(nspikes,1)*(T-4*dt)+2*dt);  % uniform distribution
            sx_temp=zeros(2,kc_global*numel(tempspikes));
            for j=1:numel(tempspikes)
                sx_temp(1,(j-1)*kc_global+1:j*kc_global)=tempspikes(j)+randn(kc_global,1)*p.taucorr;
                
                %%% Global correlations
                spikeinds=randi(ceil(Nx/Nsource),1,kc_global)*Nsource-(ns-1); 
                % mixed in space, mod(sx1,Nsource)=0,
                % mod(sx2,Nsource)=1, etc
                sx_temp(2,(j-1)*kc_global+1:j*kc_global)=spikeinds;
            end
            sx=[sx, sx_temp];
        end
    case 'spatialInput' % sigmac: spatial spread, centered at [.5 .5] 
         % peak rate is rX
%         fr=@(x,y) rX*exp(-((x-.5).^2+(y-.5).^2)/(2*sigmac^2))/(2*pi*sigmac^2);
%         %mean rate is rX
%         fr=@(x,y) rX*exp(-((x-.5).^2+(y-.5).^2));
%         center=[.5 .5]+[0.02, 0];
        rX=p.rX;
        center=p.center;
        Nx1=round(sqrt(Nx));
        CircRandN=@(mu,sigma,min,max,n)(mod(round(sigma*randn(n,1)+mu)-min,max-min+1)+min);
        mrate=rX*(2*pi*p.sigmac^2);
        tempspikes=cumsum(-log(rand(1,round(mrate*Nx*T*1.2)))/(mrate*Nx));
        tempspikes=tempspikes(tempspikes<T&tempspikes>0);
        sx=zeros(2,numel(tempspikes));
        sx(1,:)=tempspikes;
        Ix=CircRandN(center(1)*Nx1,(p.sigmac*Nx1),1,Nx1,numel(tempspikes));
        Iy=CircRandN(center(2)*Nx1,(p.sigmac*Nx1),1,Nx1,numel(tempspikes));
        sx(2,:)=(Ix-1)*Nx1+Iy;
    case 'OriMap'
        FR=p.FR';
        sx=cell(1,Nx); 
        for k=1:Nx
            tempspikes=cumsum(-log(rand(1,round(FR(k)*T*1.2)))/FR(k));
            tempspikes=tempspikes(tempspikes<T&tempspikes>0);
            sx{k}=[tempspikes; k*ones(1,length(tempspikes))];
        end
        sx=[sx{:}];
    case 'dynamicRate'
        dt=1; nspks=0;
        sx=zeros(2,p.rX*T*5*Nx);
        for t=dt:dt:T
            FR=p.FR(t)';
            spk=(rand(Nx,1)<FR(:)*dt);
            ns_temp=nnz(spk);
            if ns_temp
                sx(1,nspks+(1:ns_temp))=t;
                sx(2,nspks+(1:ns_temp))=find(spk>0);
                nspks=nspks+ns_temp;
            end    
        end
        sx=sx(:,1:nspks);
    case 'OriMap_2input'
        FR1=p.FR1';T1=p.T1;
        FR2=p.FR2';T2=p.T2;
        ISI=p.ISI;
        sx=cell(1,Nx);
        for k=1:Nx
            tempspikes=cumsum(-log(rand(1,round(FR1(k)*T1*1.2)))/FR1(k));
            tempspikes=tempspikes(tempspikes<T1&tempspikes>0);
            tempspikes2=cumsum(-log(rand(1,round(FR2(k)*T2*1.2)))/FR2(k));
            tempspikes=[tempspikes, tempspikes2(tempspikes2<T2&tempspikes2>0)+(T1+ISI)];
            sx{k}=[tempspikes; k*ones(1,length(tempspikes))];
        end
        sx=[sx{:}];
    case 'OriMap_gabor'
        dt=1; nspks=0;
        sx=zeros(2,p.rX*T*5*Nx);
%         noise=zeros(p.NI,1); 
        noise=p.sigma_n/sqrt(2*p.tau_n)*randn(p.NI,1); 
        for t=dt:dt:T
            noise=noise+(-noise*dt+p.sigma_n*randn(p.NI,1)*sqrt(dt))/p.tau_n;
            %     noise=sigma_n*randn(NI,1);
            %     FR=fr+F*noise;
            FR=p.fr+p.F*noise;
            spk=(rand(Nx,1)<FR*dt);
            ns_temp=nnz(spk);
            if ns_temp
                sx(1,nspks+(1:ns_temp))=t;
                sx(2,nspks+(1:ns_temp))=find(spk>0);
                nspks=nspks+ns_temp;
            end
        end
    case 'OriMap_gabor_Tseg'
        dt=1; nspks=0;
        sx=zeros(2,p.rX*T*5*Nx);
        n_seg=0; 
        for t=dt:dt:T
            
            if mod(t, p.T_on+p.T_off)<p.T_off+.1*dt
                FR=ones(size(Nx,1))*p.rX_off;
            elseif mod(t, p.T_on+p.T_off)<p.T_off+1.5*dt
                n_seg=n_seg+1;
                noise=p.sigma_n/sqrt(2*p.tau_n)*randn(p.NI,1);
                FR=p.fr(:,p.th_id(n_seg))+p.F*noise;
            else
                noise=noise+(-noise*dt+p.sigma_n*randn(p.NI,1)*sqrt(dt))/p.tau_n;
                FR=p.fr(:,p.th_id(n_seg))+p.F*noise;
            end
            spk=(rand(Nx,1)<FR*dt);
            ns_temp=nnz(spk);
            if ns_temp
                sx(1,nspks+(1:ns_temp))=t;
                sx(2,nspks+(1:ns_temp))=find(spk>0);
                nspks=nspks+ns_temp;
            end
        end
end

[~,J]=sort(sx(1,:));
sx=sx(:,J);
sx=sx(:,sx(1,:)>0&sx(1,:)<=T);

