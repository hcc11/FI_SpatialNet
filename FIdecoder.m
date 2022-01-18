function varargout=FIdecoder(datafname,fnamesave,Nid,NR,Nfile,Ntr,ds)
data=load(datafname(1));
Nstim=size(data.X,2);
ns=Nstim*Nfile;
N=length(Nid); 
X=zeros(ns,N,'int8');
th=zeros(ns,1,'int8');

for ID=1:Nfile
    data=load(datafname(ID));
    th((1:Nstim)+(ID-1)*Nstim)=data.th_id;
    X((1:Nstim)+(ID-1)*Nstim,:)=data.X(Nid,:)';
end
unique(th),

%%%%%%%%%%%%%%%%%%%%%%%%%%%
D1_idx = find(th==1);
D2_idx = find(th==2);
ns1=size(D1_idx,1),
ns2=size(D2_idx,1),
Ntr=min([Ntr,ns1,ns2]),
D1_idx=D1_idx(randsample(1:ns1,Ntr));
D2_idx=D2_idx(randsample(1:ns2,Ntr));

fracTR=1/3;
fracTE=1/3;

opt_save=1;

fm1=mean(X(D1_idx,:),1);
fm2=mean(X(D2_idx,:),1);
df=(fm2-fm1)/ds;
if NR==0 
    clear fm1 fm2; 
end

COV=(cov(double(X(D1_idx,:)))+cov(double(X(D2_idx,:))))/2;
disp('computed COV') 

T=size(D1_idx,1);
if NR==0 
    clear X D1_idx D2_idx; 
end

% save([fnamesave '_cov'],'COV','df','T','N','ds','-v7.3'); 
% disp('saved COV') 

FInaive=df*pinv(COV)*df',
FI_BC = FInaive*(2*T-2-N-1)/(2*T-2) - 2*N/(T*ds^2),
save(fnamesave,'Nid','datafname','ds','Nfile','FInaive','FI_BC','Ntr')

%%%%%%%%%% train linear decoder with early stopping %%%%%%%%%%%%%%%%%%%
if NR 
tic
maxiters=100000*10;
ns1=size(D1_idx,1);
ns2=size(D2_idx,1);
weights=zeros(N,NR);
Iters=zeros(NR,1);

% Generate Training and Testing sets
nsTR1=floor(ns1*fracTR);
nsTR2=floor(ns2*fracTR);
nsTE1=floor(ns1*fracTE);
nsTE2=floor(ns1*fracTE);
nsVAL1=ns1-nsTE1-nsTR1;
nsVAL2=ns2-nsTE1-nsTR2;

% Initialize output variables
FIVAL0 = NaN(1,NR);
FITR0 = NaN(1,NR);
FI_w = NaN(1,NR);

for k=1:NR
    idx1=randperm(ns1);
    DTR1_idx = D1_idx(idx1(1:nsTR1));
    DTE1_idx = D1_idx(idx1(nsTR1+1:nsTR1+nsTE1));
    DVAL1_idx = D1_idx(idx1(nsTR1+nsTE1+1:ns1));

    idx2=randperm(ns2);
    DTR2_idx = D2_idx(idx2(1:nsTR2));  % nsTR2 x N
    DTE2_idx = D2_idx(idx2(nsTR2+1:nsTR2+nsTE2));
    DVAL2_idx = D2_idx(idx2(nsTR2+nsTE2+1:ns2));

    %  Remove the mean of Training set from both Training and Test sets
    fmTR1=mean(X(DTR1_idx,:));
    fmTR2=mean(X(DTR2_idx,:));
    muTR = (fmTR1+fmTR2)/2;

    % Compute residuals for the training sets
    sbarTR = ds/2*(nsTR2-nsTR1)/(nsTR2+nsTR1);
    mupTR = (fmTR2*nsTR2-fmTR1*nsTR1)/(nsTR2+nsTR1)*ds/2;

    muTE = (sum(X(DTE1_idx,:)) + sum(X(DTE2_idx,:)))/(nsTE1+nsTE2) - muTR;
    mupTE = (sum(X(DTE2_idx,:)) - sum(X(DTE1_idx,:)))/(nsTE1+nsTE2)*ds/2;

    tmp = double(X(DTR1_idx,:))-ones(nsTR1,1)*muTR;
    COVTR = tmp'*tmp;
    tmp = double(X(DTR2_idx,:))-ones(nsTR2,1)*muTR;
    COVTR = (COVTR + tmp'*tmp)/(nsTR1+nsTR2);

    tmp = double(X(DTE1_idx,:))-ones(nsTE1,1)*muTR;
    M2TE = tmp'*tmp;
    tmp = double(X(DTE2_idx,:))-ones(nsTE2,1)*muTR;
    M2TE = (M2TE + tmp'*tmp)/(nsTR1+nsTR2);

    %********** Optimization loop
    dt=1/10/max(eig(COVTR)); % define stepsize
    dETEdt = -Inf;
    iters=0;
    w=zeros(N,1);  % initial vector of weights
    dETRdw = COVTR*w - mupTR';
    while(dETEdt<0 && iters < maxiters)
       iters=iters+1;
       w = w - dt*dETRdw; % update w

       dETRdw = COVTR*w - mupTR';
       dETEdw = M2TE*w - mupTE' + sbarTR*muTE';

       dETEdt = -dETRdw'*dETEdw;
    end
    if(iters==maxiters)
       fprintf('Max iters reached -- run %d\n',k),
    else
       fprintf('iter=%d -- run %d\n',iters,k),
       dETEdt,
    end
    Iters(k)=iters;
    %********** End of optimization loop

    weights(:,k)=w;

    % Estimate Fisher Information
    biasTR = (fmTR2-fmTR1)*w/ds;
    varTR = w'* (cov(double(X(DTR1_idx,:)))+cov(double(X(DTR2_idx,:)))) * w/2;
    FITR0(k) = biasTR^2/(varTR*(nsTR1+nsTR2-2)/(nsTR1+nsTR2-4)) - 2/(0.5*(nsTR1+nsTR2)*ds^2), 

    biasVAL = (mean(X(DVAL1_idx,:))-mean(X(DVAL2_idx,:)))*w/ds;
    varVAL = w'* (cov(double(X(DVAL1_idx,:)))+cov(double(X(DVAL2_idx,:)))) *w/2;
    FIVAL0(k) = biasVAL^2/(varVAL*(nsVAL1+nsVAL2-2)/(nsVAL1+nsVAL2-4)) - 2/(0.5*(nsVAL1+nsVAL2)*ds^2), %*** 2014.04 RCC added: bias correction

    if opt_save
        sprintf('saving, N=%d, nr=%d', N,k)
        save(fnamesave,'FIVAL0', 'FITR0','Iters','-append')
    end

end
end 

clear X;

if nargout==1
    varargout{1}=COV; 
end

end

