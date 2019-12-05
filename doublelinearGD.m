function [fractioncorrect,trainMSE,testMSE,P,M,stimests]=doublelinearGD(X,Y,P,M,niter);

% This function performs gradient descent to find best weight vectors P & M
% for y=P.X.M+y where X is a matrix and y and b are scalars.

%addpath('C:\MATLAB7\tprod')

nsamples=size(X,3)/2;

Pinit=P;
Minit=M;

TrX=X(:,:,1:nsamples);
TrY=Y(:,:,1:nsamples);

TeX=X(:,:,nsamples+1:end);
TeY=Y(:,:,nsamples+1:end);

Pdim=size(X,1);  %length of P vector.
Mdim=size(X,2);  %length of M vector.


if size(M,1)==0;M=zeros(Mdim,1);end

Xmr=TrX-repmat(mean(TrX,3),[1,1,nsamples]);
Ymr=TrY-repmat(mean(TrY,3),[1,1,nsamples]);

XmrTe=TeX-repmat(mean(TeX,3),[1,1,nsamples]);
YmrTe=TeY-repmat(mean(TeY,3),[1,1,nsamples]);

    
CovXX=tprod(Xmr,[1,2,-1],Xmr,[3,4,-1]);  % Covariance matrix of TrX shaped as a 4 tensor
CovXY=tprod(Ymr,[3,4,-1],Xmr,[1,2,-1]);  % Covariance of TrX with TrY.

dt=1;


errors=inf(2,niter+1);
Marray=zeros(Mdim,niter);
Parray=zeros(niter,Pdim);

dts(1)=dt;

Preg=0;  %amount of regularisation to apply to P
Mreg=0;  %amount of regularisation to apply to M

Sc=2/nsamples;

tic

errors(1,1:20)=Sc*sum((tprod(P,[1,-1],tprod(Xmr,[1,-1,3],M,[-1,2]),[-1,2,3])-Ymr).^2,3);
errors(2,1:20)=Sc*sum((tprod(P,[1,-1],tprod(XmrTe,[1,-1,3],M,[-1,2]),[-1,2,3])-YmrTe).^2,3);

prevtrainerror=errors(1,1);
minfound=0;

iter=20;

while minfound<10 & iter<niter;
    
    dEdP=Sc*(tprod(P,[1,-1],tprod(M,[-1,3],tprod(M,[-1,4],CovXX,[1,2,3,-1]),[1,-1,2]),[-1,2])-tprod(CovXY,[1,-1],M,[-1,2])');
    dEdM=Sc*(tprod(P,[2,-1],tprod(P,[3,-1],tprod(M,[-1,4],CovXX,[1,2,3,-1]),[1,2,-1]),[-1,1])-tprod(P,[1,-1],CovXY,[-1,2])');

    P=P-dt*(dEdP+Preg*2*P);
    M=M-dt*(dEdM+Mreg*2*M);
    
    errors(1,iter)=Sc*sum((tprod(P,[1,-1],tprod(Xmr,[1,-1,3],M,[-1,2]),[-1,2,3])-Ymr).^2,3)+Preg*sum(P.^2)+Mreg*sum(M.^2); %Training errors
    errors(2,iter)=Sc*sum((tprod(P,[1,-1],tprod(XmrTe,[1,-1,3],M,[-1,2]),[-1,2,3])-YmrTe).^2,3)+Preg*sum(P.^2)+Mreg*sum(M.^2); %Test errors
    
    trainerrorincrease=prevtrainerror<errors(1,iter);
    
    if trainerrorincrease==1;     
        P=P+dt*(dEdP+Preg*2*P);
        M=M+dt*(dEdM+Mreg*2*M);
        dt=dt/2;
    else
        prevtrainerror=errors(1,iter);
        prevtesterror=mean(errors(2,iter-19:iter-10));
        currtesterror=mean(errors(2,iter-9:iter));
        testerrorincrease=prevtesterror<currtesterror;
        minfound=minfound+testerrorincrease;
        Parray(iter,:)=P;
        Marray(:,iter)=M;
        iter=iter+1;
    end  
    
    
end

minfound

dt

iter

gdectime=toc

[minerror,minpos]=min(errors(2,:));

P=Parray(minpos,:);
M=Marray(:,minpos);

EstYTe=tprod(tprod(P,[1,-1],TeX,[-1,2,3]),[3,-1,1],M,[-1,2]);
EstYTr=tprod(tprod(P,[1,-1],TrX,[-1,2,3]),[3,-1,1],M,[-1,2]);

[EstYTe,ScaleFac]=scaleests(EstYTe,TeY);
[EstYTr,ScaleFac]=scaleests(EstYTr,TrY);

stimests=[EstYTr;EstYTe];

P=P*ScaleFac;

ErrorsTr=EstYTr-squeeze(TrY);
ErrorsTe=EstYTe-squeeze(TeY);

fractioncorrect=(sum((squeeze(TeY>0)==(EstYTe>0)))+sum((squeeze(TrY>0)==(EstYTr>0))))/(2*nsamples);  

testMSE=mean(ErrorsTe.^2);
trainMSE=mean(ErrorsTr.^2);
meansqarederror=0.5*(testMSE+trainMSE)

figure(2)
clf()
subplot(3,1,1)
plot(errors(1,:),'b')
hold on
plot(errors(2,:),'r')
plot(minpos,minerror,'ko')
subplot(3,1,2)
hold on
spp=max(iter-50,1);
epp=min(iter+50,niter);
plot(errors(2,spp:epp),'r')
plot(minpos-spp,minerror,'ko')
subplot(3,2,5)
hold on
plot([P(end),P])
plot([Pinit(end),Pinit],'g')
%plot(Pr,P,'o')
%plot([min(Pr),max(Pr)],[min(Pr),max(Pr)],'k')
subplot(3,2,6)
hold on
plot(M)
plot(Minit,'g')
%plot(Mr,M,'o')
%plot([min(Mr),max(Mr)],[min(Mr),max(Mr)],'k')

end

function [stimests,ScalFac]=scaleests(Yest,Ytrue);
ScalFac=(mean(Ytrue(Ytrue>0))-mean(Ytrue(Ytrue<0)))/(mean(Yest(Ytrue>0))-mean(Yest(Ytrue<0)));
stimests=(Yest-mean(Yest))*ScalFac+mean(Ytrue);
end
