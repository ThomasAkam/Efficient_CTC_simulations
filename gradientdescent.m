function varargout=gradientdescent(TrX,TeX,TrY,TeY,penaliseweights,toplot,niter,W);

% This function does gradient descent to find the weights W that do the best
% linear prediction of Y as Y=WX+B.

if nargin<8;W=zeros(1,size(TrX,2));end
if nargin<7;niter=500;end
if nargin<6;toplot=0;end
if nargin<5;penaliseweights=0;end

Xdim=size(TrX,2);
ntraining=size(TrX,1);


Trainingset=[TrX,TrY];
testset=[TeX,TeY];

%do gradient decent

CovTraining=cov(Trainingset,1);     %covariance matrix of training set normalised by number of training examples
Cov_RR=CovTraining(1:Xdim,1:Xdim);  %covariance of counts over training set.
Cov_RS=CovTraining(1:Xdim,Xdim+1);  %covariance of counts with stimuli.


XmrTr=TrX-repmat(mean(TrX),ntraining,1);
YmrTr=TrY-repmat(mean(TrY),ntraining,1);

XmrTe=TeX-repmat(mean(TeX),ntraining,1);
YmrTe=TeY-repmat(mean(TeY),ntraining,1);

dt= 0.9/max(eig(Cov_RR));  %set step size

errors=inf(2,niter);
Weights=zeros(niter,Xdim);

preverrors=errors(:,1);

iter=11;
minfound=0;

while minfound<10 & iter<niter;

    dEdW=Cov_RR*W'-Cov_RS;
    dEdWsmoothing=(W-(circshift(W,[0,1])+circshift(W,[0,-1]))/2)';
    %dEdW=Cov_RR*W'-Cov_RS+penaliseweights*W';

    W=W-dt*(dEdW'+0*dEdWsmoothing');

    errors(1,iter)=mean((XmrTr*W'-YmrTr).^2);
    errors(2,iter)=mean((XmrTe*W'-YmrTe).^2);
    
    errorincrease=preverrors<errors(:,iter);
    
    minfound=minfound+errorincrease(2);

    preverrors=errors(:,iter-10);

    iter=iter+1;

    Weights(iter,:)=W;
    
end

[minerror,minpos]=min(errors(2,:));
bestweights=Weights(minpos,:);

EstYTr=TrX*bestweights';
EstYTe=TeX*bestweights';

onesYs=[ones(size(TrY)),TrY];
coefs=onesYs\EstYTr;

bestweightsscaled=bestweights/coefs(2);
Bscaled=coefs(1);

estYsscaled=TeX*bestweights'+Bscaled;


if toplot~=0;
figure(toplot)
clf()
subplot(5,1,1)
plot(errors(1,:))
hold on
plot(errors(2,:),'r')
plot(minpos,minerror,'ko')
subplot(5,1,2)
hold on
spp=max(iter-50,1);
epp=min(iter+50,niter);
plot(errors(2,spp:epp),'r')
plot(minpos-spp,minerror,'ko')
subplot(5,1,3)
%plot([bestweightsscaled(end),bestweightsscaled])
%xlim([1,size(bestweightsscaled,2)+1])
plot(bestweightsscaled)
%ylim([min(bestweightsscaled),max(bestweightsscaled)])
subplot(5,1,4)
hold on
plot([0,size(bestweights,2)+1],[0,0],'k')
hold on
plot(TeY,estYsscaled,'or')
plot([min(TeY),max(TeY)],[min(TeY),max(TeY)],'k')

subplot(5,1,5)
hist(estYsscaled-TeY,100)
xlim([-pi,pi])

end


if nargout == 1
varargout={estYsscaled};
elseif nargout==3
varargout={estYsscaled bestweightsscaled Bscaled};
else 
YTrSc=scaleests(EstYTr,TrY);
YTeSc=scaleests(EstYTe,TeY);
stimests=[YTrSc;YTeSc];
ErrorsTr=YTrSc-TrY;
ErrorsTe=YTeSc-TeY;
fractioncorrect=(sum(((TeY>0)==(YTeSc>0)))+sum(((TrY>0)==(YTrSc>0))))/(2*size(TrY,1));  
testMSE=mean(ErrorsTe.^2);
trainMSE=mean(ErrorsTr.^2);
meansqarederror=0.5*(testMSE+trainMSE)
varargout={fractioncorrect,trainMSE,testMSE,W,[],stimests};  
end

end

function stimests=scaleests(Yest,Ytrue);
ScalFac=(mean(Ytrue(Ytrue>0))-mean(Ytrue(Ytrue<0)))/(mean(Yest(Ytrue>0))-mean(Yest(Ytrue<0)));
stimests=(Yest-mean(Yest))*ScalFac+mean(Ytrue);
end
