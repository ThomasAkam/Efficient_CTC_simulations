function simultaneouslinearfiltering(paramfile,savefolder);
% This function takes a parameter files which defines a set of parameters for a simulation,
% runs the simulations, and saves the data generated to the specified folder, plots the results.

addpath('C:\Program Files\MATLAB\tprod')

clear global sp ip expStruct;

global sp ip expStruct;

loadparams=0;

toPlot=1;

if loadparams==0;  %specify parameters below.

    expStruct.expName='test';
    savefolder='';

    % simulation parameters

    sp.nSamplesFinal=2000;
    sp.nSamplesInit=1000;
    sp.DefaultStartSeperation=pi/4;
    sp.dt=1;
    sp.method='rreg';  %'gdec'=gradient descent.  'rreg'=ridge regression.
    sp.FFTcutfreq=300
    sp.nSpatialbins=8;

    %Input network parameters

    ip.poissonrate=10000;  % Total input rate of each sender population.
    ip.senderfreq=50;
    ip.senderfreqSD=0.2;
    ip.senderamp=2;
    ip.senderampSD=0;
    ip.distractorfreqs=[50,50,50,50];%[12,20,110,145];%
    ip.distrangfreqSD=0.2;  %determines how broadly the frequency is smeared out
    ip.distractoramp=2;
    ip.distampSD=0;
    ip.ndistractors=4;
    ip.windowlength=80/sp.dt;
    ip.detailedbalance=0;
    ip.syncp=3; 
    ip.oscshape='VonMises'%'Sin';  %
    ip.coherentdistractors=0;
    ip.knownmodulation=1;
    ip.RelativeCutFreq=0.5;

    paramsouter={'ip.senderamp'}
    paramsoutervalues=[2];

    paramsinner={'ip.poissonrate'}
    paramsinnervalues=logspace(log10(5000),log10(100000),2)/sp.nSpatialbins;%

else
    
    %load('data\AsyncDistractors\ADeffectofAmplitude.mat')   %Or use previously saved parameters.
    load(paramfile)
    
    expStruct.expName=ParamStruct.Name;
    sp=ParamStruct.simparams;
    ip=ParamStruct.inputparams;
    paramsouter=ParamStruct.paramsouter;
    paramsoutervalues=ParamStruct.paramsoutervalues;
    paramsinner=ParamStruct.paramsinner;
    paramsinnervalues=ParamStruct.paramsinnervalues;
    
end
    
sp.toPlot=toPlot;

expStruct.params.simParams=sp;expStruct.params.inputParams=ip;
expStruct.params.paramsouter=paramsouter;expStruct.params.paramsinner=paramsinner;
expStruct.params.outervalues=paramsoutervalues;expStruct.params.innervalues=paramsinnervalues;

ParamSequencer(paramsinner,paramsinnervalues,paramsouter,paramsoutervalues,savefolder);

end

function ParamSequencer(paramsinner,paramsinnervalues,paramsouter,paramsoutervalues,savefolder);
% This function Iterates through the different paramters, calling the
% simulator and plotting the results.

global sp ip ppos expStruct;

for outerpos=1:size(paramsoutervalues,2);
    ppos.outer=outerpos;
    plotcolor=rand(1,3);

    for pout=1:size(paramsouter,2); % Iterate through outer paramter.
        eval([paramsouter{pout},'=',num2str(paramsoutervalues(pout,ppos.outer))])
    end
    
    sp.startseperation=sp.DefaultStartSeperation;
    Meansquarederrors=[];
    FinalFractionsCorrect=[];

    for innerpos=1:size(paramsinnervalues,2);  %Iterate through inner paramter.
        ppos.inner=innerpos;
        for pin=1:size(paramsinner,2);
            eval([paramsinner{pin},'=',num2str(paramsinnervalues(pin,ppos.inner))])
        end
        
        [meansqarederror,endseperation,FinalFractionCorrect]=FindCorrectSeperation;
        Meansquarederrors=[Meansquarederrors,meansqarederror]
        FinalFractionsCorrect=[FinalFractionsCorrect,FinalFractionCorrect]
        sp.startseperation=endseperation;
        
        save(strcat(savefolder,expStruct.expName,'Data.mat'),'expStruct')
        

        figure(95) %Plot performance as fn of paramters.
        subplot(2,1,1)
        hold on
        plot(paramsinnervalues(1,1:size(Meansquarederrors,2)),sqrt(Meansquarederrors),'color',plotcolor);
        subplot(2,1,2)
        hold on
        plot(paramsinnervalues(1,1:size(Meansquarederrors,2)),1./Meansquarederrors,'color',plotcolor);  

    end
    
    figure(93)  %Plot performance as fn of paramters. 
    subplot(2,1,1)
    hold on
    plot(paramsinnervalues(1,1:size(Meansquarederrors,2)),sqrt(Meansquarederrors),'color',plotcolor);
    subplot(2,1,2)
    hold on
    plot(paramsinnervalues(1,1:size(Meansquarederrors,2)),1./Meansquarederrors,'color',plotcolor);
end
end

function [meansqarederror,seperation,finalfractioncorrect]=FindCorrectSeperation;
%  This function finds the appropriate seperation to give 75-80% correct
%  discrimination.

global sp;

minSep=0;
maxSep=pi;

seperation=sp.startseperation;  %seperation of the two stimuli used to train the LOLE.

[fractioncorrect,meansqarederror]=sequentialweightsdecoder(seperation,sp.nSamplesInit);
fractionscorrect=fractioncorrect;seperations=seperation;

seperationok=0;

while seperationok==0;
    
    if fractioncorrect>=0.75 & fractioncorrect<=0.8;
        disp('Correct seperation')
        [finalfractioncorrect,meansqarederror]=sequentialweightsdecoder(seperation,sp.nSamplesFinal);
        break
    elseif fractioncorrect<0.75;
        minSep=seperation;
        seperation=0.5*(seperation+maxSep);
    if seperation>3.05;
        disp('Widest separation insuficient')
        [finalfractioncorrect,meansqarederror]=sequentialweightsdecoder(seperation,sp.nSamplesFinal);
        break
    end
    elseif fractioncorrect>0.8;
        maxSep=seperation;
        seperation=0.5*(seperation+minSep);
    end
    
    [fractioncorrect,meansqarederror]=sequentialweightsdecoder(seperation,sp.nSamplesInit);
    fractionscorrect=[fractionscorrect,fractioncorrect]
    seperations=[seperations,seperation]

    figure(234)    % plot fraction correct against separation.
    subplot(1,2,2);
    cla()
    sortedSepFrac=sortrows([seperations',fractionscorrect'],1);
    semilogx(sortedSepFrac(:,1),sortedSepFrac(:,2),'o-')

    
end
end


function [fractioncorrect,meansqarederror]=sequentialweightsdecoder(seperation,nsamples);

% This version iteritavely fits wieghts both in time and for the LOLE.

global sp ip ppos expStruct;

[allspikearrays,alltargetmodulations,stimangles]=GenerateInputActivity(seperation,nsamples);

modffts=fft(repmat(hanning(ip.windowlength)',nsamples,1).*alltargetmodulations,[],2);

spikeffts=fft(allspikearrays,[],2);

HzperFFTpoint=1000/size(modffts,2);
cutfreq=sp.FFTcutfreq;
Nfftpoints=ceil(cutfreq/HzperFFTpoint);

absmodffts=abs(modffts(:,1:Nfftpoints));   
absspikeffts=abs(spikeffts(:,1:Nfftpoints,:));
spikemodphasedifs=angle(spikeffts(:,1:Nfftpoints,:))-repmat(permute(angle(modffts(:,1:Nfftpoints)),[3,2,1]),[sp.nSpatialbins,1,1]);

fouriercoefs=[absspikeffts.*repmat(permute(absmodffts,[3,2,1]),[sp.nSpatialbins,1,1]).*cos(spikemodphasedifs),absspikeffts.*repmat(permute(absmodffts,[3,2,1]),[sp.nSpatialbins,1,1]).*sin(spikemodphasedifs)];
fouriercoefs(:,1+size(fouriercoefs,2)/2,:)=[];
fouriercoefs=fouriercoefs/mean(mean(mean(fouriercoefs)));

clear modffts spikeffts absmodffts absspikeffts spikemodphasedifs

disp('Fitting initial temporal weights')

Weightst=initialtemporalweightgd(fouriercoefs,stimangles(:,1),alltargetmodulations,sp.sbinangles);
Weightss=zeros(1,sp.nSpatialbins);

disp('Fitting all weights')
[fractioncorrect,trainMSE,testMSE,Weightss,Weightst,stimests]=doublelinearGD(fouriercoefs,permute(stimangles(:,1),[3,2,1]),Weightss,Weightst',10000);

meansqarederror=0.5*(trainMSE+testMSE)

if nsamples==sp.nSamplesFinal
    
 expStruct.data.stimAngles{ppos.outer,ppos.inner}=stimangles;
 expStruct.data.stimests{ppos.outer,ppos.inner}=stimests;
 expStruct.data.seperation(ppos.outer,ppos.inner)=seperation;
 expStruct.data.MeanSquaredErrorTrain(ppos.outer,ppos.inner)=trainMSE;
 expStruct.data.MeanSquaredErrorTest(ppos.outer,ppos.inner)=testMSE;
 expStruct.data.fractioncorrect(ppos.outer,ppos.inner)=fractioncorrect;
 expStruct.data.Spatialweights{ppos.outer,ppos.inner}=Weightss;
 expStruct.data.temporalweights{ppos.outer,ppos.inner}=Weightst;
 

end
end

function [allspikearrays,alltargetmodulations,stimangles]=GenerateInputActivity(seperation,nsamples);
% This function generates the input activity.

global sp ip;

sp.sbinangles=[1:sp.nSpatialbins]'*2*pi/sp.nSpatialbins;

targetstimangles=reshape([ones(1,nsamples/2);-ones(1,nsamples/2)],nsamples,1)*seperation/2;
stimangles=[targetstimangles,rand(nsamples,ip.ndistractors)*2*pi]; %stimulus values

modulations=zeros(ip.ndistractors+1,ip.windowlength*nsamples);

[modulations(1,:),senderphase]=oscillationdistorter(ip.senderfreq,ip.senderamp,ip.senderampSD,ip.senderfreqSD,ip.windowlength*nsamples,sp.dt,ip.syncp,ip.oscshape,ip.RelativeCutFreq);

dphases=[1:ip.ndistractors]*2*pi/(ip.ndistractors+1);

if size(ip.distractorfreqs,2)<ip.ndistractors;
    ip.distractorfreqs=repmat(ip.distractorfreqs(1),1,ip.ndistractors);
end    
    
for d=1:ip.ndistractors;
    if ip.coherentdistractors==0;
    modulations(d+1,:)=oscillationdistorter(ip.distractorfreqs(d),ip.distractoramp,ip.distampSD,ip.distrangfreqSD,ip.windowlength*nsamples,sp.dt,ip.syncp,ip.oscshape,ip.RelativeCutFreq);
    else  
    oscampvec=ip.distractoramp*(1+ip.distampSD*filteredwhitenoise(ip.windowlength*nsamples,sp.dt,ip.RelativeCutFreq*ip.senderfreq));
    if ip.distractorphase==-1;
        modulations(d+1,:)=exp(oscampvec.*cos(senderphase+dphases(d)))./besseli(0,oscampvec);  
    else
        modulations(d+1,:)=exp(oscampvec.*cos(senderphase+ip.distractorphase))./besseli(0,oscampvec);     
    end
    end
end

allspikearrays=zeros(sp.nSpatialbins,ip.windowlength,nsamples);
alltargetmodulations=zeros(nsamples,ip.windowlength);


for s=1:nsamples;
    if rem(s,5000)==0;
        s
        
        
    end
    samplestimangles=stimangles(s,:);
    spatialrates=srates(sp.sbinangles,samplestimangles);
    rate=ip.poissonrate*(spatialrates*modulations(:,((s-1)*ip.windowlength)+1:s*ip.windowlength));
    if ip.detailedbalance==1;
        ispatialrates=sum(spatialrates(:,2:end),2);
        irate=ip.poissonrate*repmat(ispatialrates,1,ip.windowlength);
        spikes=poissrnd(rate*sp.dt/1000)-poissrnd(irate*sp.dt/1000);
    else
        spikes=poissrnd(rate*sp.dt/1000);
    end
    if ip.knownmodulation==1;
    alltargetmodulations(s,:)=modulations(1,((s-1)*ip.windowlength)+1:s*ip.windowlength);
    else
    alltargetmodulations(s,:)=sum(spikes,1);
    end
    allspikearrays(:,:,s)=spikes;
end

if sp.toPlot==1;
figure(233)    
 clf()
 subplot(1,2,1)
 plot(modulations(1,1:ip.windowlength),'r')
 hold on
 for d=1:ip.ndistractors; 
 plot((2.5*d)+modulations(d+1,1:ip.windowlength),'k')
 end
 subplot(2,2,2)
plot([1:8]',spatialrates);


[powerspec,freqvec] = pwelch(modulations(1,:),1000/sp.dt,[],2^nextpow2(1000/sp.dt),1000/sp.dt);
subplot(2,2,4)
hold on
plot(freqvec,sqrt(powerspec),'r')
xlim([0,400])

[powerspec,freqvec] = pwelch(modulations(2,:),1000/sp.dt,[],2^nextpow2(1000/sp.dt),1000/sp.dt);
ylim([0,0.3])
plot(freqvec,sqrt(powerspec),'k')
xlim([0,400])
end
end


function varargout=oscillationdistorter(centerfreq,amplitude,ampSD,angfreqSD,duration,dt,syncp,oscshape,RelCutFreq);

%This function generates irregular oscillation waveforms.

phase=rand(1)*2*pi;
meanangfreq=2*pi*centerfreq/1000;
cutofffreq=centerfreq*RelCutFreq;
angfreqvec=meanangfreq*(1+angfreqSD*filteredwhitenoise(duration,dt,cutofffreq));
oscampvec=amplitude*(1+ampSD*filteredwhitenoise(duration,dt,cutofffreq));
osctime=cumsum(angfreqvec*dt,2)+repmat(phase',1,size(angfreqvec,2));  %this is operational time in which each period of oscillation has same duration.
oscphase=rem(osctime,2*pi);

if strcmp('Sin',oscshape);
allsyncparams=[[1;1],[2/3;2],[8/35;4],[16/231;6],[128/6435;8]];
syncparams=allsyncparams(:,syncp);
modulation=(1-amplitude)+oscampvec.*(syncparams(1)*(1+sin(osctime)).^syncparams(2));
elseif strcmp('VonMises',oscshape);
modulation=exp(oscampvec.*cos(osctime))./besseli(0,oscampvec);
else
    error('Unrecognised oscillation waveform name ip.oscshape')
end

if min(modulation)<0;
    error('Minimum value of firing rate modulation less than zero, check your parameters!')
end

if nargout == 1
varargout={modulation};
else
varargout={modulation,oscphase};
end

end

function Weightst=initialtemporalweightgd(fouriercoefs,targetstimangles,alltargetmodulations,sbinangles);

nsamples=size(targetstimangles,1);
nspatialbins=size(sbinangles,1);
windowlength=size(alltargetmodulations,2);

traincoefs=fouriercoefs(:,:,1:nsamples/2);
testcoefs=fouriercoefs(:,:,1+nsamples/2:end);

trainstimangles=targetstimangles(1:nsamples/2,:);
teststimangles=targetstimangles(1+nsamples/2:end,:);

%targetmodulations=reshape(repmat(permute(alltargetmodulations,[3,1,2]),nspatialbins,1),[],windowlength,1);

trainfouriercoefs=reshape(permute(traincoefs,[2,1,3]),size(traincoefs,2),[])';
testfouriercoefs=reshape(permute(testcoefs,[2,1,3]),size(testcoefs,2),[])';

trainfouriercoefs=trainfouriercoefs/mean(mean(trainfouriercoefs));
testfouriercoefs=testfouriercoefs/mean(mean(testfouriercoefs));

traintargetrates=reshape(2/3*(1+cos(repmat(sbinangles,1,nsamples/2)+repmat(trainstimangles',size(sbinangles)))).^2,[],1);
testtargetrates=reshape(2/3*(1+cos(repmat(sbinangles,1,nsamples/2)+repmat(teststimangles',size(sbinangles)))).^2,[],1);


[ests,Weightst,B]=gradientdescent(trainfouriercoefs,testfouriercoefs,traintargetrates,testtargetrates,0,12,6000);
end

function spatialrates=srates(sbinangles,samplestimangles);

spatialrates=2/3*(1+cos(repmat(sbinangles,1,size(samplestimangles,2))+repmat(samplestimangles,size(sbinangles)))).^2;
%spatialrates=8/35*(1+cos(repmat(sbinangles,1,size(samplestimangles,2))+repmat(samplestimangles,size(sbinangles)))).^4;

end

function [stimests,ScalFac]=scaleests(Yest,Ytrue);
ScalFac=(mean(Ytrue(Ytrue>0))-mean(Ytrue(Ytrue<0)))/(mean(Yest(Ytrue>0))-mean(Yest(Ytrue<0)));
stimests=(Yest-mean(Yest))*ScalFac+mean(Ytrue);
end


function signal=filteredwhitenoise(duration,dt,cutofffreq);
%  This function generates normally distributed zero mean filtered white
%  noise with a mean power (variance) of 1 and a flat spectrum up to cutofffreq and zero power above that.

sectionlength=duration/dt;
Hzperfftpoint=1000/duration;
Nnonzeropoints=floor(cutofffreq/Hzperfftpoint);

fftbit=0.5*randn(1,Nnonzeropoints)+i*0.5*randn(1,Nnonzeropoints);
ftranform=(sqrt(1000*duration/cutofffreq)/dt)*[randn(1),fftbit,zeros(1,sectionlength-1-2*Nnonzeropoints),conj(fliplr(fftbit))];
signal=ifft(ftranform);

end


