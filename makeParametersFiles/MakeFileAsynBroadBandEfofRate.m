function MakeParameterFile();

% Make parameter file to assess affect of input networks rate using broad band asynchronous disctractors.

Name='AsynBroadBandEfofRate';

% simulation parameters

sp.nSamplesFinal=10000;
sp.nSamplesInit=2000;
sp.DefaultStartSeperation=pi/4;
sp.dt=1;
sp.method='gdec'; 
sp.FFTcutfreq=500;
sp.nSpatialbins=8;

%Input network parameters

ip.poissonrate=10000;
ip.senderfreq=50;
ip.senderfreqSD=0.3;
ip.senderamp=1.16;
ip.senderampSD=0.1;

ip.distractorfreqs=10;
ip.distrangfreqSD=0.1;  %determines how broadly the frequency is smeared out
ip.distractoramp=0;
ip.distampSD=0;

ip.ndistractors=3;

ip.windowlength=100/sp.dt;
ip.detailedbalance=0;
ip.syncp=1; 
ip.oscshape='VonMises';%'Sin';  %
ip.coherentdistractors=0;
ip.knownmodulation=1;
ip.RelativeCutFreq=0.5;

sp
ip

paramsouter=[{'ip.distractoramp'}]
paramsoutervalues=[0]

paramsinner=[{'ip.poissonrate'}]
paramsinnervalues=logspace(log10(1),log10(20),10)*10000/sp.nSpatialbins
  

ParamStruct.Name=Name;
ParamStruct.simparams=sp;
ParamStruct.inputparams=ip;
ParamStruct.paramsouter=paramsouter;
ParamStruct.paramsoutervalues=paramsoutervalues;
ParamStruct.paramsinner=paramsinner;
ParamStruct.paramsinnervalues=paramsinnervalues;


save(strcat(ParamStruct.Name,'.mat'),'ParamStruct')
end