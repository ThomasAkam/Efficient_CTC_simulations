function MakeParameterFile();

% Make parameter file to assess affect of integration time using phase division distractors.

Name='PhaseSepEfofIntTime';

% simulation parameters

sp.nSamplesFinal=10000;
sp.nSamplesInit=2000;
sp.DefaultStartSeperation=pi;
sp.dt=1;
sp.method='gdec';  %'gdec'=gradient descent.  'rreg'=ridge regression.
sp.FFTcutfreq=500;
sp.nSpatialbins=8;

%Input network parameters

ip.poissonrate=50000/sp.nSpatialbins;
ip.senderfreq=50;
ip.senderfreqSD=0.1;
ip.senderamp=5.3;
ip.senderampSD=0.1;

ip.distractorfreqs=50;
ip.distrangfreqSD=0.1;  %determines how broadly the frequency is smeared out
ip.distractoramp=5.3;
ip.distampSD=0.1;
ip.distractorphase=-1;

ip.ndistractors=3;

ip.windowlength=100/sp.dt;
ip.detailedbalance=0;
ip.syncp=1; 
ip.oscshape='VonMises';%'Sin';  %
ip.coherentdistractors=1;
ip.knownmodulation=1;
ip.RelativeCutFreq=0.5;

sp
ip

paramsouter=[{'ip.senderfreq'}]
paramsoutervalues=[12.5,25,50,100]

paramsinner=[{'ip.windowlength'}]
paramsinnervalues=[1,2,4,6,10,16,25,40,63,100,158,251,398,631,1000];
  

ParamStruct.Name=Name;
ParamStruct.simparams=sp;
ParamStruct.inputparams=ip;
ParamStruct.paramsouter=paramsouter;
ParamStruct.paramsoutervalues=paramsoutervalues;
ParamStruct.paramsinner=paramsinner;
ParamStruct.paramsinnervalues=paramsinnervalues;


save(strcat(ParamStruct.Name,'.mat'),'ParamStruct')
end