function [rawDataParametric, parametrizationParams] = parameterizeDataForGLM(rawData)
%   Convert covariate data to a parametric format for GLM fitting
%
%   rawData : rawData obtained by running rawDataForGLM
%
%   rawDataParametric : data in parametric form ready for GLM fit
%   parametrizationParams : parameters used for paremerization
%
%% Set some parameters

threshContext(1)=366;
threshContext(2)=515;

binSizePos=16;
maxPos=657;

factorResp=5;
n_resp_bins=20;
maxResp=8*factorResp;


factorSpeed=1;
binSizeSpeed=1;

nBinsTrialEffect=14;

binSize = 1000/rawData.sampleRate;


maxOdorInhalations=12;

 
%%

rawDataParametric=rawData;

nTrials=rawData.nTrials;

posVec=1:binSizePos:maxPos;
n_pos_bins=length(posVec);
indsAssoc=find(posVec>=threshContext(1) & posVec<=threshContext(2));
posAssoc=posVec(indsAssoc);


respFreqVec=linspace(0,maxResp,n_resp_bins);

temp=[];
for i=1:nTrials
    temp=[temp; rawData.trial(i).speed ];
end
mxSpeed=prctile(temp(temp>0),99);
clear temp

maxSpeed=mxSpeed*factorSpeed;
speedVec=0:binSizeSpeed:maxSpeed;
n_speed_bins=length(speedVec);


trialEffectVec=linspace(1,nTrials,nBinsTrialEffect);
binSizeTrialEffect=trialEffectVec(2)-trialEffectVec(1);



for jTrial=1:length(rawData.trial)
    
    XTemp=zeros(size(rawData.trial(jTrial).position,1),length(posVec));
    for i=1:size(rawData.trial(jTrial).position,1)
        
        [XTemp(i,:)]= hist(find(rawData.trial(jTrial).position(i,:)),posVec);
    end
    
    
    if rawData.trial(jTrial).rewCtxt
        rawDataParametric.trial(jTrial).positionRew=XTemp;
        rawDataParametric.trial(jTrial).positionRew=conv2(ones(1,3),gausswin(3),XTemp,'same')/3;
        rawDataParametric.trial(jTrial).positionUnrew=zeros(size(XTemp));
    else
        rawDataParametric.trial(jTrial).positionUnrew=XTemp;
        rawDataParametric.trial(jTrial).positionUnrew=conv2(ones(1,3),gausswin(3),XTemp,'same')/3;
        rawDataParametric.trial(jTrial).positionRew=zeros(size(XTemp));
    end
    
    
    
    if rawData.trial(jTrial).rewCtxt & rawData.trial(jTrial).rewOdor
        rawDataParametric.trial(jTrial).association=rawData.trial(jTrial).firstOdorInhal;
    else
        rawDataParametric.trial(jTrial).association=[];
    end
    
   
    if rawData.trial(jTrial).rewCtxt & rawData.trial(jTrial).rewOdor
        
        rawDataParametric.trial(jTrial).modulationRewOdor=rawData.trial(jTrial).odor;
        rawDataParametric.trial(jTrial).modulationUnrewOdor=zeros(size(rawData.trial(jTrial).odor));
        
    elseif ~rawData.trial(jTrial).rewCtxt & rawData.trial(jTrial).rewOdor
        rawDataParametric.trial(jTrial).modulationRewOdor=zeros(size(rawData.trial(jTrial).odor));
        rawDataParametric.trial(jTrial).modulationUnrewOdor=zeros(size(rawData.trial(jTrial).odor));
        
    elseif rawData.trial(jTrial).rewCtxt & ~rawData.trial(jTrial).rewOdor
        
        rawDataParametric.trial(jTrial).modulationRewOdor=zeros(size(rawData.trial(jTrial).odor));
        rawDataParametric.trial(jTrial).modulationUnrewOdor=rawData.trial(jTrial).odor;
        
    elseif ~rawData.trial(jTrial).rewCtxt & ~rawData.trial(jTrial).rewOdor
        
        rawDataParametric.trial(jTrial).modulationRewOdor=zeros(size(rawData.trial(jTrial).odor));
        rawDataParametric.trial(jTrial).modulationUnrewOdor=zeros(size(rawData.trial(jTrial).odor));
        
    end
    
    speedTemp=zeros(size(rawData.trial(jTrial).speed,1),length(speedVec));
    for i=1:size(rawData.trial(jTrial).speed,1)
        [speedTemp(i,:)]= hist(abs(rawData.trial(jTrial).speed(i,:)*factorSpeed),speedVec);
    end
    speedTemp=conv2(ones(1,3),gausswin(3),speedTemp,'same');
    speedTemp(:,end)=0;  % SOLUTION FOR RANK DEFFICIENT
    speedTemp=speedTemp/3;
    rawDataParametric.trial(jTrial).speed=speedTemp;
    
    
    
    trialEffectTemp=zeros(size(rawData.trial(jTrial).speed,1),length(trialEffectVec));
    for i=1:size(rawData.trial(jTrial).speed,1)
        [trialEffectTemp(i,:)]=hist(jTrial,trialEffectVec);
    end
    rawDataParametric.trial(jTrial).trialEffect=conv2(ones(1,5),gausswin(5),trialEffectTemp,'same'); % SOLUTION FOR RANK DEFFICIENT
    rawDataParametric.trial(jTrial).trialEffect=rawDataParametric.trial(jTrial).trialEffect/10;
    
    numInhals=sum(rawData.trial(jTrial).odor);
    inhalSamp=find(rawData.trial(jTrial).odor);
    
    odorMatrix=zeros(size(rawData.trial(jTrial).odor,1),1); 
    odorMatrix(inhalSamp(1:min(numInhals,maxOdorInhalations)),1)= 1;
    
    odorAdapt=zeros(size(rawData.trial(jTrial).odor,1),maxOdorInhalations);
    for indInhal=2:min(numInhals,maxOdorInhalations)
        odorAdapt(inhalSamp(indInhal),indInhal)= 1;
    end
    
    if rawData.trial(jTrial).rewOdor
        
        rawDataParametric.trial(jTrial).odorRew=odorMatrix;
        rawDataParametric.trial(jTrial).odorUnrew=zeros(size(odorMatrix));
        
        rawDataParametric.trial(jTrial).odorRewAdapt=odorAdapt;
        rawDataParametric.trial(jTrial).odorUnrewAdapt=zeros(size(odorAdapt));
        
        rawDataParametric.trial(jTrial).firstInhalOdorRew=rawData.trial(jTrial).firstOdorInhal;
        rawDataParametric.trial(jTrial).firstInhalOdorUnrew=[];
        
    elseif ~rawData.trial(jTrial).rewOdor
        rawDataParametric.trial(jTrial).odorRew=zeros(size(odorMatrix));
        rawDataParametric.trial(jTrial).odorUnrew=odorMatrix;
        
        rawDataParametric.trial(jTrial).odorRewAdapt=zeros(size(odorAdapt));
        rawDataParametric.trial(jTrial).odorUnrewAdapt=odorAdapt;
        
        rawDataParametric.trial(jTrial).firstInhalOdorUnrew=rawData.trial(jTrial).firstOdorInhal;
        rawDataParametric.trial(jTrial).firstInhalOdorRew=[];
    else
        rawDataParametric.trial(jTrial).odorRew=zeros(size(odorMatrix));
        rawDataParametric.trial(jTrial).odorUnrew=zeros(size(odorMatrix));
        
        rawDataParametric.trial(jTrial).odorRewAdapt=zeros(size(odorAdapt));
        rawDataParametric.trial(jTrial).odorUnrewAdapt=zeros(size(odorAdapt));
        
        rawDataParametric.trial(jTrial).firstInhalOdorUnrew=[];
        rawDataParametric.trial(jTrial).firstInhalOdorRew=[];
    end
    
    if rawData.trial(jTrial).GO
        rawDataParametric.trial(jTrial).lastOdorInhalG=rawData.trial(jTrial).lastOdorInhal;
        rawDataParametric.trial(jTrial).lastOdorInhalN=[];
    else
        rawDataParametric.trial(jTrial).lastOdorInhalG=[];
        rawDataParametric.trial(jTrial).lastOdorInhalN=rawData.trial(jTrial).lastOdorInhal;
    end
end
parametrizationParams.binSize=binSize;
parametrizationParams.maxPos=maxPos;
parametrizationParams.n_pos_bins=n_pos_bins;
parametrizationParams.binSizePos=binSizePos;
parametrizationParams.maxSpeed=maxSpeed;
parametrizationParams.n_speed_bins=n_speed_bins;
parametrizationParams.binSizeSpeed=binSizeSpeed;
parametrizationParams.numTrials=nTrials;
parametrizationParams.nBinsTrialEffect=nBinsTrialEffect;
parametrizationParams.binSizeTrialEffect=binSizeTrialEffect;
parametrizationParams.threshContext=threshContext;
parametrizationParams.totalNeurons=sum(contains(fieldnames(rawData.trial(1)),'sptrain'));
