function parameterContribution(fnBase,neuronsToFit)
%% Load the raw data

cut=strfind(fnBase,'_modelFit');
fn=[fnBase(1:cut-1) '.mat'];
rawData = load(fn);

nTrials = rawData.nTrials;
unitOfTime = 'ms';
binSize = 1000/rawData.sampleRate;
%% Convert data to parametric form
[rawDataParametric, parametrizationParams] = parameterizeDataForGLM(rawData);

%% Build data structure
clear expt
expt = buildGLM.initExperiment(unitOfTime, binSize, [], []);

expt = buildGLM.registerContinuous(expt, 'positionRew', 'Position along rewarded corridor', size(rawData.trial(1).position,2));
expt = buildGLM.registerContinuous(expt, 'positionUnrew', 'Position along unrewarded corridor', size(rawData.trial(1).position,2));

expt = buildGLM.registerContinuous(expt, 'inhalations', 'Animal respiration', size(rawData.trial(1).inhalations,2));

expt = buildGLM.registerContinuous(expt, 'odorRew','Rewarded odor response', size(rawData.trial(1).odor,2));
expt = buildGLM.registerContinuous(expt, 'odorUnrew', 'Unrewarded odor response', size(rawData.trial(1).odor,2));

expt = buildGLM.registerTiming(expt, 'modulationRewOdor', 'Modulation of rewarded odor response when in rewarded context');
expt = buildGLM.registerTiming(expt, 'modulationUnrewOdor', 'Modulation of unrewarded odor response when in rewarded context');

expt = buildGLM.registerContinuous(expt, 'speed', 'Animal running speed', size(rawData.trial(1).speed,2));

expt = buildGLM.registerTiming(expt, 'firstInhalOdorRew', 'First inhalation of rewarded odor');
expt = buildGLM.registerTiming(expt, 'firstInhalOdorUnrew', 'First inhalation of unrewarded odor');

expt = buildGLM.registerTiming(expt, 'licks', 'Licking');

expt = buildGLM.registerTiming(expt, 'reward', 'Animal receives reward');

expt = buildGLM.registerContinuous(expt, 'trialNumber', 'Modulation along trials', 1);

expt = buildGLM.registerTiming(expt, 'lastOdorInhalG', 'Decision phase after last odor inhalation in GO trials');
expt = buildGLM.registerTiming(expt, 'lastOdorInhalN', 'Decision phase after last odor inhalation in NO-GO trials');

expt = buildGLM.registerTiming(expt, 'firstLick', 'Modulation before first lick');


expt.trial = rawDataParametric.trial;
%%

cut=strfind(fnBase,'/');
Files=dir([fnBase(1:cut(end)) '*.mat']);
for kkk=1:length(neuronsToFit)
    numNeuron=neuronsToFit(kkk);
    
    for jj=1:length(Files)
        if contains(Files(jj).name,['fit_Neuron' sprintf('%02i', numNeuron)])
            fileNum=jj;
            continue
        end
    end
    if exist('fileNum')
        load([Files(fileNum).folder '/' Files(fileNum).name]);
        factor=5;
        if isempty(selectInd)
            bestModelLLH=[];
            save([fnBase(1:cut(end)) 'parameterContribNeuron' num2str(numNeuron) '.mat' ],'bestModelLLH');
        else
            
            kernLabs={'positionRew', 'positionUnrew','modulationRewOdor','modulationUnrewOdor','odorRew','odorUnrew','odorRewAdapt','odorUnrewAdapt','licks','inhalations','reward','speed','preGO','trialEffect','hist','population'};
            modelLabelsVar={'X','X','M','M','O','O','O','O','L','I','R','S','G','T','H','P'};
            
            
            selected_model = bestModel{selectInd};
            
            disp(['The selected model is ' [modelLabels{selected_model}]])
            
            modelUsed=[modelLabels{selected_model}];
            
            bestModelLLH=LLH(selectInd,:);
            
            llhVariable=nan(length(modelLabels),1);
            
            
            wmlFolds=wmlMat{selectInd}{maxComb(selectInd)};
            
            clear dspec dm dmTest
            dspec = buildGLM.initDesignSpec(expt);
            if ~contains(([Files(fileNum).folder '/' Files(fileNum).name]),'HistAndPop')
                dspec=dspecBuild_lab(dspec,['T' modelUsed],numNeuron);
            else
                dspec=dspecBuild_lab(dspec,['THP' modelUsed],numNeuron);
            end
            
            startIdxs=[1 (cumsum([dspec.covar(:).edim]) + 1)];
            
            
            covarLabels={dspec.covar.label};
            dt=dspec.expt.binSize/1000;
            yTest=[];
            yHatFullModel=[];
            yHatNoCovar=cell(length(modelLabels),1);
            tSamplesIn=cell(length(modelLabels),1);
            for k=1:numFolds
                sections=numFolds*factor; 
                testTrialIndices=[];
                for i=0:factor-1
                    testTrialIndices=[testTrialIndices edges(k+i*numFolds):edges(k+i*numFolds+1)-1 ];
                end
                yTest = [yTest ; buildGLM.getBinnedSpikeTrain(dspec.expt, ['sptrain' num2str(numNeuron)], testTrialIndices)];
                
                if ~contains(([Files(fileNum).folder '/' Files(fileNum).name]),'HistAndPop')
                    dmTest=buildGLM.buildMyDesingMatrix(dspec,['T' modelUsed],testTrialIndices,parametrizationParams,0,1);
                else
                    dmTest=buildGLM.buildMyDesingMatrix(dspec,['THP' modelUsed],testTrialIndices,parametrizationParams,0,1);
                end
                
                yHatFullModel= [yHatFullModel ; exp(dmTest.X*wmlFolds(k,:)')*dt];
                
                for i=1:length(modelUsed)
                    % First calculate prediction for model without variable of interest
                    if ~contains(([Files(fileNum).folder '/' Files(fileNum).name]),'HistAndPop')
                        modelToTest=setdiff(['T' modelUsed],modelUsed(i));
                    else
                        modelToTest=setdiff(['THP' modelUsed],modelUsed(i));
                    end
                    kernsIn=[];
                    for ii=1:length(modelToTest)
                        kernsIn=[kernsIn {kernLabs{find(strcmp(modelToTest(ii),modelLabelsVar))}}];
                    end
                    indsIN=[];
                    for ii=1:length(kernsIn)
                        indCovarIn=find(strcmp(covarLabels,kernsIn{ii}));
                        indsIN=[indsIN, startIdxs(indCovarIn)+1:startIdxs(indCovarIn+1)];
                        
                    end
                    indsIN=sort(indsIN);
                    indsIN=[1 indsIN]; % force in bias term
                    
                    offset=length(yHatNoCovar{selected_model(i)});
                    yHatNoCovar{selected_model(i)} = [yHatNoCovar{selected_model(i)}; exp(dmTest.X(:,indsIN) * wmlFolds(k,indsIN)')*dt];
                    
                    % Now check at what time samples that variable was active
                    
                    kernsIn={kernLabs{find(strcmp(modelUsed(i),modelLabelsVar))}};
                    
                    indsIN=[];
                    for ii=1:length(kernsIn)
                        indCovarIn=find(strcmp(covarLabels,kernsIn{ii}));
                        indsIN=[indsIN, startIdxs(indCovarIn)+1:startIdxs(indCovarIn+1)];
                        
                    end
                    indsIN=sort(indsIN);
                    
                    tSamplesIn{selected_model(i)}=[tSamplesIn{selected_model(i)}; offset+find(sum(logical(dmTest.X(:,indsIN)),2))];
                    
                    
                end
                
                
            end
            for n=1:length(modelUsed)
                
                llhVariable(selected_model(n))=sum(yTest(tSamplesIn{selected_model(n)}).*(log(yHatFullModel(tSamplesIn{selected_model(n)})) - log(yHatNoCovar{selected_model(n)}(tSamplesIn{selected_model(n)}))) - (yHatFullModel(tSamplesIn{selected_model(n)}) - yHatNoCovar{selected_model(n)}(tSamplesIn{selected_model(n)}))) / (sum(yTest(tSamplesIn{selected_model(n)}))+1);
                
                
            end
            llhVariable=llhVariable/log(2); % in bits per spike during variable action
            
            save([fnBase(1:cut(end)) 'parameterContribNeuron' num2str(numNeuron) '.mat' ],'bestModelLLH','llhVariable','modelLabels');
        end
        
        
        
        clear fileNum
    else
        bestModelLLH=[];
        save([fnBase(1:cut(end)) 'parameterContribNeuron' num2str(numNeuron) '.mat' ],'bestModelLLH');
    end
end
end

