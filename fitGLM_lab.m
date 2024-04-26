function fitGLM_lab(fnRaw,neuronsToFit,fitModel)

%   Fit GLM models and select best model
%   
%   fnRaw : ouptut file obtained by running rawDataForGLM
%   neuronsToFit :  list of neuron numbers to fit. 
%   fitModel :  string specifying specific model to fit, e.g., 'XAOL', according to modelLabels.
%               Otherwise 'all' will fit all possible models
%
%% Set model labels

modelLabels{1}='X'; % Position in virtual corridor
modelLabels{2}='O'; % Odor response
modelLabels{3}='L'; % Licking
modelLabels{4}='I'; % Inhalation
modelLabels{5}='R'; % Reward
modelLabels{6}='S'; % Running speed
modelLabels{7}='G'; % Pre-GO
modelLabels{8}='M'; % Contextual modulation of odor response

%% Set cross-validation parameters
numFolds=10; factor=5; % 10-fold cross-validation

%% Load the raw data
rawData = load(fnRaw); 

numTrials = rawData.nTrials; % number of trials
unitOfTime = 'ms';
binSize = 1000/rawData.sampleRate;

sections=numFolds*factor; edges=(round(linspace(1,numTrials+1,sections+1)));

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

expt = buildGLM.registerTiming(expt, 'modulationRewOdor', 'Contextual modulation of rewarded odor response when in rewarded context');
expt = buildGLM.registerTiming(expt, 'modulationUnrewOdor', 'Contextual modulation of unrewarded odor response when in rewarded context');

expt = buildGLM.registerContinuous(expt, 'speed', 'Animal running speed', size(rawData.trial(1).speed,2));

expt = buildGLM.registerTiming(expt, 'firstInhalOdorRew', 'First inhalation of rewarded odor');
expt = buildGLM.registerTiming(expt, 'firstInhalOdorUnrew', 'First inhalation of unrewarded odor');

expt = buildGLM.registerTiming(expt, 'licks', 'Licking');

expt = buildGLM.registerTiming(expt, 'reward', 'Animal receives reward');

expt = buildGLM.registerContinuous(expt, 'trialNumber', 'Modulation along trials', 1); 

expt = buildGLM.registerTiming(expt, 'lastOdorInhalG', 'Decision phase after last odor inhalation in GO trials');
expt = buildGLM.registerTiming(expt, 'lastOdorInhalN', 'Decision phase after last odor inhalation in NO-GO trials');

expt = buildGLM.registerTiming(expt, 'firstLick', 'First lick response'); 


expt.trial = rawDataParametric.trial;

cut=strfind(fnRaw,'.mat');
if ~exist([fnRaw(1:cut-1) '_exptData.mat'],'file')
    save([fnRaw(1:cut-1) '_exptData.mat'],'expt')
end
%% RUN GLM FIT AND TEST PROCEDURE

% Define all possible model combinations
totalVars=length(modelLabels);
clear modelType
for numVars=1:totalVars
    modelType{numVars}=nchoosek(1:totalVars,numVars);
end
numModels=0;
for i=1:length(modelType)
    numModels=numModels+size(modelType{i},1);
end

binfun = expt.binfun; 

if strcmp(fitModel,'all') % FIT ALL POSSIBLE MODELS
    %%
    for kkk=1:length(neuronsToFit) % RUN THROUGH NEURONS TO FIT
        numNeuron=neuronsToFit(kkk);
               
        wmlMat=cell(1,totalVars);
        wmlMean=cell(1,totalVars);
        hessian=cell(1,totalVars);
        testFit=cell(1,totalVars);
        trainFit=cell(1,totalVars);
        rho=cell(1,totalVars);
        
        LLH=nan(totalVars,numFolds);
        AIKs=nan(totalVars,numFolds);
        bestModel=cell(totalVars,1);
        maxComb=nan(totalVars,1);
        p_llh=nan(totalVars,1);
        
        tic;
        for numVars=1:totalVars
            count=1;
            totalCombs=size(modelType{numVars},1);
            wmlMat{numVars}=cell(totalCombs,1);
            
            if numVars>1
                included=logical(sum(ismember(modelType{numVars},modelType{numVars-1}(maxComb(numVars-1),:)),2)==numVars-1);
            else
                included=1:totalVars;
            end
            
            %% 
            % FOR NEURON numNeuron, RUN THROUGH ALL MODEL COMBINATIONS WITH A TOTAL NUMBER OF VARIABLES numVars
            for varComb=1:totalCombs 
                
                varsToInclude=modelType{numVars}(varComb,:);
                if ~included(varComb)
                    clear dspec dm dmTest
                    dspec = buildGLM.initDesignSpec(expt);
                    dspec=dspecBuild_lab(dspec,['T' modelLabels{varsToInclude}],numNeuron);
                    dmDummy=buildGLM.compileSparseDesignMatrix(dspec, 1);
                    numCol=size(dmDummy.X,2);
                    clear dmDummy
                    wmlMat{numVars}{varComb}=nan(numFolds,numCol+1);
                    hessian{numVars}{varComb}=nan(numFolds,numCol+1,numCol+1);
                    wmlMean{numVars}{varComb}=nan(1,numCol+1);
                    testFit{numVars}{varComb}=nan(numFolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of test data, AIC
                    trainFit{numVars}{varComb}=nan(numFolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of train data, AIC
                    continue
                end
                disp(['Calculating models with ' num2str(numVars) ' variable(s) (out of ' num2str(totalVars) '): going through combination #' num2str(count) ' (out of ' num2str(length(find(included))) ')'])
                disp(['Testing model ' [modelLabels{varsToInclude}]])
                count=count+1;
                
                clear dspec dm dmTest
                dspec = buildGLM.initDesignSpec(expt);
                
                dspec=dspecBuild_lab(dspec,['T' modelLabels{varsToInclude}],numNeuron);
                
                dmDummy=buildGLM.compileSparseDesignMatrix(dspec, 1);
                numCol=size(dmDummy.X,2);
                clear dmDummy
                wmlMat{numVars}{varComb}=nan(numFolds,numCol+1);
                hessian{numVars}{varComb}=nan(numFolds,numCol+1,numCol+1);
                wmlMean{numVars}{varComb}=nan(1,numCol+1);
                testFit{numVars}{varComb}=nan(numFolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of test data, AIC
                trainFit{numVars}{varComb}=nan(numFolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of train data, AIC
                
                 % ESTIMATE RIDGE REGRESSION HYPERPARAMETER rho ON ALL TRIALS THROUGH EVIDENCE MAXIMIZATION
                dmValidation=buildGLM.buildMyDesingMatrix(dspec,['T' modelLabels{varsToInclude}],1:numTrials,parametrizationParams,0,1);
                yValidation = buildGLM.getBinnedSpikeTrain(dspec.expt, ['sptrain' num2str(numNeuron)], 1:numTrials);
                wInit=(dmValidation.X'*dmValidation.X + eye(size(dmValidation.X,2)))\(dmValidation.X'*yValidation); % Initial guess is regularized least squares with rho=1
                nlfun = @nlfuns.exp;
                rhoGrid=[1 10 100];
                [wRidgeValidation,rho{numVars}{varComb}] = autoRegress_PoissonRidge(dmValidation.X,yValidation,nlfun,2:(size(dmValidation.X,2)),.1,rhoGrid,wInit);
                                
                % RUN RIDGE REGRESSION WITH PARAMETER rho ON ALL TRAIN AND TEST FOLDS
                for k = 1 :numFolds
                    disp(['Fold ' num2str(k)])
                    
                    if k==1
                        wInit=[];
                    else
                        temp=wmlMat{numVars}{varComb}(k-1,:);
                        wInit=wmlMat{numVars}{varComb}(k-1,~isnan(temp));
                    end
                    [wmlMat{numVars}{varComb}(k,:),trainFit{numVars}{varComb}(k,:),testFit{numVars}{varComb}(k,:)]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,[modelLabels{varsToInclude}],wInit,k,numFolds,factor,rho{numVars}{varComb});
                    
                end
                wmlMean{numVars}{varComb} = nanmean(wmlMat{numVars}{varComb});
                
            end % End of model combinations loop
            
            %% 
            % LOOK FOR BEST MODEL WITH numVars VARIABLES, AND STOP IF THERE WAS NO IMPROVEMENT WITH RESPECT TO THE BEST MODEL WITH numVars-1 VARIABLES
            if numVars==1
                temp=cell2mat(testFit{numVars}');
                
                startComb=1:numFolds:size(temp,1);
                endComb=numFolds:numFolds:size(temp,1)+numFolds-1;
                
                meanLH=nan(length(startComb),1);
                for i=1:length(startComb)
                    meanLH(i)=nanmean(temp(startComb(i):endComb(i),3));
                end
                [maxVal,maxComb(numVars)]=max(meanLH);
                
                p_llhTest=nan(length(startComb),1);
                for i=1:length(startComb)
                    [p_llhTest(i),~] =signrank(temp(startComb(i):endComb(i),3),zeros(numFolds,1),'tail','right');
                end
                [minVal,maxComb(numVars)]=nanmin(p_llhTest);
                
                LLH(numVars,:)=temp(startComb(maxComb(numVars)):endComb(maxComb(numVars)),3);
                AIKs(numVars,:)=temp(startComb(maxComb(numVars)):endComb(maxComb(numVars)),7);
                
                bestModel{numVars}=modelType{numVars}(maxComb(numVars),:);
                
                if minVal > 0.05
                    selectInd=[];
                    selected_model = [];
                    clear dspec dm dmTest
                    dspec = buildGLM.initDesignSpec(expt);
                    dspec=dspecBuild_lab(dspec,['T' modelLabels{selected_model}],numNeuron);
                    dm=buildGLM.buildMyDesingMatrix(dspec,['T'  modelLabels{selected_model}],1,parametrizationParams,1,1);
                    disp(['Model could not be fit to any variable'])
                    break
                else
                    disp(['Temporary model selection for neuron ' num2str(numNeuron) ': ' [modelLabels{bestModel{1}}]])
                end
                
            else
                totalCombs=size(modelType{numVars},1);
                
                temp=[];
                for kk=1:totalCombs
                    if isempty(testFit{numVars}{kk})
                        temp=[temp; zeros(numFolds,7)];
                    else
                        temp=[temp; testFit{numVars}{kk}];
                    end
                    
                end
                
                startComb=1:numFolds:size(temp,1);
                endComb=numFolds:numFolds:size(temp,1)+numFolds-1;
             
                indModels=find(included);
                p_llhTest=nan(length(startComb),1);
                for i=1:length(indModels)
                    llhNew=temp(startComb(indModels(i)):endComb(indModels(i)),3);
                    llhOld=LLH(numVars-1,:)';
                    includeFold=~isoutlier(llhNew-llhOld);
                    [p_llhTest(indModels(i)),~] =signrank(llhNew(includeFold),llhOld(includeFold),'tail','right');
                end
                [minVal,maxComb(numVars)]=nanmin(p_llhTest);
                
                LLH(numVars,:)=temp(startComb(maxComb(numVars)):endComb(maxComb(numVars)),3);
                AIKs(numVars,:)=temp(startComb(maxComb(numVars)):endComb(maxComb(numVars)),7);
                bestModel{numVars}=modelType{numVars}(maxComb(numVars),:);
                
                % STATISTICAL COMPARISON
                llhNew=LLH(numVars,:);
                llhOld=LLH(numVars-1,:);
                includeFold=~isoutlier(llhNew-llhOld);
                [p_llh(numVars),~] = signrank(llhNew(includeFold),llhOld(includeFold),'tail','right');
                
                if p_llh(numVars) < 0.05
                    selectInd=numVars;
                    selected_model = bestModel{selectInd};
                    if numVars == totalVars
                        clear dspec dm dmTest
                        dspec = buildGLM.initDesignSpec(expt);
                        dspec=dspecBuild_lab(dspec,['T' modelLabels{selected_model}],numNeuron);
                        dm=buildGLM.buildMyDesingMatrix(dspec,['T'  modelLabels{selected_model}],1,parametrizationParams,1,1);
                        disp(['Stopping optimization at model with ' num2str(numVars) ' variables'])
                        
                        disp(['Running final fit with all trials...'])
                        temp=wmlMat{selectInd}{maxComb(selectInd)}(1,:);
                        wInit=wmlMat{selectInd}{maxComb(selectInd)}(1,~isnan(temp));
                        [wmlFit,~,~]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['T'  modelLabels{selected_model}],wInit,0,0,0,rho{selectInd}{maxComb(selectInd)});
                        
                    else
                        disp(['Temporary model selection for neuron ' num2str(numNeuron) ': ' [modelLabels{selected_model}]])
                    end
                else
                    selectInd=numVars-1;
                    selected_model = bestModel{selectInd};
                    clear dspec dm dmTest
                    dspec = buildGLM.initDesignSpec(expt);
                    dspec=dspecBuild_lab(dspec,['T' modelLabels{selected_model}],numNeuron);
                    dm=buildGLM.buildMyDesingMatrix(dspec,['T'  modelLabels{selected_model}],1,parametrizationParams,1,1);
                    disp(['Stopping optimization at model with ' num2str(numVars) ' variables'])
                    
                    disp(['Running final fit with all trials...'])
                    temp=wmlMat{selectInd}{maxComb(selectInd)}(1,:);
                    wInit=wmlMat{selectInd}{maxComb(selectInd)}(1,~isnan(temp));
                    [wmlFit,~,~]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['T'  modelLabels{selected_model}],wInit,0,0,0,rho{selectInd}{maxComb(selectInd)});
                    
                    break
                end
            end % End of model selection loop
            %% 
            
        end % End of numVars loop
        toc
        
        disp('%%%%%%%%%')
        disp(['RESULT: For neuron' num2str(numNeuron) ' the selected model is ' [modelLabels{selected_model}] ])
        disp('%%%%%%%%%')
        
        cut=findstr(fnRaw,'.mat');
        if ~ isfolder([ fnRaw(1:cut-1) '_modelFit'])
            mkdir(sprintf([ fnRaw(1:cut-1) '_modelFit']',k))
        end
        if isempty(selected_model)
            wmlFit=[];
        end
        dm.dspec=rmfield(dm.dspec,'expt'); % REMOVE expt TO REDUCE FILE SIZE
        save([ fnRaw(1:cut-1) '_modelFit/fit_Neuron' sprintf('%02i', numNeuron) '_' [modelLabels{selected_model}] '.mat'],'wmlFit','wmlMat','wmlMean','testFit','trainFit','rho','totalVars','modelLabels','modelType','numTrials','numFolds','edges','binfun','LLH','AIKs','p_llh','bestModel','maxComb','selectInd','selected_model','dm','parametrizationParams','fitModel')
        
        if ~isempty(selected_model)
            %Adding Spike history and population coupling kernels to the selected model
            disp('%%%%%%%%%')
            disp('Adding Spike history and population coupling kernels to the model')
            disp('%%%%%%%%%')
            
            numVars=selectInd;
            varsToInclude=selected_model;
            clear dspec dm dmTest
            dspec = buildGLM.initDesignSpec(expt);
            
            dspec=dspecBuild_lab(dspec,['THP' modelLabels{varsToInclude}],numNeuron);
            
            dmDummy=buildGLM.compileSparseDesignMatrix(dspec, 1);
            numCol=size(dmDummy.X,2);
            clear dmDummy
            varComb=maxComb(selectInd);
            wmlMat{numVars}{varComb}=nan(numFolds,numCol+1);
            wmlMean{numVars}{varComb}=nan(1,numCol+1);
            testFit{numVars}{varComb}=nan(numFolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of test data, AIC
            trainFit{numVars}{varComb}=nan(numFolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of train data, AIC
            for k = 1 :numFolds
                disp(['Fold ' num2str(k)])
                
                if k==1
                    wInit=[];
                else
                    temp=wmlMat{numVars}{varComb}(k-1,:);
                    wInit=wmlMat{numVars}{varComb}(k-1,~isnan(temp));
                end
                [wmlMat{numVars}{varComb}(k,:),trainFit{numVars}{varComb}(k,:),testFit{numVars}{varComb}(k,:)]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['HP' modelLabels{varsToInclude}],wInit,k,numFolds,factor,rho{numVars}{varComb});
                  
            end
            wmlMean{numVars}{varComb} = nanmean(wmlMat{numVars}{varComb});
            LLH(numVars,:)=testFit{numVars}{varComb}(:,3)';
            AIKs(numVars,:)=testFit{numVars}{varComb}(:,7)';
            
            clear dspec dm dmTest
            dspec = buildGLM.initDesignSpec(expt);
            dspec=dspecBuild_lab(dspec,['THP' modelLabels{selected_model}],numNeuron);
            dm=buildGLM.buildMyDesingMatrix(dspec,['THP'  modelLabels{selected_model}],1,parametrizationParams,1,1);
            
            disp(['Running final fit with all trials...'])
            temp=wmlMat{selectInd}{maxComb(selectInd)}(1,:);
            wInit=wmlMat{selectInd}{maxComb(selectInd)}(1,~isnan(temp));
            [wmlFit,~,~]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['THP'  modelLabels{selected_model}],wInit,0,0,0,rho{selectInd}{maxComb(selectInd)});
            
                        
            disp('%%%%%%%%%')
            disp('Finished')
            disp('%%%%%%%%%')
            
            cut=findstr(fnRaw,'.mat');
            if ~ isfolder([ fnRaw(1:cut-1) '_modelFit/HistAndPop'])
                mkdir(sprintf([ fnRaw(1:cut-1) '_modelFit/HistAndPop']',k))
            end
            dm.dspec=rmfield(dm.dspec,'expt'); % REMOVE expt TO REDUCE FILE SIZE
            save([ fnRaw(1:cut-1) '_modelFit/HistAndPop/fit_Neuron' sprintf('%02i', numNeuron) '_' [modelLabels{selected_model}] '.mat'],'wmlFit', 'wmlMat','wmlMean','testFit','trainFit','rho','totalVars','modelLabels','modelType','numTrials','numFolds','edges','binfun','LLH','AIKs','p_llh','bestModel','maxComb','selectInd','selected_model','dm','parametrizationParams','fitModel')
        end
        
    end % End of neuron number loop, when fitting all models
    %%

else % FIT ONLY fitModel MODEL
    cut=findstr(fnRaw,'.mat');
    fnBase=[ fnRaw(1:cut-1) '_modelFit/fit_Neuron'] ;
    cut=strfind(fnBase,'/');
    if contains(fitModel,'H')
        HPincluded=1;
        Files=dir([fnBase(1:cut(end)) '/HistAndPop/*.mat']);
    else
        HPincluded=0;
        Files=dir([fnBase(1:cut(end)) '/*.mat']);
    end
    
    
    fitModelReal=fitModel;
    for kkk=1:length(neuronsToFit) % RUN THROUGH NEURONS TO FIT
        numNeuron=neuronsToFit(kkk);
        clear fileNum
        for jj=1:length(Files)
            if contains(Files(jj).name,['fit_Neuron' sprintf('%02i', numNeuron)])
                fileNum=jj;
                continue
            end
        end
        if exist('fileNum')
            load([Files(fileNum).folder '/' Files(fileNum).name]);
            fitModel=fitModelReal;
            
            if HPincluded
               fitModel=erase(fitModel,{'H','P'});
            end
               
            if HPincluded
                backUpdir='/oldHyP';
                
            else
                backUpdir='/old';
                
            end
            if ~ isfolder([Files(fileNum).folder backUpdir])
                mkdir([Files(fileNum).folder backUpdir])
            end
            movefile([Files(fileNum).folder '/' Files(fileNum).name],[Files(fileNum).folder backUpdir '/' Files(fileNum).name]);
            
            
        end
        tic;
        
        varsToInclude=nan(length(fitModel),1);
        for jj=1:length(fitModel)
            [~,varsToInclude(jj)]=ismember(fitModel(jj),[modelLabels]);
        end
        varsToInclude=sort(varsToInclude,'ascend')';
        numVars=length(varsToInclude);
        varComb=find(sum(ismember(modelType{numVars},varsToInclude),2)==numVars);
        disp(['Testing model ' [modelLabels{varsToInclude}]])
        clear dspec dm dmTest
        dspec = buildGLM.initDesignSpec(expt);
        dspec=dspecBuild_lab(dspec,['T' modelLabels{varsToInclude}],numNeuron);
  
        
        clear dmDummy rho
              
        
        % ESTIMATE RIDGE REGRESSION HYPERPARAMETER rho ON ALL TRIALS THROUGH EVIDENCE MAXIMIZATION
        dmValidation=buildGLM.buildMyDesingMatrix(dspec,['T' modelLabels{varsToInclude}],1:numTrials,parametrizationParams,0,1);
        yValidation = buildGLM.getBinnedSpikeTrain(dspec.expt, ['sptrain' num2str(numNeuron)], 1:numTrials);
        wInit=(dmValidation.X'*dmValidation.X + eye(size(dmValidation.X,2)))\(dmValidation.X'*yValidation); % Initial guess is regularized least squares with rho=1
        nlfun = @nlfuns.exp;
        rhoGrid=[1 10 100];
        [wRidgeValidation,rho{numVars}{varComb}] = autoRegress_PoissonRidge(dmValidation.X,yValidation,nlfun,2:(size(dmValidation.X,2)),.1,rhoGrid,wInit);
%        rho{numVars}{varComb}=5;
        clear dspec dm dmTest
        dspec = buildGLM.initDesignSpec(expt);
        if HPincluded
            dspec=dspecBuild_lab(dspec,['THP' modelLabels{varsToInclude}],numNeuron);
        else
            dspec=dspecBuild_lab(dspec,['T' modelLabels{varsToInclude}],numNeuron);
        end
        dmDummy=buildGLM.compileSparseDesignMatrix(dspec, 1);
        numCol=size(dmDummy.X,2);
        wmlMat{numVars}{varComb}=nan(numFolds,numCol+1);
        wmlMean{numVars}{varComb}=nan(1,numCol+1);
        
        % RUN RIDGE REGRESSION WITH PARAMETER rho ON ALL TRAIN AND TEST FOLDS 
        for k = 1 :numFolds
            disp(['Fold ' num2str(k)])
            
            if k==1
                wInit=[];
            else
                temp=wmlMat{numVars}{varComb}(k-1,:);
                wInit=wmlMat{numVars}{varComb}(k-1,~isnan(temp));
            end
            if HPincluded
                [wmlMat{numVars}{varComb}(k,:),trainFit{numVars}{varComb}(k,:),testFit{numVars}{varComb}(k,:)]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['HP' modelLabels{varsToInclude}],wInit,k,numFolds,factor,rho{numVars}{varComb});
            else
                [wmlMat{numVars}{varComb}(k,:),trainFit{numVars}{varComb}(k,:),testFit{numVars}{varComb}(k,:)]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,[modelLabels{varsToInclude}],wInit,k,numFolds,factor,rho{numVars}{varComb});
            end
        end
        
        wmlMean{numVars}{varComb} = nanmean(wmlMat{numVars}{varComb});
        LLH(numVars,:)=testFit{numVars}{varComb}(:,3)';
        AIKs(numVars,:)=testFit{numVars}{varComb}(:,7)';
        
        bestModel{numVars}=modelType{numVars}(varComb,:);
        
        % STATISTICAL COMPARISON
        llhNew=LLH(numVars,:);
        if numVars==1
            llhOld=zeros(1,numFolds);
        else
            llhOld=LLH(numVars-1,:);
        end
        includeFold=~isoutlier(llhNew-llhOld);
        [p_llh(numVars),~] = signrank(llhNew(includeFold),llhOld(includeFold),'tail','right');
       
        maxComb(numVars)=varComb;
        selectInd=numVars;
        selected_model = bestModel{selectInd};
        clear dspec dm dmTest
        dspec = buildGLM.initDesignSpec(expt);
        if HPincluded
            dspec=dspecBuild_lab(dspec,['THP' modelLabels{selected_model}],numNeuron);
            dm=buildGLM.buildMyDesingMatrix(dspec,['THP'  modelLabels{selected_model}],1,parametrizationParams,1,1);
            disp(['Running final pass with all trials...'])
            temp=wmlMat{selectInd}{maxComb(selectInd)}(1,:);
            wInit=wmlMat{selectInd}{maxComb(selectInd)}(1,~isnan(temp));
            [wmlFit,~,~]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['THP'  modelLabels{selected_model}],wInit,0,0,0,rho{selectInd}{maxComb(selectInd)});
        else
            dspec=dspecBuild_lab(dspec,['T' modelLabels{selected_model}],numNeuron);
            dm=buildGLM.buildMyDesingMatrix(dspec,['T'  modelLabels{selected_model}],1,parametrizationParams,1,1);
            disp(['Running final pass with all trials...'])
        temp=wmlMat{selectInd}{maxComb(selectInd)}(1,:);
        wInit=wmlMat{selectInd}{maxComb(selectInd)}(1,~isnan(temp));
        [wmlFit,~,~]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['T'  modelLabels{selected_model}],wInit,0,0,0,rho{selectInd}{maxComb(selectInd)});
        end
        
        cut=findstr(fnRaw,'.mat');
        disp('%%%%%%%%%')
        disp(['RESULT: For neuron' num2str(numNeuron) ' the selected model is ' [modelLabels{selected_model}] ])
        disp('%%%%%%%%%')
        
        cut=findstr(fnRaw,'.mat');
        if ~ isfolder([fnRaw(1:cut-1) '_modelFit'])
            mkdir(sprintf([fnRaw(1:cut-1) '_modelFit']',k))
        end
        dm.dspec=rmfield(dm.dspec,'expt'); % REMOVE expt TO REDUCE FILE SIZE
        if HPincluded
            save([ fnRaw(1:cut-1) '_modelFit/HistAndPop/fit_Neuron' sprintf('%02i', numNeuron) '_' [modelLabels{selected_model}] '.mat'],'wmlFit','wmlMat','wmlMean','testFit','trainFit','rho','totalVars','modelLabels','modelType','numTrials','numFolds','edges','binfun','LLH','AIKs','p_llh','bestModel','maxComb','selectInd','selected_model','dm','parametrizationParams','fitModel')
        else
            save([ fnRaw(1:cut-1) '_modelFit/fit_Neuron' sprintf('%02i', numNeuron) '_' [modelLabels{selected_model}] '.mat'],'wmlFit','wmlMat','wmlMean','testFit','trainFit','rho','totalVars','modelLabels','modelType','numTrials','numFolds','edges','binfun','LLH','AIKs','p_llh','bestModel','maxComb','selectInd','selected_model','dm','parametrizationParams','fitModel')
        end
        
        if ~isempty(selected_model) & ~HPincluded
            %Adding Spike history and population coupling kernels to the selected model
            disp('%%%%%%%%%')
            disp('Adding Spike history and population coupling kernels to the model')
            disp('%%%%%%%%%')
            
            numVars=selectInd;
            varsToInclude=selected_model;
            clear dspec dm dmTest
            dspec = buildGLM.initDesignSpec(expt);
            
            dspec=dspecBuild_lab(dspec,['THP' modelLabels{varsToInclude}],numNeuron);
            
            dmDummy=buildGLM.compileSparseDesignMatrix(dspec, 1);
            numCol=size(dmDummy.X,2);
            clear dmDummy
            varComb=maxComb(selectInd);
            wmlMat{numVars}{varComb}=nan(numFolds,numCol+1);
            wmlMean{numVars}{varComb}=nan(1,numCol+1);
            testFit{numVars}{varComb}=nan(numFolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of test data, AIC
            trainFit{numVars}{varComb}=nan(numFolds,7); % var ex, correlation, llh increase, mse, # of spikes, length of train data, AIC
            for k = 1 :numFolds
                disp(['Fold ' num2str(k)])
                
                if k==1
                    wInit=[];
                else
                    temp=wmlMat{numVars}{varComb}(k-1,:);
                    wInit=wmlMat{numVars}{varComb}(k-1,~isnan(temp));
                end
                
                [wmlMat{numVars}{varComb}(k,:),trainFit{numVars}{varComb}(k,:),testFit{numVars}{varComb}(k,:)]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['HP' modelLabels{varsToInclude}],wInit,k,numFolds,factor,rho{numVars}{varComb});
            end
            wmlMean{numVars}{varComb} = nanmean(wmlMat{numVars}{varComb});
            LLH(numVars,:)=testFit{numVars}{varComb}(:,3)';
            AIKs(numVars,:)=testFit{numVars}{varComb}(:,7)';
            
            clear dspec dm dmTest
            dspec = buildGLM.initDesignSpec(expt);
            dspec=dspecBuild_lab(dspec,['THP' modelLabels{selected_model}],numNeuron);
            dm=buildGLM.buildMyDesingMatrix(dspec,['THP'  modelLabels{selected_model}],1,parametrizationParams,1,1);
            
            disp(['Running final fit with all trials...'])
            temp=wmlMat{selectInd}{maxComb(selectInd)}(1,:);
            wInit=wmlMat{selectInd}{maxComb(selectInd)}(1,~isnan(temp));
            [wmlFit,~,~]= fitData2FoldMAPRidge(dspec,parametrizationParams,numNeuron,['THP'  modelLabels{selected_model}],wInit,0,0,0,rho{selectInd}{maxComb(selectInd)});
                    
            
            disp('%%%%%%%%%')
            disp('Finished')
            disp('%%%%%%%%%')
            
            cut=findstr(fnRaw,'.mat');
            if ~ isfolder([ fnRaw(1:cut-1) '_modelFit/HistAndPop'])
                mkdir(sprintf([ fnRaw(1:cut-1) '_modelFit/HistAndPop']',k))
            end
            fnBase=[ fnRaw(1:cut-1) '_modelFit/HistAndPop/fit_Neuron'] ;
            cut=strfind(fnBase,'/');
            Files=dir([fnBase(1:cut(end)) '*.mat']);
            clear fileNum
            for jj=1:length(Files)
                if contains(Files(jj).name,['fit_Neuron' sprintf('%02i', numNeuron)])
                    fileNum=jj;
                    continue
                end
            end
            if exist('fileNum')
                if ~ isfolder([Files(fileNum).folder '/old'])
                    mkdir([Files(fileNum).folder '/old'])
                end
                movefile([Files(fileNum).folder '/' Files(fileNum).name],[Files(fileNum).folder '/old/' Files(fileNum).name]);
            end
            cut=findstr(fnRaw,'.mat');
            dm.dspec=rmfield(dm.dspec,'expt'); % REMOVE expt TO REDUCE FILE SIZE
            save([ fnRaw(1:cut-1) '_modelFit/HistAndPop/fit_Neuron' sprintf('%02i', numNeuron) '_' [modelLabels{selected_model}] '.mat'],'wmlFit','wmlMat','wmlMean','testFit','trainFit','rho','totalVars','modelLabels','modelType','numTrials','numFolds','edges','binfun','LLH','AIKs','p_llh','bestModel','maxComb','selectInd','selected_model','dm','parametrizationParams','fitModel')
        end
        
    end
    %%
end

end