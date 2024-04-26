function [wmlMat,trainFit,testFit]=fitData2FoldMAPRidge(dspec,parametrizationParms,numNeuron,modelVars,wInit,k,numFolds,factor,rho)
% k=0 means fit all trials, no test dataset is created


numTrials=length(dspec.expt.trial);
dt=dspec.expt.binSize/1000;

if k>0
    sections=numFolds*factor; edges=(round(linspace(1,numTrials+1,sections+1)));
    testTrialIndices=[];
    for i=0:factor-1
        testTrialIndices=[testTrialIndices edges(k+i*numFolds):edges(k+i*numFolds+1)-1 ];
    end
    trainTrialIndices = setdiff(1:numTrials,testTrialIndices);
else
    trainTrialIndices = 1:numTrials;
end

if k>0
    dmTest=buildGLM.buildMyDesingMatrix(dspec,['T' modelVars],testTrialIndices,parametrizationParms,0,1);
    yTest = buildGLM.getBinnedSpikeTrain(dspec.expt, ['sptrain' num2str(numNeuron)], testTrialIndices);
end

dm=buildGLM.buildMyDesingMatrix(dspec,['T' modelVars],trainTrialIndices,parametrizationParms,0,1);
y = buildGLM.getBinnedSpikeTrain(dspec.expt, ['sptrain' num2str(numNeuron)], trainTrialIndices);


nlfun = @nlfuns.exp;

wInit=(dm.X'*dm.X + eye(size(dm.X,2)))\(dm.X'*y); % Initial guess is regularized least squares with rho=1
opts = optimoptions(@fminunc, 'Algorithm', 'trust-region', 'GradObj', 'on', 'Hessian','on');
mstruct.neglogli = @neglogli_poiss;
mstruct.logprior = @logprior_ridge;
mstruct.liargs = {dm.X,y,nlfun,dt};
mstruct.priargs = {[2:size(dm.X,2)]',0};
lfpost = @(w)(neglogpost_GLM(w,rho,mstruct));
%[wRidge2, fval, ~, ~, ~, H] = fminunc(lfpost, k0, opts);
[wml, fval, ~, ~, ~, H] = fminunc(lfpost, wInit, opts);

if k>0
    % METRICS FOR TEST DATA
    
    fr_hat_test = exp(dmTest.X * wml);
    fr_test=full(yTest)/dt;
    filter=[0 ones(1,0.1/dt) 0];filter = filter/sum(filter); % 100-ms boxcar filter
    fr_test=conv(fr_test,filter,'same');
    
    % compute variance explained test data
    sse = sum((fr_hat_test-fr_test).^2);
    sst = sum((fr_test-mean(fr_test)).^2);
    varExplain_test = 1-(sse/sst);
    
    % compute correlation with test data
    correlation_test = corr(fr_test,fr_hat_test,'type','Pearson');
    
    % compute llh increase from "mean firing rate model"
    spikecount_hat_test = fr_hat_test*dt;
    spikecount_test= full(yTest);
    log_llh_test_model = nansum( spikecount_test.*log(spikecount_hat_test) - spikecount_hat_test ) / sum(spikecount_test); % Albeit constant independent of filter parameters ( log(factorial(spikecount_test) )
    log_llh_test_mean = nansum( spikecount_test.*log(mean(spikecount_test)) - mean(spikecount_test)) / sum(spikecount_test); % Albeit constant independent of filter parameters ( log(factorial(spikecount_test) )
    log_llh_test = (log_llh_test_model - log_llh_test_mean);
    log_llh_test = log(2)*log_llh_test; % in bits of information
    
    % compute MSE
    mse_test = nanmean((fr_hat_test-fr_test).^2);
    
    % AKAIKE
    AICtest = 2*nansum( spikecount_hat_test - spikecount_test.*log(spikecount_hat_test)) + 2*length(wml);
    
    % fill in all the relevant values for the test fit cases
    testFit = [varExplain_test correlation_test log_llh_test mse_test sum(spikecount_test) numel(testTrialIndices) AICtest];
    
    % METRICS FOR TRAINING DATA
    
    fr_hat_train = exp(dm.X * wml);
    fr_train=full(y)/dt;
    filter=[0 ones(1,0.1/dt) 0];filter = filter/sum(filter); % 100-ms boxcar filter
    fr_train=conv(fr_train,filter,'same');
    
    % compute variance explained train data
    sse = sum((fr_hat_train-fr_train).^2);
    sst = sum((fr_train-mean(fr_train)).^2);
    varExplain_train = 1-(sse/sst);
    
    % compute correlation with train data
    correlation_train = corr(fr_train,fr_hat_train,'type','Pearson');
    
    % compute llh increase from "mean firing rate model"
    spikecount_hat_train = fr_hat_train*dt;
    spikecount_train= full(y);
    log_llh_train_model = nansum( spikecount_train.*log(spikecount_hat_train) - spikecount_hat_train ) / sum(spikecount_train); % Albeit constant independent of filter parameters ( log(factorial(spikecount_train) )
    log_llh_train_mean = nansum( spikecount_train.*log(mean(spikecount_train)) - mean(spikecount_train)) / sum(spikecount_train); % Albeit constant independent of filter parameters ( log(factorial(spikecount_train) )
    log_llh_train = (log_llh_train_model - log_llh_train_mean);
    log_llh_train = log(2)*log_llh_train; % in bits of information
    
    % compute MSE
    mse_train = nanmean((fr_hat_train-fr_train).^2);
    
    % AKAIKE
    AICtrain = 2*nansum( spikecount_hat_train - spikecount_train.*log(spikecount_hat_train)) + 2*length(wml);
    
    % fill in all the relevant values for the train fit cases
    trainFit = [varExplain_train correlation_train log_llh_train mse_train sum(spikecount_train) numel(trainTrialIndices) AICtrain];
    
    % save the parameters
    wmlMat(1,1:length(wml)) = wml;
else
    wmlMat(1,1:length(wml)) = wml;
    trainFit = [];
    testFit = [];
end