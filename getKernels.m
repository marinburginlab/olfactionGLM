function wKerns=getKernels(fnFit,includeHP,plotFlag)

if includeHP==1
    cut=strfind(fnFit,'/');
    fnFit= [fnFit(1:cut(end)) 'HistAndPop' fnFit(cut(end):end)];
end
cut=strfind(fnFit,'_');
numNeuron=str2num(fnFit(cut(end)-2:cut(end)-1));

load(fnFit);


disp(['The selected model is ' [modelLabels{selected_model}]])


cut=strfind(fnFit,'_modelFit');
load([fnFit(1:cut-1) '_exptData.mat']);
dm.dspec.expt=expt;
clear expt;

numTrials=length(dm.dspec.expt.trial);


dspec = buildGLM.initDesignSpec(dm.dspec.expt);
if contains(fnFit,'HistAndPop')
    dspec=dspecBuild_lab(dspec,['THP' modelLabels{selected_model}],numNeuron);
    dm=buildGLM.buildMyDesingMatrix(dspec,['THP' modelLabels{selected_model}],1:numTrials,parametrizationParams,1,1);
else
    dspec=dspecBuild_lab(dspec,['T' modelLabels{selected_model}],numNeuron);
    dm=buildGLM.buildMyDesingMatrix(dspec,['T' modelLabels{selected_model}],1:numTrials,parametrizationParams,1,1);
end

wsk = buildGLM.combineWeights(dm, wmlFit);
   

if contains(fnFit,'HistAndPop')
    wKerns=parseKernels(wsk,parametrizationParams,['HP' modelLabels{selected_model}],plotFlag);
else
    wKerns=parseKernels(wsk,parametrizationParams,[modelLabels{selected_model}],plotFlag);
end