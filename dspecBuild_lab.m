function dspec=dspecBuild_lab(dspec,includeInModel,numNeuron)

durationH=200;
durationP=200;
durationI=460;
durationO=460;
durationL=360;
durationR=360;
durationG=2000;
durationM=1500;


if  ismember('H',includeInModel)
    bs = basisFactory.makeNonlinearRaisedCosMine(10, dspec.expt.binSize, [0 durationH], 1);
    dspec = buildGLM.addCovariateSpiketrain(dspec, 'hist', ['sptrain' num2str(numNeuron)], 'History filter',bs);
end

if  ismember('P',includeInModel)
    bs = basisFactory.makeNonlinearRaisedCosMine(10, dspec.expt.binSize, [0 durationP], 1);
    dspec = buildGLM.addCovariateSpiketrain(dspec, 'population', ['popActNot' num2str(numNeuron)], 'Coupling to rest of population',bs);

end

if  ismember('X',includeInModel)
    bs = basisFactory.makeSmoothTemporalBasis('boxcar', 1,1, @(x)(x==0)+ceil(x/1));
    offset=0;
    dspec = buildGLM.addCovariateRaw(dspec, 'positionRew','Position along rewarded corridor', bs,offset,@(trial) (trial.rewCtxt == 1));
    dspec = buildGLM.addCovariateRaw(dspec, 'positionUnrew','Position along unrewarded corridor', bs,offset,@(trial) (trial.rewCtxt == 0));
    
end

if ismember('G',includeInModel)
    downSampleFactor=4;
    nBases=floor(durationG/dspec.expt.binSize)/(2*downSampleFactor)-3;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine',durationG,nBases, dspec.expt.binfun);
    bs.B=bs.B(1:floor(bs.centers(end)),:);
    bs.tr=bs.tr(1:floor(bs.centers(end)),:);
    offset=-floor(bs.tr(end,1))-1;
    dspec = buildGLM.addCovariateTiming(dspec,  'preGO', 'firstLick', 'Modulation before GO response',bs,offset,@(trial) (trial.GO == 1));
end

if ismember('S',includeInModel)
    bs = basisFactory.makeSmoothTemporalBasis('boxcar', 1,1, @(x)(x==0)+ceil(x/1));
    dspec = buildGLM.addCovariateRaw(dspec, 'speed',[],bs);
end

if ismember('I',includeInModel)
    nBases=floor(durationI/dspec.expt.binSize)/2-3;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine',durationI,nBases, dspec.expt.binfun);
    offset=0;
    dspec = buildGLM.addCovariateRaw(dspec, 'inhalations','Mouse inhalations', bs,offset);
end

if ismember('O',includeInModel)
    nBases=floor(durationO/dspec.expt.binSize)/2-3;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine',durationO,nBases, dspec.expt.binfun);
    offset=0;
    dspec = buildGLM.addCovariateRawMatrix(dspec, 'odorRew','Rewarded odor inhalations', bs,offset,@(trial) (trial.rewOdor == 1));
    dspec = buildGLM.addCovariateRawMatrix(dspec, 'odorUnrew','Unrewarded odor inhalations', bs,offset,@(trial) (trial.rewOdor == 0));
    
    nBases=6;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine',durationO,nBases, dspec.expt.binfun);
    offset=0;
    dspec = buildGLM.addCovariateRawMatrix(dspec, 'odorRewAdapt','Rewarded odor inhalations', bs,offset,@(trial) (trial.rewOdor == 1));
    dspec = buildGLM.addCovariateRawMatrix(dspec, 'odorUnrewAdapt','Unrewarded odor inhalations', bs,offset,@(trial) (trial.rewOdor == 0));
    
end

if ismember('M',includeInModel)
    offset=0;
    bs = basisFactory.makeNonlinearRaisedCosMine(20, dspec.expt.binSize, [0 durationM], 100);    
    dspec = buildGLM.addCovariateTiming(dspec,  'modulationRewOdor', 'firstInhalOdorRew', 'Contextual modulation of rewarded odor response when in rewarded context',bs,offset,@(trial) (trial.rewOdor == 1 & trial.rewCtxt==1));
    dspec = buildGLM.addCovariateTiming(dspec,  'modulationUnrewOdor', 'firstInhalOdorUnrew', 'Contextual modulation of unrewarded odor response when in rewarded context',bs,offset,@(trial) (trial.rewOdor == 0 & trial.rewCtxt==1));
end

if ismember('T',includeInModel)
    bs = basisFactory.makeSmoothTemporalBasis('boxcar', 1,1, @(x)(x==0)+ceil(x/1));
    dspec = buildGLM.addCovariateRaw(dspec, 'trialEffect',[],bs);
end

if ismember('L',includeInModel)
    nBases=floor(durationL/dspec.expt.binSize)/2-3;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine',durationL,nBases, dspec.expt.binfun);
    dspec = buildGLM.addCovariateTiming(dspec, 'licks','','Mouse licks',bs);
end

if ismember('R',includeInModel)
    nBases=floor(durationR/dspec.expt.binSize)/2-3;
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine',durationR,nBases, dspec.expt.binfun);
    dspec = buildGLM.addCovariateTiming(dspec, 'reward','','Mouse receives reward',bs);
end

