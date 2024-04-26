function dm=buildMyDesingMatrix(dspec,vars,trialIndices,params,remConst,addBias)
% BUILD DESIGN MATRIX
dm = buildGLM.compileSparseDesignMatrix(dspec, trialIndices);
if  ismember('X',vars)
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', params.maxPos, params.n_pos_bins, @(x)(x==0)+ceil(x/params.binSizePos));
    nCols=size(bs.B,1);
    a=zeros(nCols,nCols);
    for i=1:nCols
        temp=zeros(1,nCols);
        temp(i)=1;
        a(i,:)=conv(temp,gausswin(5),'same');
        
    end
    bs.B=a;
    bs.B=bs.B./sum(bs.B,2);
    
    ind=find(strcmp({dm.dspec.covar.label},'positionRew'));
    dm.dspec.covar(ind).basis=bs;
    
    ind=find(strcmp({dm.dspec.covar.label},'positionUnrew'));
    dm.dspec.covar(ind).basis=bs;
    
end

if  ismember('S',vars)
    %bs = basisFactory.makeSmoothTemporalBasis('raised cosine', params.maxSpeed+1, params.n_speed_bins, @(x)(x==0)+ceil(x/(params.binSizeSpeed)));
    bs = basisFactory.makeSmoothTemporalBasis('raised cosine', floor(params.maxSpeed) +1, params.n_speed_bins, @(x)(x==0)+ceil(x/(params.binSizeSpeed)));
    
    nCols=size(bs.B,1);
    a=zeros(nCols,nCols);
    for i=1:nCols
        temp=zeros(1,nCols);
        temp(i)=1;
        a(i,:)=conv(temp,gausswin(5),'same');
        
    end
    bs.B=a;
    bs.B=bs.B./sum(bs.B,2);
    
    indSpeed=find(strcmp({dm.dspec.covar.label},'speed'));
    dm.dspec.covar(indSpeed).basis=bs;
    
end

if  ismember('T',vars)
    bs = basisFactory.makeSmoothTemporalBasis('boxcar', params.numTrials, params.nBinsTrialEffect, @(x)(x==0)+ceil(x/(params.binSizeTrialEffect)));
    
    nCols=size(bs.B,1);
    a=zeros(nCols,nCols);
    for i=1:nCols
        temp=zeros(1,nCols);
        temp(i)=1;
        a(i,:)=conv(temp,gausswin(3),'same');
        
    end
    bs.B=a;
    bs.B=bs.B./sum(bs.B,2);
    
    indTrEff=find(strcmp({dm.dspec.covar.label},'trialEffect'));
    dm.dspec.covar(indTrEff).basis=bs;
    
end
if  ismember('P',vars)
    lims=cumsum([dm.dspec.covar.edim]);
    indPop=find(strcmp({dm.dspec.covar.label},'population'));
    if indPop==1
        dm.X(:,1:lims(indPop))=dm.X(:,1:lims(indPop))/(params.totalNeurons-1);
    else
        dm.X(:,lims(indPop-1)+1:lims(indPop))=dm.X(:,lims(indPop-1)+1:lims(indPop))/(params.totalNeurons-1);
    end
end
if remConst
    dm = buildGLM.removeConstantCols(dm);
end
if addBias
    dm = buildGLM.addBiasColumn(dm);
end