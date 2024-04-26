function  wKerns=parseKernels(wsk,parametrizationParams,modelToPlot,plotFlag)

doneO=0;
doneI=0;

% PLOT KERNELS
if  ismember('H',modelToPlot)
    wKerns.hist.axis=wsk.hist.tr;
    wKerns.hist.gain=exp(wsk.hist.data);
    if plotFlag(1)
        if ~plotFlag(2)
        figure; plot(wKerns.hist.axis,wKerns.hist.gain,'lineWidth',1)
        title('History')
        xlabel('Time (ms)')
        set(gcf,'Name','History')
        else
            fh = findobj( 'Type', 'Figure', 'Name', 'History' );
            if isempty(fh)
                figure('Name','History')
                fh = findobj( 'Type', 'Figure', 'Name', 'History' );
                title('History')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh); plot(wKerns.hist.axis,wKerns.hist.gain,'lineWidth',1) 
        end
    end
end

if  ismember('P',modelToPlot)
    wKerns.population.axis=wsk.population.tr;
    wKerns.population.gain=exp(wsk.population.data);
    if plotFlag(1)
        if ~plotFlag(2)
        figure; plot(wKerns.population.axis,wKerns.population.gain,'lineWidth',1)
        title('Population')
        xlabel('Time (ms)')
        set(gcf,'Name','Population')
        else
            fh = findobj( 'Type', 'Figure', 'Name', 'Population' );
            if isempty(fh)
                figure('Name','Population')
                fh = findobj( 'Type', 'Figure', 'Name', 'Population' );
                title('Population')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh); plot(wKerns.population.axis,wKerns.population.gain,'lineWidth',1)
        end
    end
    
end

if ismember('M',modelToPlot)
    wKerns.modulationRewOdor.axis=wsk.modulationRewOdor.tr;
    wKerns.modulationRewOdor.gain=exp(wsk.modulationRewOdor.data);
    wKerns.modulationUnrewOdor.axis=wsk.modulationUnrewOdor.tr;
    wKerns.modulationUnrewOdor.gain=exp(wsk.modulationUnrewOdor.data);
    
    if plotFlag(1)
        if ~plotFlag(2)
        figure; plot(wKerns.modulationRewOdor.axis, wKerns.modulationRewOdor.gain,'lineWidth',1)
        hold on
        plot(wKerns.modulationUnrewOdor.axis, wKerns.modulationUnrewOdor.gain,'lineWidth',1)
        legend('Contextual modulation of rewarded odor in rewarded context','Contextual modulation of unrewarded odor in rewarded context')
        title('Modulation')
        xlabel('Time (ms)')
        set(gcf,'Name','Contextual modulation')
        else
            fh1 = findobj( 'Type', 'Figure', 'Name', 'Modulation Rewarded Odor' );
            fh2 = findobj( 'Type', 'Figure', 'Name', 'Modulation Unrewarded Odor' );
            if isempty(fh1)
                figure('Name','Modulation Rewarded Odor')
                fh1 = findobj( 'Type', 'Figure', 'Name', 'Modulation Rewarded Odor' );
                title('Modulation Rewarded Odor')
                xlabel('Time (ms)')
                hold on
                
                figure('Name','Modulation Unrewarded Odor')
                fh2 = findobj( 'Type', 'Figure', 'Name', 'Modulation Unrewarded Odor' );
                title('Modulation Unrewarded Odor')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh1); plot(wKerns.modulationRewOdor.axis, wKerns.modulationRewOdor.gain,'lineWidth',1)
            figure(fh2); plot(wKerns.modulationUnrewOdor.axis, wKerns.modulationUnrewOdor.gain,'lineWidth',1)
        end
    end
    
end
if ismember('X',modelToPlot)
    wKerns.positionRew.axis=wsk.positionRew.tr*parametrizationParams.binSizePos/parametrizationParams.binSize;
    wKerns.positionRew.gain=exp(wsk.positionRew.data);
    wKerns.positionUnrew.axis=wsk.positionUnrew.tr*parametrizationParams.binSizePos/parametrizationParams.binSize;
    wKerns.positionUnrew.gain=exp(wsk.positionUnrew.data);
    
    if plotFlag(1)
        if ~plotFlag(2)
        figure; plot(wKerns.positionRew.axis,wKerns.positionRew.gain,'lineWidth',1)
        hold on; plot(wKerns.positionUnrew.axis,wKerns.positionUnrew.gain,'lineWidth',1)
        plot([parametrizationParams.threshContext(1) parametrizationParams.threshContext(1)],[min(wKerns.positionRew.gain) max(wKerns.positionRew.gain)],':k')
        plot([parametrizationParams.threshContext(2) parametrizationParams.threshContext(2)],[min(wKerns.positionRew.gain) max(wKerns.positionRew.gain)],':k')
        legend('Pos in rewarded context','Pos in unrewarded context')
        title('Position'); xlabel('Pos along corridor')
        set(gcf,'Name','Position')
        else
            fh1 = findobj( 'Type', 'Figure', 'Name', 'Position in rewarded context' );
            fh2 = findobj( 'Type', 'Figure', 'Name', 'Position in unrewarded context' );
            if isempty(fh1)
                figure('Name','Position in rewarded context')
                fh1 = findobj( 'Type', 'Figure', 'Name','Position in rewarded context' );
                hold on
                plot([parametrizationParams.threshContext(1) parametrizationParams.threshContext(1)],[0.5 2],':k')
                plot([parametrizationParams.threshContext(2) parametrizationParams.threshContext(2)],[0.5 2],':k')
                title('Position in rewarded context')
                xlabel('Pos along corridor')
               
                
                figure('Name','Position in unrewarded context')
                fh2 = findobj( 'Type', 'Figure', 'Name', 'Position in unrewarded context' );
                hold on
                plot([parametrizationParams.threshContext(1) parametrizationParams.threshContext(1)],[0.5 2],':k')
                plot([parametrizationParams.threshContext(2) parametrizationParams.threshContext(2)],[0.5 2],':k')
                title('Position in unrewarded context')
                xlabel('Pos along corridor')
               
            end
            figure(fh1); plot(wKerns.positionRew.axis,wKerns.positionRew.gain,'lineWidth',1)
            figure(fh2); plot(wKerns.positionUnrew.axis,wKerns.positionUnrew.gain,'lineWidth',1)
        end
    end
    
end

if ismember('S',modelToPlot)
    wKerns.speed.axis=wsk.speed.tr*parametrizationParams.binSizeSpeed/parametrizationParams.binSize;
    wKerns.speed.gain=exp(wsk.speed.data);
  
   
    if plotFlag(1)
        if ~plotFlag(2)
        figure; plot(wKerns.speed.axis,wKerns.speed.gain,'lineWidth',1)
        title('Speed'); xlabel('Cm per seg')
        set(gcf,'Name','Speed')
        else
            fh = findobj( 'Type', 'Figure', 'Name', 'Speed' );
            if isempty(fh)
                figure('Name','Speed')
                fh = findobj( 'Type', 'Figure', 'Name', 'Speed' );
                title('Speed')
                xlabel('Cm per seg')
                hold on
            end
            figure(fh); plot(wKerns.speed.axis,wKerns.speed.gain,'lineWidth',1)
        end
    end
    
end

if ismember('O',modelToPlot)
    wKerns.odorRew.axis=wsk.odorRew.tr;
    wKerns.odorRew.gain=exp(wsk.odorRew.data);
    wKerns.odorUnrew.axis=wsk.odorUnrew.tr;
    wKerns.odorUnrew.gain=exp(wsk.odorUnrew.data);
    
    wKerns.odorRewAdapt.axis=wsk.odorRewAdapt.tr;
    wKerns.odorRewAdapt.gain=exp(wsk.odorRewAdapt.data);
    wKerns.odorUnrewAdapt.axis=wsk.odorUnrewAdapt.tr;
    wKerns.odorUnrewAdapt.gain=exp(wsk.odorUnrewAdapt.data);
    
     if plotFlag(1)
        if ~plotFlag(2)
        figure; plot(wKerns.odorRew.axis,wKerns.odorRew.gain,'linewidth',1)
        hold on
        plot(wKerns.odorUnrew.axis,wKerns.odorUnrew.gain,'linewidth',1)
        legend('Rewarded odor','Unrewarded odor')
        title('Odors'); xlabel('Time (ms)')
        set(gcf,'Name','Odors')
        
        
        cmap=colormap(cbrewer('seq', 'Reds',1.5*size(wKerns.odorRewAdapt.gain,2)));
        figure;
        subplot(2,1,1)
        hold on
        for i=1:size(wKerns.odorRewAdapt.gain,2)
            plot(wKerns.odorRewAdapt.axis(:,1),wKerns.odorRewAdapt.gain(:,i),'linewidth',1,'color',cmap(.5*size(wKerns.odorRewAdapt.gain,2)+i,:))
            
        end
        title('Rewarded Odor Adaptation'); xlabel('Time (ms)');
        subplot(2,1,2)
        hold on
        for i=1:size(wKerns.odorUnrewAdapt.gain,2)
            plot(wKerns.odorUnrewAdapt.axis(:,1),wKerns.odorUnrewAdapt.gain(:,i),'linewidth',1,'color',cmap(.5*size(wKerns.odorRewAdapt.gain,2)+i,:))
            
        end
        title('Unrewarded Odor Adaptation'); xlabel('Time (ms)');
        set(gcf,'Name','Odor Adaptation')

        else
            fh1 = findobj( 'Type', 'Figure', 'Name', 'Rewarded Odor' );
            fh2 = findobj( 'Type', 'Figure', 'Name', 'Unrewarded Odor' );
            if isempty(fh1)
                figure('Name','Rewarded Odor')
                fh1 = findobj( 'Type', 'Figure', 'Name', 'Rewarded Odor' );
                title('Rewarded Odor')
                xlabel('Time (ms)')
                hold on
                
                figure('Name','Unrewarded Odor')
                fh2 = findobj( 'Type', 'Figure', 'Name', 'Unrewarded Odor' );
                title('Unrewarded Odor')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh1); plot(wKerns.odorRew.axis,wKerns.odorRew.gain,'linewidth',1)
            figure(fh2); plot(wKerns.odorUnrew.axis,wKerns.odorUnrew.gain,'linewidth',1)
        end
    end
end

if ismember('I',modelToPlot)
    wKerns.inhalations.axis=wsk.inhalations.tr;
    wKerns.inhalations.gain=exp(wsk.inhalations.data);
    
    if plotFlag(1)
        if ~plotFlag(2)
        figure;
        plot(wKerns.inhalations.axis,wKerns.inhalations.gain,'linewidth',1)
        
        title('Inhalations'); xlabel('Time (ms)')
        set(gcf,'Name','Inhalations')
        else
            fh = findobj( 'Type', 'Figure', 'Name', 'Inhalations' );
            if isempty(fh)
                figure('Name','Inhalations')
                fh = findobj( 'Type', 'Figure', 'Name', 'Inhalations' );
                title('Inhalations')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh);  plot(wKerns.inhalations.axis,wKerns.inhalations.gain,'linewidth',1)
        end
    end
    
end

if ismember('L',modelToPlot)
    wKerns.licks.axis=wsk.licks.tr;
    wKerns.licks.gain=exp(wsk.licks.data);
    
    
    if plotFlag(1)
        if ~plotFlag(2)
        figure;
        plot(wKerns.licks.axis,wKerns.licks.gain,'linewidth',1)
        title('Licking'); xlabel('Time (ms)')
        set(gcf,'Name','Licking')
        else
            fh = findobj( 'Type', 'Figure', 'Name', 'Licking' );
            if isempty(fh)
                figure('Name','Licking')
                fh = findobj( 'Type', 'Figure', 'Name', 'Licking' );
                title('Licking')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh);  plot(wKerns.licks.axis,wKerns.licks.gain,'linewidth',1)
        end
    end
end

if ismember('R',modelToPlot)
    wKerns.reward.axis=wsk.reward.tr;
    wKerns.reward.gain=exp(wsk.reward.data);
    
    if plotFlag(1)
        if ~plotFlag(2)
        figure;
        plot(wKerns.reward.axis,wKerns.reward.gain,'linewidth',1)
        title('Reward'); xlabel('Time (ms)')
        set(gcf,'Name','Reward')
        else
            fh = findobj( 'Type', 'Figure', 'Name', 'Reward' );
            if isempty(fh)
                figure('Name','Reward')
                fh = findobj( 'Type', 'Figure', 'Name', 'Reward' );
                title('Reward')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh);  plot(wKerns.reward.axis,wKerns.reward.gain,'linewidth',1)
        end
    end
end

if ismember('T',modelToPlot)
    wKerns.trialEffect.axis=wsk.trialEffect.tr/parametrizationParams.binSize*parametrizationParams.binSizeTrialEffect;
    wKerns.trialEffect.gain=exp(wsk.trialEffect.data);
    
     if plotFlag(1)
        if ~plotFlag(2)
        figure;
        plot(wKerns.trialEffect.axis,wKerns.trialEffect.gain,'linewidth',1)
        title('Trial effect'); xlabel('Trial #')
        set(gcf,'Name','Trial effect')
        else
            fh = findobj( 'Type', 'Figure', 'Name', 'Trial effect' );
            if isempty(fh)
                figure('Name','Trial effect')
                fh = findobj( 'Type', 'Figure', 'Name', 'Trial effect' );
                title('Trial effect')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh);  plot(wKerns.trialEffect.axis,wKerns.trialEffect.gain,'linewidth',1)
        end
    end
end


if ismember('G',modelToPlot)
    wKerns.preGO.axis=wsk.preGO.tr;
    wKerns.preGO.gain=exp(wsk.preGO.data);
    
    if plotFlag(1)
        if ~plotFlag(2)
        figure;
        plot(wKerns.preGO.axis,wKerns.preGO.gain,'linewidth',1)
        legend('PreGO')
        title('Pre GO response'); xlabel('Time (ms)')
        set(gcf,'Name','PreGO')
        else
            fh = findobj( 'Type', 'Figure', 'Name', 'PreGO' );
            if isempty(fh)
                figure('Name','PreGO')
                fh = findobj( 'Type', 'Figure', 'Name', 'PreGO' );
                title('PreGO')
                xlabel('Time (ms)')
                hold on
            end
            figure(fh);  plot(wKerns.preGO.axis,wKerns.preGO.gain,'linewidth',1)
        end
    end
    
end



