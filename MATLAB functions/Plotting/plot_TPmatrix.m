function plot_TPmatrix(hmm,Work_dir,write_output)
    % Plot transition probability matrix
    
    imagesc(getTransProbs(hmm)); axis square % square image
    
    % find maximum prob to set colormap limit
    maxP = max(getTransProbs(hmm),[],'all');
    UL = floor(maxP) + ceil((maxP-floor(maxP))/0.05) * 0.05; % round up to nearest 0.05
    % colormap
    colormap(brewermap([],'OrRd')); caxis([0,UL]); cb = colorbar;
    
    % axis labels
    xlabel('To state'), ylabel('From state')
    % colorbar label
    ylabel(cb,'Transition probability','FontSize',16,'Rotation',270)
    cb.Label.Position(1) = 4.5; % move label to the right
    
    % set font size and tick labels
    set(gca,'FontSize',16,'linewidth',1,'xtick',1:6,'xticklabel',1:6)
    
    % remove ticks
    ax = gca;
    ax.XAxis.TickLength = [0 0];
    ax.YAxis.TickLength = [0 0];
    
    % save figure
    if write_output; if hmm.train.id_mixture == 0; model = 'HMM'; else; model = 'GMM'; end
        print(gcf,[Work_dir '/Figures/' model '_TP.png'],'-dpng','-r300'); end
    
end
