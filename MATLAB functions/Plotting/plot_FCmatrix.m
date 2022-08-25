function plot_FCmatrix(results,s,Work_dir,write_output) 
    % Plot functional connectivity matrices per state
    
    if isstruct(results); [~,FC] = getFuncConn(results,s); 
    else; FC = results(:,:,s); end
    
    imagesc(FC); axis square % display square image
    
    % colormap
    colormap(brewermap([],'*RdBu')); caxis([-1,1]); cb = colorbar; % set color map
    
    % axis labels
    ylabel(cb,'Pearson correlation (r)','FontSize',18,'Rotation',270) % colorbar label
    cb.Label.Position(1) = 4.5; % move label to the right
    
    set(gca,'FontSize',18,'linewidth',1,'YDir','normal','XTick',[],'YTick',[]) % reverse y axis
    
    % remove ticks
    ax = gca;
    ax.XAxis.TickLength = [0 0];
    ax.YAxis.TickLength = [0 0];
    
    % add network labels a la Yeo 2011 (Table 2.2)
    text(0,4,'VIS','Color',[0.635 0.322 0.678],'FontSize',16,'HorizontalAlignment','right')
    text(1,0,'VIS','Color',[0.635 0.322 0.678],'FontSize',16,'HorizontalAlignment','right','Rotation',90)
    xline(4.5); yline(4.5)
    text(0,8,'SOM','Color',[0.471 0.604 0.753],'FontSize',16,'HorizontalAlignment','right')
    text(5,0,'SOM','Color',[0.471 0.604 0.753],'FontSize',16,'HorizontalAlignment','right','Rotation',90)
    xline(8.5); yline(8.5)
    text(0,12,'DOR','Color',[0.251 0.596 0.196],'FontSize',16,'HorizontalAlignment','right') 
    text(9,0,'DOR','Color',[0.251 0.596 0.196],'FontSize',16,'HorizontalAlignment','right','Rotation',90)
    xline(12.5); yline(12.5)
    text(0,17,'VEN','Color',[1 0.392 0.996],'FontSize',16,'HorizontalAlignment','right')
    text(13,0,'VEN','Color',[0.8 0.4 0.9],'FontSize',16,'HorizontalAlignment','right','Rotation',90)
    xline(17.5); yline(17.5)
    text(0,21,'FRO','Color',[0.941 0.726 0.267],'FontSize',16,'HorizontalAlignment','right')
    text(18,0,'FRO','Color',[0.941 0.726 0.267],'FontSize',16,'HorizontalAlignment','right','Rotation',90)
    xline(21.5); yline(21.5)
    text(0,27,'DEF','Color',[0.851 0.443 0.486],'FontSize',16,'HorizontalAlignment','right')
    text(22,0,'DEF','Color',[0.851 0.443 0.486],'FontSize',16,'HorizontalAlignment','right','Rotation',90)
    
    % save figure
    if write_output
        if isstruct(results); if results.train.id_mixture == 0; model = 'HMM'; else; model = 'GMM'; end
        else; model = 'SWK'; end
        print(gcf,[Work_dir '/Figures/' model '_FC_s' num2str(s) '.png'],'-dpng','-r300'); 
    end
    
end
