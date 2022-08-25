function plot_corrmatrix(r,p,type,partial,Work_dir,write_output)
    % Plot correlation matrix of the MW scores
    
    imagesc(tril(r)); axis square % square image
    
    % colormap
    colormap(brewermap([],'*RdBu')); caxis([-1,1]); cb = colorbar;
    
    % axis labels
    %xlabel('Item'), yl = ylabel('Item');
    %yl.Position(1) = yl.Position(1) - 3.6;
    % colorbar label
    if strcmp(type, 'Pearson')
        notation = 'r';
    elseif strcmp(type, 'Spearman')
        notation = '\rho';
    elseif strcmp(type, 'Kendall')
        notation = '\tau';
    end
    if partial; type = ['Partial ' type]; end
    ylabel(cb,[type ' correlation (' notation ')'],'FontSize',16,'Rotation',270) 
    cb.Label.Position(1) = 4.5; % move label to the right
    
    % tick labels
    xticks(1:12) 
    xticklabels({'Pos.','Neg.','Fut.','Past','Myself','Others','Surr.','Vigil.','Imag.','Words','Spec.','Intru.'})
    %xticklabels({'Positive','Negative','Future','Past','Myself','Others','Surrounding','Vigilance','Images','Words','Specific','Intrusive'})
    xtickangle(90)
    yticks(1:12)
    yticklabels({'Pos.','Neg.','Fut.','Past','Myself','Others','Surr.','Vigil.','Imag.','Words','Spec.','Intru.'})
    %yticklabels({'Positive','Negative','Future','Past','Myself','Others','Surrounding','Vigilance','Images','Words','Specific','Intrusive'})
    
    % add values (and asterisks)
    values = tril(r);
    for x = 1:size(values,1)
        for y = 1:size(values,2)
            value = values(x,y);
            if value ~= 0 && value ~= 1
                if p(x,y) < 0.001
                    text(y,x,[num2str(round(value,2)) '**'],'FontSize',6,'HorizontalAlignment','center');
                elseif p(x,y) < 0.025/size(values,1)
                    text(y,x,[num2str(round(value,2)) '*'],'FontSize',6,'HorizontalAlignment','center');
                else
                    text(y,x,num2str(round(value,2)),'FontSize',6,'HorizontalAlignment','center'); 
                end
            end
        end
    end
    
    % set font sizes
    set(gca,'FontSize',15,'linewidth',1) % for cb ticks (cb label above)
    ax = gca;
    ax.XAxis.FontSize = 12; % x ticks
    ax.YAxis.FontSize = 12; % y ticks
    ax.XLabel.FontSize = 15; % x label
    ax.YLabel.FontSize = 15; % y label
    
    % remove ticks
    ax.XAxis.TickLength = [0 0];
    ax.YAxis.TickLength = [0 0];
    
    % save figure
    if write_output; print(gcf,[Work_dir '/Figures/MW_corr.png'],'-dpng','-r300'); end
    
end
