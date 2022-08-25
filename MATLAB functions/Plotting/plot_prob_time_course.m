function plot_prob_time_course(Gamma,period,Work_dir,write_output)
    % Plot state activation probability time course
    % Code adapted from plot_Gamma.m of HMM-MAR
    
    t = length(period); K = size(Gamma,2); 
    
    colors = [0.75 0.18 0.29 % reddish
              0.87 0.73 0.42 % orange
              0.51 0.80 0.29 % green
              0.31 0.66 0.81 % light blue
              0.49 0.18 0.56 % purple
              0 0 1]; % deep blue
    if size(colors,1) < K
        rng(111)
        n = K - size(colors,1);
        colors = [colors; zeros(n,3)];
        for i = (K-n+1):K; colors(i,:) = rand(3,1)'; end
    end
    colors = colors(1:K,:);
    
    set(gca,'ColorOrder',colors,'NextPlot','replacechildren')
    area(Gamma(period,:)); ylim([0 1]); xlim([0 t])
    
    % axis labels
    set(gca,'FontSize',12,'linewidth',1);
    xlabel('Time point'), ylabel('State activation probability')
    
    % add legend
    hold on
    h = zeros(K,1);
    names = char(zeros(K,7));
    for i = 1:K
        h(i) = scatter(NaN,NaN,100,colors(i,:),'square','filled');
        names(i,:) = ['State ' sprintf('%d',i)];
    end
    legend(h,names)
    
    % save figure
    if write_output; print(gcf,[Work_dir '/Figures/HMM_prob_course.png'],'-dpng','-r300'); end
    
end
