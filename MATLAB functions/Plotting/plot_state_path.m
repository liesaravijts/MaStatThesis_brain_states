function plot_state_path(state_path,K,period,Work_dir,write_output)
    % Plot (most likely) state path in step plot
    % Code adapted from https://stackoverflow.com/questions/38843676/adding-colors-to-lines-of-a-stairstep-plot
    
    x = 1:length(period);
    y = state_path(period,:);
    states = 1:K;
    
    % set colors
    colors = [0.75 0.18 0.29 % reddish
              0.87 0.73 0.42 % orange
              0.51 0.80 0.29 % green
              0.31 0.66 0.81 % light blue
              0.49 0.18 0.56 % purple
              0 0 1]; % deep blue
    if size(colors,1) < K
        rng(111)
        n = length(states) - size(colors,1);
        colors = [colors; zeros(n,3)];
        for i = (length(states)-n+1):length(states); colors(i,:) = rand(3,1)'; end
    end
    colors = colors(1:K,:);
    
    idx = find([1;diff(y)]);  % find state change
    idx(end+1) = length(x)+1; % add last point
    
    % left end state
    k = 1;
    % find current state level
    c = states == y(idx(k));
    % plot bold line
    plot(x([idx(k),idx(k+1)-1]),y(idx(k))*ones(2,1),'color',colors(c,:),'linewidth',5);
    hold on
    for k = 2 : length(idx)-1
        % find current state level
        c = states == y(idx(k));
        % plot dashed line from left state to current state
        plot(x([idx(k)-1,idx(k)]),[y(idx(k-1));y(idx(k))],'color',[0 0 0],'linewidth',0.5);
        % plot bold line for current state with specified color
        plot(x([idx(k),idx(k+1)-1]),y(idx(k))*ones(2,1),'color',colors(c,:),'linewidth',5);
    end
    xlim([0 length(x)]); ylim([1 K])
    
    % axis labels
    set(gca,'FontSize',12,'linewidth',1);
    xlabel('Time point'), ylabel('State active')
    
    % add legend
    h = zeros(max(states),1);
    names = char(zeros(max(states),7));
    for i = 1:max(states)
        h(i) = scatter(NaN,NaN,100,colors(i,:),'square','filled');
        names(i,:) = ['State ' sprintf('%d',i)];
    end
    legend(h,names)
    
    % save figure
    if write_output; print(gcf,[Work_dir '/Figures/HMM_state_path.png'],'-dpng','-r300'); end
    
end
