%% grid of parameter densities

xlabels = {'\mu_{enter,pre}', '\sigma_{enter,pre}', 'm_{\mu,stay,pre}'...
    'm_{\sigma,stay,pre}','\mu_s(V=15)', '\mu_s(V=30)', '\mu_s(V=60)',...
    'H_0(V=0)','H_0(V=15)','H_0(V=30)','H_0(V=60)', 'H',...
    'T_a','T_m','\sigma_s(V=15)','\sigma_s(V=30)','\sigma_s(V=60)',...
    'unstable variability'};
figure('Position',[100 100 600 600])
tlo = tiledlayout(D, D, 'TileSpacing', 'none', 'Padding', 'none');

%plot prob densities along diagonal
xlimits = zeros(D,2);
for i=1:D
    nexttile(i+(i-1)*D);
    [f,xi] = ksdensity(thetas_all_linear{1}(:,i),'Weights',weights_all_linear{1},'BoundaryCorrection','reflection');
    plot(xi,f,'Color','k','LineWidth',0.5)
    xlimits(i,:) = xlim; % this will set it for the rest of the column
    axis square
    set(gca,'YTickLabel',[])
    % if i~=D
        set(gca,'XTickLabel',[])
    % end
end
for i=1:D
    for j=(i+1):D
        disp([i j])
        nexttile(i+(j-1)*D);
        binscatter(thetas_all_linear{1}(:,i),thetas_all_linear{1}(:,j),[25 25],'XLimits',xlimits(i,:),'YLimits',xlimits(j,:))
        colormap(parula);
        colorbar off
        axis square
        set(gca,'YTickLabel',[])
        set(gca,'XTickLabel',[])
        % if j==D
        %     xlabel(xlabels{i})
        % end
        % ax = gca;
        % if i~=1
        %     set(gca,'YTickLabel',[])
        %     ax.YAxis.TickLength = [0 0];
        % else
        %     ax.YAxis.TickLength = [0.03 0.03];
        % end
        % if j~=D
        %     set(gca,'XTickLabel',[])
        %     ax.XAxis.TickLength = [ 0 0];
        % else
        %     ax.XAxis.TickLength = [0.03 0.03];
        % end
        % set(gca,'TickDir','out')

    end
end
% fontsize(6,'points')
exportgraphics(gcf, 'correlationFig.pdf', 'ContentType', 'vector')

