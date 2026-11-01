colors = [0 0 0;
    0.102, 0.373, 0.478;
    0.739, 0.390, 0.571;
    1, 0.647, 0];

% based on the results of Fig 7
Tas = [0.029, 0.0215];
Tms = [0, 0.1710];
x=0:0.01:0.25;

figure
% j=1: Tm=0
% j=2: Tm>0
for j=1:2
    Ta = Tas(j);
    Tm = Tms(j);
    subplot(1,2,j)
    ss = [0 0.3 0.4 0.6]; % no shock, weak shock, medium shock, strong shock
    hold on
    for i=1:length(ss)
        plot(x,normcdf(ss(i).*x./sqrt( (Ta + Tm*x).^2 + (Ta + Tm*(1-ss(i))*x).^2 )),'Color',colors(i,:),'LineWidth',1)
    end
    xlim([0 0.25])
    ylim([0.5 1])
    set(gca,'linewidth',0.5)
    ylabel('Learned P(avoid)')
    xlabel('Fraction KCs responding')
    axis square
    xticks([])
    yticks([])
    set(gca,'units','centimeters','Position',[1+(j-1)*4,1,2.5,2.5],'FontSize',6,'FontName','Arial')
end

