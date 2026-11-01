function [] = plotPerfVsOverlap(overlaps, performance, fractionKCsActive, maxActivity, whichOdor)
hold on
total_rescaled = rescale(fractionKCsActive,1,256*max(fractionKCsActive,[],'all')/maxActivity);
[~,reorder] = sort(fractionKCsActive,'descend');
cmap = colormap;
for odor1=reorder(:).' % force reorder to be a row vector
    [~,odor2_reordered] = sort(overlaps(odor1,:));
    odor2_reordered(odor2_reordered==odor1) = [];
    % plot(overlaps(odor1,odor2_reordered),performance(odor1,odor2_reordered),'Color',[cmap(ceil(total_rescaled(odor1)),:), 0.2])
    switch whichOdor
        case 1
            thisOdor = odor1;
        case 2
            thisOdor = odor2_reordered;
    end

    scatter(overlaps(odor1,odor2_reordered),performance(odor1,odor2_reordered),1.5,cmap(ceil(total_rescaled(thisOdor)),:),'filled','MarkerFaceAlpha',0.1)

end

end
