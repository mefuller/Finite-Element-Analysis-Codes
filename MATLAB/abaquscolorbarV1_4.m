function cbar = abaquscolorbarV1_4(Data,BottomLimit,TopLimit,levels)
%20200629 - Condensed function into a single algorithm to handle all cases.
% Tick scaling and labeling have been fixed.  Additional inputs for
% specified limits added to this function and contourPlot.  If limits are
% not supplied, the min and max of the data range is used.

    top2 = max(max(Data));
    bot2 = min(min(Data));
    if isempty(BottomLimit) == 1
        botlim = bot2;
    else
        botlim = BottomLimit;
    end
    if isempty(TopLimit) == 1
        toplim = top2;
    else
        toplim = TopLimit;
    end
    tickw = (toplim-botlim)/(levels);

    if bot2 < botlim && toplim < top2
       cmap_bt = colormap([0.15 0.15 0.15;jet(levels);.85 .85 .85]);
       cbar = colorbar;
       tickw = (toplim-botlim)/(levels);
       cbar.Ticks = [botlim-tickw:tickw:toplim+tickw];
       caxis([botlim-tickw toplim+tickw]);
       for i = 1:length(cbar.Ticks)
           cbar.TickLabels(i) = {cbar.Ticks(i)};
       end
       cbar.TickLabels(1) = {bot2};
       cbar.TickLabels(levels+3) = {top2};
    elseif bot2 < botlim && toplim >= top2
       cmap_b = colormap([0.15 0.15 0.15;jet(levels)]);
       cbar = colorbar;
       tickw = (toplim-botlim)/(levels);
       cbar.Ticks = [botlim-tickw:tickw:toplim];
       caxis([botlim-tickw toplim]);
       for i = 1:length(cbar.Ticks)
           cbar.TickLabels(i) = {cbar.Ticks(i)};
       end
       cbar.TickLabels(1) = {bot2};
    elseif bot2 >= botlim && toplim < top2
       cmap_t = colormap([jet(levels);.85 .85 .85]);
       cbar = colorbar;
       tickw = (toplim-botlim)/(levels);
       cbar.Ticks = [botlim:tickw:toplim+tickw];
       caxis([botlim toplim+tickw]);
       for i = 1:length(cbar.Ticks)
           cbar.TickLabels(i) = {cbar.Ticks(i)};
       end
       cbar.TickLabels(levels+2) = {top2};
    else
       cmap = colormap(jet(levels));
       cbar = colorbar;
       tickw = (toplim-botlim)/(levels);
       cbar.Ticks = [botlim:tickw:toplim];
       caxis([botlim toplim]);
       for i = 1:length(cbar.Ticks)
           cbar.TickLabels(i) = {cbar.Ticks(i)};
       end
    end
    for i = 1:length(cbar.TickLabels)
        %Sets the tick labels to engineering notation.
        cbar.TickLabels{i} = num2str(str2num(cbar.TickLabels{i}),'%1.3e');
    end
end


