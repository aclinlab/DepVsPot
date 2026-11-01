function rasterizeScatterLayer(ax, dpi)
% Replace scatter objects in AX with an opaque RGB raster image.
% Capture via PRINT (full canvas, always RGB) — fixes both the fill-to-box
% crop bug and the yellow-square colormap bug from exportgraphics.
    if nargin < 2, dpi = 600; end
    sc = findobj(ax, 'Type', 'scatter');
    if isempty(sc), return; end

    xl = ax.XLim; yl = ax.YLim;
    xd = ax.XDir; yd = ax.YDir;

    oldU = ax.Units; ax.Units = 'pixels';
    pos = ax.Position; ax.Units = oldU;
    wpx = max(round(pos(3)),100); hpx = max(round(pos(4)),100);

    tf = figure('Visible','off','Units','pixels', ...
                'Position',[100 100 wpx hpx],'Color','w', ...
                'PaperPositionMode','auto');   % <-- no paper margins
    ta = axes(tf,'Units','normalized','Position',[0 0 1 1],'Color','w');
    copyobj(sc, ta);
    ta.Colormap = ax.Colormap; ta.CLim = ax.CLim;
    ta.XLim = xl; ta.YLim = yl; ta.XDir = xd; ta.YDir = yd;
    ta.DataAspectRatioMode = 'auto';           % stretch to fill [0 0 1 1]
    axis(ta,'off');

    img = print(tf, '-RGBImage', sprintf('-r%d', dpi));  % in-memory, always H×W×3 RGB
    close(tf);
    if strcmp(yd,'normal'),  img = flipud(img); end
    if strcmp(xd,'reverse'), img = fliplr(img); end

    delete(sc);
    h = image(ax, 'XData', xl, 'YData', yl, 'CData', img);
    uistack(h, 'bottom');

    ax.XLim = xl; ax.YLim = yl; ax.XDir = xd; ax.YDir = yd;
end