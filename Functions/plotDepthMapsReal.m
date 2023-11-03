function plotDepthMapsReal(g1Depth)
    I = length(g1Depth);
    subplot = @(m,n,p) subtightplot (m, n, p, [0.0001 0.01], [0.05 0.05], [0.01 0.01]);
    Filt={'Guidance','Init. Depth','JBF','MRFF','JBU','NAF','LSF','WMF','GIF','GAD','CWMF','FWMF','WGIF','GGIF','MSJF','FWL1','FWLS','SDF','EGIF','MGSD','MGDD'};
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i=1:I
        subplot(3,7,i);
        imshow(g1Depth{i},[]);
        colormap gray;
        grid off; axis off;
        title([Filt{i}],'FontSize', 22);
        h = get(gca, 'Title');
        set(h, 'Units', 'normalized');
    end
end