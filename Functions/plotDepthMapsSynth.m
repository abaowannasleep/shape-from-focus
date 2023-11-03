function plotDepthMapsSynth(g1Depth,D,p1,p2)
    ofst = 5; % offset for trimming 5 boundry pixels
    D=D(ofst+1:end-ofst,ofst+1:end-ofst);
    I = length(g1Depth);
    subplot = @(m,n,p) subtightplot (m, n, p, [0.0001 0.01], [0.05 0.05], [0.01 0.01]);
    Filt={'Init. Depth','GT','JBF','MRFF','JBU','NAF','LSF','WMF','GIF','GAD','CWMF','FWMF','WGIF','GGIF','MSJF','FWL1','FWLS','SDF','EGIF','MGSD','MGDD'};
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i=1:I
        subplot(3,7,i);
        if i==1
            mesh(g1Depth{2});
        elseif i==2
            mesh(D);
        else
            mesh(g1Depth{i});
        end
            colormap jet;
            grid off; axis off; axis tight;
            title([Filt{i}],'FontSize', 25);
            h = get(gca, 'Title');
            set(h, 'Units', 'normalized');
            set(h, 'Position', [p1,p2]);
    end
end