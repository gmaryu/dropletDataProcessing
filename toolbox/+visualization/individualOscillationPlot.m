function individualOscillationPlot(oscillationData, frameToMin, savePath)

%%
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

positions = unique(oscillationData.timeSeries.POS_ID)';
ts = oscillationData.timeSeries;
cycles = oscillationData.cycle;
for p = 1:length(positions)
    pos = positions(p);
    
    tmpOsci = ts(ts.POS_ID == pos,:);
    tmpCycs = cycles(cycles.POS_ID == pos, :);
    ids = unique(tmpOsci.TRACK_ID);
    for i = 1:length(ids)
        % temporal signal
        t = tmpOsci(tmpOsci.TRACK_ID==ids(i),:).POSITION_T;
        signal = tmpOsci(tmpOsci.TRACK_ID==ids(i),:).MAIN_SIGNAL;

        % peaks and troughs
        peaks = tmpCycs(tmpCycs.TRACK_ID==ids(i),:).START_FRAME;
        troughs = tmpCycs(tmpCycs.TRACK_ID==ids(i),:).TROUGH_FRAME;

        % signal value at peaks and troughs
        peakIdx = ismember(t, peaks);
        peakT = t(peakIdx);
        peakV = signal(peakIdx);

        troughIdx = ismember(t, troughs);
        troughT = t(troughIdx);
        troughV = signal(troughIdx);

        % Create an invisible figure.
        f = figure('Visible','off');
        %f = figure();
        % (Optional) Set the figure size for printing.
        % Here we set PaperUnits to inches and define a PaperPosition vector:
        % [left, bottom, width, height]
        f.PaperUnits = 'inches';      % Units for printing
        f.PaperPosition = [0 0 6 2];    % Size: 6 inch by 4 inch

        % Plot something, e.g., a random line plot.
        hold on
        plot(t*frameToMin, signal, 'black')
        scatter(peakT*frameToMin, peakV, 'red','fill',"v");
        scatter(troughT*frameToMin, troughV, 'blue','fill',"^");
        hold off
        title(sprintf('Pos:%d DropletID:%d', pos,ids(i)));
        xlabel('Time (min)')
        ylabel('FRET/CFP Ratio')

        % Save the figure as a PNG file. The -dpng flag tells MATLAB to create a PNG.
        % The -r300 flag sets the resolution to 300 dpi.
        fn = sprintf('Pos%d_DropletID%d.png', pos,ids(i));
        print(fullfile(savePath, fn), '-dpng', '-r300');

        % Close the figure to free up resources.
        close(f);

        
    end
end

end