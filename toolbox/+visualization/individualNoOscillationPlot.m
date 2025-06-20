function individualNoOscillationPlot(oscillationData, frameToMin, savePath)

%%
if ~exist(savePath, 'dir')
    mkdir(savePath);
end

positions = unique(oscillationData.noOscillation.POS_ID)';
ts = oscillationData.noOscillation;

for p = 1:length(positions)
    pos = positions(p);
    
    tmpData = ts(ts.POS_ID == pos,:);
    
    ids = unique(tmpData.TRACK_ID);
    for i = 1:length(ids)
        % temporal signal
        t = tmpData(tmpData.TRACK_ID==ids(i),:).POSITION_T;
        signal = tmpData(tmpData.TRACK_ID==ids(i),:).MAIN_SIGNAL;

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
        plot(t*frameToMin, signal, 'blue')
        
        hold off
        title(sprintf('Pos:%d DropletID:%d', pos,ids(i)));
        xlabel('Time (min)')
        ylabel('FRET/CFP Ratio')

        % Save the figure as a PNG file. The -dpng flag tells MATLAB to create a PNG.
        % The -r300 flag sets the resolution to 300 dpi.
        savePath_p = fullfile(savePath,sprintf('Pos%d',pos));
        if ~exist(savePath_p, 'dir')
            mkdir(savePath_p);
        end
        fn = sprintf('Pos%d_DropletID%03d.png', pos,ids(i));
        %fn = sprintf('Pos%d_DropletID%d.png', pos,ids(i));
        print(fullfile(savePath_p, fn), '-dpng', '-r300');

        % Close the figure to free up resources.
        close(f);

        
    end
end

end