function peakMatrix = plotRaster(dataSet, mode)

peakMatrix = [];


count = 0;
max_cyc = max(dataSet.cycle.CYCLE_ID);
total_peaks = dataSet.cycle;
unique_pos = unique(total_peaks.POS_ID);

count_vector =nan(size(unique_pos));

% all 
if strcmp(mode,'all')
    for p = 1:length(unique_pos)
        tmp_pos = total_peaks(total_peaks.POS_ID == unique_pos(p),:);
        unique_droplets = unique(tmp_pos.TRACK_ID);
        for d = 1:length(unique_droplets)
            tmp_droplet = tmp_pos(tmp_pos.TRACK_ID == unique_droplets(d),:);
            tmp_peaks = zeros(max_cyc,1)*NaN;
            tmp_peaks(1:height(tmp_droplet)) = tmp_droplet.START_FRAME;
            peakMatrix = [peakMatrix; tmp_peaks'];
        end
        count = count + length(unique_droplets);
        count_vector(p) = count;
    end
elseif strcmp(mode,'nuc')
    for p = 1:length(unique_pos)
        tmp_pos = total_peaks(total_peaks.POS_ID == unique_pos(p),:);
        unique_droplets = dataSet.info.TRACK_ID(bitand(dataSet.info.POS_ID == unique_pos(p), dataSet.info.SPERM_COUNT > 0));
        for d = 1:length(unique_droplets)
            tmp_droplet = tmp_pos(tmp_pos.TRACK_ID == unique_droplets(d),:);
            tmp_peaks = zeros(max_cyc,1)*NaN;
            tmp_peaks(1:height(tmp_droplet)) = tmp_droplet.START_FRAME;
            peakMatrix = [peakMatrix; tmp_peaks'];
        end
        count = count + length(unique_droplets);
        count_vector(p) = count;
    end
elseif strcmp(mode,'cyt')
    for p = 1:length(unique_pos)
        tmp_pos = total_peaks(total_peaks.POS_ID == unique_pos(p),:);
        dataSet.info.SPERM_COUNT(isnan(dataSet.info.SPERM_COUNT)) = 0;
        unique_droplets = dataSet.info.TRACK_ID(bitand(dataSet.info.POS_ID == unique_pos(p), dataSet.info.SPERM_COUNT == 0));
        for d = 1:length(unique_droplets)
            tmp_droplet = tmp_pos(tmp_pos.TRACK_ID == unique_droplets(d),:);
            tmp_peaks = zeros(max_cyc,1)*NaN;
            tmp_peaks(1:height(tmp_droplet)) = tmp_droplet.START_FRAME;
            peakMatrix = [peakMatrix; tmp_peaks'];
        end
        count = count + length(unique_droplets);
        count_vector(p) = count;
    end
else
    disp('Invalid option: Choose all/nuc/cyt');
end


%figure(); 
hold on
for c = 1:size(peakMatrix,2)
    tmp_c = peakMatrix(:,c);
    scatter(tmp_c*dataSet.FrameToMin,1:size(peakMatrix,1),10,'o','filled');
end
ylim([1,size(peakMatrix,1)])
yline(count_vector, '--');
xlabel('Time (min)');
hold off;

end