function [props, peakMatrix] = plotRaster(dataSet, line_option)

% return values;
peakMatrix = [];
vec_pos = [];
vec_trackid = [];
vec_maxCycles = [];
vec_medianDiameter =[];
%vec_fastestPreiods = [];
%vec_slowestPreiods = [];


count = 0;                              % counting a number of droplets
max_cyc = max(dataSet.cycle.CYCLE_ID);  % max cycle number in this condition
total_peaks = dataSet.cycle;            % peaks information
info = dataSet.info;                    % 
unique_pos = unique(total_peaks.POS_ID);% a list of position of this condition
count_vector =nan(size(unique_pos));    % number of droplet counts for each position


for p = 1:length(unique_pos)
    tmp_cycle_pos = total_peaks(total_peaks.POS_ID == unique_pos(p),:);
    tmp_info_pos = info(info.POS_ID == unique_pos(p),:);
    unique_droplets = unique(tmp_cycle_pos.TRACK_ID);
    for d = 1:length(unique_droplets)
        tmp_droplet = tmp_cycle_pos(tmp_cycle_pos.TRACK_ID == unique_droplets(d),:);
        tmp_info_droplet = tmp_info_pos(tmp_info_pos.TRACK_ID == unique_droplets(d),:);
        tmp_peaks = zeros(max_cyc+1,1)*NaN;
        allpeaks = unique([tmp_droplet.START_FRAME; tmp_droplet.END_FRAME]);
        tmp_peaks(1:height(allpeaks)) = allpeaks;
        peakMatrix = [peakMatrix; tmp_peaks'];
        vec_pos = [vec_pos, unique_pos(p)];
        vec_trackid = [vec_trackid, unique_droplets(d)];
        vec_maxCycles = [vec_maxCycles, tmp_info_droplet.CYCLE_NUMBER+1];
        vec_medianDiameter = [vec_medianDiameter, tmp_info_droplet.ORIGINAL_MED_DIAMETER];
    end
    count = count + length(unique_droplets);
    count_vector(p) = count;
end
props.POS_ID = vec_pos';
props.TRACK_ID = vec_trackid';
props.maxCycles = vec_maxCycles';
props.medianDiameter = vec_medianDiameter';
%{
% all 
if strcmp(mode,'all')
    for p = 1:length(unique_pos)
        tmp_pos = total_peaks(total_peaks.POS_ID == unique_pos(p),:);
        unique_droplets = unique(tmp_pos.TRACK_ID);
        for d = 1:length(unique_droplets)
            tmp_droplet = tmp_pos(tmp_pos.TRACK_ID == unique_droplets(d),:);
            tmp_peaks = zeros(max_cyc+1,1)*NaN;
            allpeaks = unique([tmp_droplet.START_FRAME; tmp_droplet.END_FRAME]);
            tmp_peaks(1:height(allpeaks)) = allpeaks;
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
            tmp_peaks = zeros(max_cyc+1,1)*NaN;
            allpeaks = unique([tmp_droplet.START_FRAME; tmp_droplet.END_FRAME]);
            tmp_peaks(1:height(allpeaks)) = allpeaks;
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
            tmp_peaks = zeros(max_cyc+1,1)*NaN;
            allpeaks = unique([tmp_droplet.START_FRAME; tmp_droplet.END_FRAME]);
            tmp_peaks(1:height(allpeaks)) = allpeaks;
            peakMatrix = [peakMatrix; tmp_peaks'];
        end
        count = count + length(unique_droplets);
        count_vector(p) = count;
    end
else
    disp('Invalid option: Choose all/nuc/cyt');
end
%}

%figure(); 
hold on
for c = 1:size(peakMatrix,2)
    tmp_c = peakMatrix(:,c);
    scatter(tmp_c*dataSet.FrameToMin,1:size(peakMatrix,1),10,'o','filled');
end
ylim([1,size(peakMatrix,1)])
if line_option
    yline(count_vector, '--');
end
xlabel('Time (min)');
hold off;

end