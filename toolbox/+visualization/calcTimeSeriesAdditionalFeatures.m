function dataSet = calcTimeSeriesAdditionalFeatures(dataSet)
dnarenormfactor = 1; % 1e7 factor to make numbers not too big

% --- make IGNORE column if necessary ---
if ~ismember("IGNORED", dataSet.info.Properties.VariableNames)
    disp("IGNORED is added to dataSet.Info");
    dataSet.info.IGNORED = zeros(size(dataSet.info,1),1);
end

if ~ismember("IGNORED", dataSet.cycle.Properties.VariableNames)
    disp("IGNORED is added to dataSet.Info");
    dataSet.cycle.IGNORED = zeros(size(dataSet.cycle,1),1);
end

if ~ismember("IGNORED", dataSet.cycle.Properties.VariableNames)
    disp("IGNORED is added to dataSet.Info");
    dataSet.timeSeries.IGNORED = zeros(size(dataSet.timeSeries,1),1);
end

% --- extract variables ---
tm = dataSet.timeSeries;

%% Update existing features
tm.MINUTE = tm.FRAME * dataSet.FrameToMin; % minute

% --- check if a new column has to be calculated / pixel number of droplet----

% if ~ismember('AREA_NPIXEL', string(tm.Properties.VariableNames))
%     disp('New column AREA_NPIXEL added');
%     tm.AREA_NPIXEL = tm.AREA ./ (dataSet.PixelToUm^2);
% end
if ismember("NPIXEL_DROPLET", string(tm.Properties.VariableNames))
    tm.DIAMETER = 2 * sqrt(tm.NPIXEL_DROPLET * dataSet.PixelToUm^2 / pi); 
else
    tm.DIAMETER = 2 * tm.RADIUS;
end
tm.ORIGINAL_MED_DIAMETER = postprocessing.originalDiameter(tm.DIAMETER, tm.TUBE_HEIGHT);
tm.VOLUMEUM3 = (4/3)*pi*(tm.ORIGINAL_MED_DIAMETER./2).^3; %.convertAreaPixelsToVolume(tm.DROPLET_NPIXEL, dataSet.PixelToUm); % um^3

%
% --- correct old "NPIXEL_NUC_MOD" and "SUM_NUCLEUS_HORCHST_INT_MOD"
if ismember('NPIXEL_NUC_MOD', string(tm.Properties.VariableNames)) && all(~isnan(tm.NPIXEL_NUC_MOD))
    % if any(tm.NPIXEL_NUC_MOD ~= max(tm.NPIXEL_NUC, tm.NPIXEL_DNA))
    %     tm.NPIXEL_NUC_MOD = max(tm.NPIXEL_NUC, tm.NPIXEL_DNA);
    %     disp('NPIXEL_NUC_MOD is fixed')
    % end

    tm.NUC_VOLUMEUM3 = visualization.convertAreaPixelsToVolume(tm.NPIXEL_NUC_MOD, dataSet.PixelToUm); % um^3
    tm.NCVR = tm.NUC_VOLUMEUM3 ./ tm.VOLUMEUM3;
    tm.NUC_SURF_AREA = 4.*tm.NPIXEL_NUC_MOD.*dataSet.PixelToUm^2;
elseif ismember('NPIXEL_NUC', string(tm.Properties.VariableNames))
    tm.NUC_VOLUMEUM3 = visualization.convertAreaPixelsToVolume(tm.NPIXEL_NUC, dataSet.PixelToUm); % um^3
    tm.NCVR = tm.NUC_VOLUMEUM3 ./ tm.VOLUMEUM3;
    tm.NUC_SURF_AREA = 4.*tm.NPIXEL_NUC.*dataSet.PixelToUm^2;
end

if ismember('SUM_NUCLEUS_HORCHST_INT_MOD', string(tm.Properties.VariableNames))
    % if any(tm.SUM_NUCLEUS_HORCHST_INT_MOD ~= max(tm.SUM_SPERM_HOECHST_INT, tm.SUM_NUCLEUS_HOECHST_INT))
    %     tm.SUM_NUCLEUS_HORCHST_INT_MOD = max(tm.SUM_SPERM_HOECHST_INT, tm.SUM_NUCLEUS_HOECHST_INT);
    %     disp('SUM_NUCLEUS_HORCHST_INT_MOD is fixed')
    % end

    tm.DNACR_NUC = tm.SUM_NUCLEUS_HORCHST_INT_MOD ./ tm.VOLUMEUM3 / dnarenormfactor; % a.u (sum px) / um^3
    tm.DNACR_DNA = tm.SUM_SPERM_HOECHST_INT ./ tm.VOLUMEUM3 / dnarenormfactor;
end

dataSet.timeSeries = tm;
end