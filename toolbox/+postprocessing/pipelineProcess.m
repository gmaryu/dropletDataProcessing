function data = pipelineProcess(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                                 spermCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
                                 automaticSpermCount, hoechstoffset, FRETNumerator, FRETDenominator)
% pipelineProcess Quantifies fluorescent signals and updates the tracking database.
%
%   data = pipelineProcess(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
%             spermCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, automaticSpermCount, ...
%             hoechstoffset, FRETNumerator, FRETDenominator)
%
% This function processes each position (provided in the cell array "database") for which 
% the position ID is in totalPositions. For each such position, it:
%   - Calls analyzeTrackMate2 to process the tracking data.
%   - Initializes sperm count information if needed.
%   - For each droplet (unique TRACK_ID with peaks), it appends measurements such as median diameter,
%     nuclear data (if spermCondition is true), and cycle metrics. It also applies logic for 
%     automatic sperm count based on nuclear/DNA ratios.
%
% The results for each position are stored back into the database structure and then merged
% across positions into the output structure "data", which also carries conversion factors.
%
% Inputs:
%   database            - Cell array of database structures.
%   totalPositions      - Array of valid position IDs.
%   frameToMin          - Scalar conversion factor from frame number to minutes.
%   pixelToUm           - Scalar conversion factor from pixels to micrometers.
%   initialPeakTimeBound- Scalar bound (in minutes) beyond which oscillations are ignored.
%   forceIgnore         - Table with fields PosID and DropID for droplets to ignore.
%   spermCondition      - Logical; if true, process nuclear/DNA data.
%   nucChannel          - String; channel name used for nuclear staining.
%   dnaChannel          - String; channel name used for DNA staining.
%   overwriteNucMask    - Logical; if true, recalculate nuclear mask.
%   overwriteDNAInfo    - Logical; if true, recalculate DNA quantification.
%   automaticSpermCount - Logical; if true, perform automatic sperm count.
%   hoechstoffset       - Logical; toggle for Hoechst offset correction.
%   FRETNumerator       - String; field name for FRET numerator.
%   FRETDenominator     - String; field name for FRET denominator.
%
% Output:
%   data - Structure with merged tables:
%          data.info       : Summary table with droplet information.
%          data.timeSeries : Appended time series data.
%          data.cycle      : Cycle information.
%          Also includes conversion factors.
%
% Example:
%   data = postprocessing.pipelineProcess(database, totalPositions, 6, 2.649, 100, forceIgnore, true, ...
%             "CFP", "DAPI", true, true, true, true, "MEAN_INTENSITY_CH5", "MEAN_INTENSITY_CH3");

    arguments
        database {iscell(database)}
        totalPositions (:,1) double
        frameToMin (1,1) double {mustBePositive}
        pixelToUm (1,1) double {mustBePositive}
        initialPeakTimeBound (1,1) double {mustBePositive}
        forceIgnore table
        spermCondition logical
        nucChannel (1,1) string
        dnaChannel (1,1) string
        overwriteNucMask logical
        overwriteDNAInfo logical
        automaticSpermCount logical
        hoechstoffset logical
        FRETNumerator (1,1) string
        FRETDenominator (1,1) string
    end

    % Process each position in the database.
    for posIdx = 1:length(database)
        db = database{posIdx};
        if ismember(db.posId, totalPositions)
            % Analyze tracking data and quantify fluorescent signals.
            [trackMate, trackPeaks, trackNoPeaks] = postprocessing.analyzeTrackMate2(db, FRETNumerator, FRETDenominator, frameToMin, forceIgnore);
            
            % Initialize sperm count CSV if it does not exist.
            if ~isfile(db.spermCountCsv)
                spc_table = table([], [], 'VariableNames', {'DropID','Count'});
                writetable(spc_table, db.spermCountCsv);
            end
            spermRef = readtable(db.spermCountCsv);
            
            trackMateAppended = table();
            trackPeaksAppended = table();
            idswPeaks = unique(trackPeaks.TRACK_ID);
            spermCount = [];
            medianDiameter = [];
            ignoredDroplet = [];
            
            % Get droplets to ignore for this position.
            fi = forceIgnore(forceIgnore.PosID == db.posId, :).DropID;
            
            for dropletIdx = 1:length(idswPeaks)
                dropletID = idswPeaks(dropletIdx);
                tm = trackMate(trackMate.TRACK_ID == dropletID, :);
                tp = trackPeaks(trackPeaks.TRACK_ID == dropletID, :);
                
                % Retrieve sperm count for the droplet.
                spermData = spermRef(spermRef.DropID == dropletID, :).Count;
                if isempty(spermData)
                    spermCount(end+1,1) = nan;
                else
                    spermCount(end+1,1) = spermData;
                end
                
                medianDiameter(end+1,1) = median(tm.RADIUS * 2);
                fprintf(" - Droplet %d of Pos %d", dropletID, db.posId);
                
                if ismember(dropletID, fi)
                    ignoredDroplet(end+1,1) = 1;
                    fprintf(" - short tracking frames\n");
                    continue;
                end
                
                if spermCondition
                    nuclearMaskFile = sprintf("%s/nuclear_%03d.mat", db.croppedImages, dropletID);
                    spermMaskFile   = sprintf("%s/dna_%03d.mat", db.croppedImages, dropletID);
                    
                    % Update nuclear and DNA information.
                    postprocessing.nuclearQuantification(database, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo);
                    
                    try
                        nuclearArea = load(nuclearMaskFile).nuclearArea;
                        idxToFrameNuc = load(nuclearMaskFile).idxToFrame;
                        fprintf(" - Nuclear mask obtained");
                        
                        hoechstSum = load(spermMaskFile).hoechstsum;
                        hoechstNPixels = load(spermMaskFile).npts;
                        idxToFrameDNA = load(spermMaskFile).idxToFrame;
                        
                        if ~hoechstoffset
                            smoothbg = load(spermMaskFile).smoothbg;
                            hoechstSum = hoechstSum + smoothbg .* hoechstNPixels;
                        end
                        
                        fprintf(" - Hoechst intensity quantified");
                        assert(all(tm.FRAME == idxToFrameNuc'));
                        assert(all(tm.FRAME == idxToFrameDNA'));
                        
                        tm.NPIXEL_NUC = nuclearArea';
                        tm.NPIXEL_DNA = hoechstNPixels';
                        tm.SUMINTENSITY_DNA = hoechstSum';
                    catch
                        fprintf(" - Ignored. .mat file not found.\n");
                        ignoredDroplet(end+1,1) = 1;
                        continue;
                    end
                else
                    fprintf(" - cytoplasm only -");
                end
                
                if tp.START_FRAME(1) * frameToMin > initialPeakTimeBound
                    fprintf(" - Ignored. Very late oscillation.\n");
                    ignoredDroplet(end+1,1) = 1;
                    continue;
                end
                
                % Process cycle information.
                interphaseStartFrame = [];
                interphaseEndFrame = [];
                CytoNPixelsMedian = [];
                NucNPixelsQ90 = [];
                DNASumIntQ90 = [];
                DNANPixelsQ90 = [];
                DNASumIntIncRate = [];
                
                for cycleIdx = 1:height(tp)
                    startidx = tp.START_INDEX(cycleIdx);
                    endidx = tp.END_INDEX(cycleIdx);
                    cycleData = tm(startidx:endidx-1, :);
                    
                    if spermCondition
                        intstart = find(cycleData.NPIXEL_NUC > 0, 1);
                        intend = length(cycleData.NPIXEL_NUC) - find(flipud(cycleData.NPIXEL_NUC > 0), 1) + 1;
                    else
                        intstart = [];
                    end
                    
                    CytoNPixelsMedian(end+1,1) = median(cycleData.AREA / pixelToUm^2);
                    t = frameToMin * cycleData.FRAME;
                    if spermCondition
                        hoechst = cycleData.SUMINTENSITY_DNA;
                        t = t(~isnan(hoechst));
                        hoechst = hoechst(~isnan(hoechst));
                    else
                        hoechst = zeros(size(cycleData.FRAME));
                    end
                    
                    if numel(t) > 5
                        p = polyfit(t, hoechst, 1);
                        if p(1) > 0
                            DNASumIntIncRate(end+1,1) = log10(p(1));
                        else
                            DNASumIntIncRate(end+1,1) = nan;
                        end
                    else
                        DNASumIntIncRate(end+1,1) = nan;
                    end
                    
                    if isempty(intstart)
                        interphaseStartFrame(end+1,1) = nan;
                        interphaseEndFrame(end+1,1) = nan;
                        NucNPixelsQ90(end+1,1) = nan;
                        DNASumIntQ90(end+1,1) = nan;
                        DNANPixelsQ90(end+1,1) = nan;
                    else
                        interphaseStartFrame(end+1,1) = cycleData.FRAME(intstart);
                        intData = cycleData(cycleData.FRAME >= interphaseStartFrame(end), :);
                        NucNPixelsQ90(end+1,1) = quantile(intData.NPIXEL_NUC, 0.9);
                        DNASumIntQ90(end+1,1) = quantile(intData.SUMINTENSITY_DNA, 0.9);
                        DNANPixelsQ90(end+1,1) = quantile(intData.NPIXEL_DNA, 0.9);
                        if intend > intstart
                            interphaseEndFrame(end+1,1) = cycleData.FRAME(intend);
                        else
                            interphaseEndFrame(end+1,1) = nan;
                        end
                    end
                end
                
                tp.INTERPHASE_START_FRAME = interphaseStartFrame;
                tp.INTERPHASE_END_FRAME = interphaseEndFrame;
                tp.AREA_NPIXELS_MEDIAN = CytoNPixelsMedian;
                tp.NUC_NPIXELS_Q90 = NucNPixelsQ90;
                tp.DNA_SUMINT_Q90 = DNASumIntQ90;
                tp.DNA_NPIXELS_Q90 = DNANPixelsQ90;
                tp.DNA_INC_RATE_COEFF = DNASumIntIncRate;
                
                % Automatic sperm count logic.
                if automaticSpermCount
                    mn = detectMultiNuclei(nuclearMask);
                    if mn > 1
                        spermCount(end) = nan;
                        fprintf(" - Multiple nuclei detected\n");
                    else
                        if height(tp) >= 3
                            if sum((tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN).^(3/2) > 0.0001) > height(tp) - 2
                                if max((tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN).^(3/2)) > 0.05 && ...
                                   max((tp.DNA_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN).^(3/2)) > 0.005
                                    spermCount(end) = 1;
                                    fprintf(" - Single nucleus\n");
                                else
                                    spermCount(end) = nan;
                                    fprintf(" - Fail type 1\n");
                                end
                            elseif sum((tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN).^(3/2) > 0.001) < 2
                                spermCount(end) = 0;
                                fprintf(" - No nucleus\n");
                            else
                                spermCount(end) = nan;
                                fprintf(" - Fail type 2\n");
                            end
                        else
                            spermCount(end) = nan;
                            fprintf(" - Number of cycles too small\n");
                        end
                    end
                end
                
                trackMateAppended = [trackMateAppended; tm];
                trackPeaksAppended = [trackPeaksAppended; tp];
                ignoredDroplet(end+1,1) = 0;
                fprintf(" - \n");
            end
            
            dbInfo = array2table([idswPeaks, ignoredDroplet, spermCount, medianDiameter], ...
                        'VariableNames', {'TRACK_ID', 'IGNORED', 'SPERM_COUNT', 'MEDIAN_DIAMETER'});
            database{posIdx}.info = [array2table(db.posId * ones(height(dbInfo),1), 'VariableNames', {'POS_ID'}), dbInfo];
            database{posIdx}.timeSeries = trackMateAppended;
            database{posIdx}.cycle = trackPeaksAppended;
        end
    end

    % Merge database across positions.
    info = table();
    timeSeries = table();
    cycle = table();
    
    for posIdx = 1:length(database)
        db = database{posIdx};
        if ismember(db.posId, totalPositions)
            posInfo = [array2table(db.posId * ones(height(db.info),1), 'VariableNames', {'POS_ID'}), db.info];
            info = [info; posInfo];
            timeSeries = [timeSeries; [array2table(db.posId * ones(height(db.timeSeries),1), 'VariableNames', {'POS_ID'}), db.timeSeries]];
            cycle = [cycle; [array2table(db.posId * ones(height(db.cycle),1), 'VariableNames', {'POS_ID'}), db.cycle]];
        end
    end
    
    data.info = info;
    data.timeSeries = timeSeries;
    data.cycle = cycle;
    data.FrameToMin = frameToMin;
    data.PixelToUm = pixelToUm;
    data.InitialPeakTimeBound = initialPeakTimeBound;
end
