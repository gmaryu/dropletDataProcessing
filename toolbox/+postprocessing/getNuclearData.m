function [tm, tp, nucleiCount] = getNuclearData(croppedImages, dropletID, tm, tp, nucleiCount, automaticNucleiCount)
% getNuclearData Load nuclear and DNA quantification data for a droplet.
%
%   nuclearData = getNuclearData(croppedImages, posId, dropletID, nucChannel, dnaChannel, automaticSpermCount, hoechstoffset)
%
% This function constructs the filenames for the nuclear and DNA (Hoechst) data based on
% the croppedImages directory, position, and droplet ID. It then loads the corresponding .mat
% files and, if hoechstoffset is false, adjusts the Hoechst sum accordingly.
%
% Inputs:
%   croppedImages - (1,1) string specifying the directory with cropped droplet images.
%   dropletID     - (1,1) numeric droplet identifier.

%   automaticSpermCount - (1,1) logical flag indicating whether to apply auto nuclei conunting.
%
% Output:
%   nuclearData - A structure with fields:
%       .nuclearArea    - Array with the nuclear area (in pixels) per frame.
%       .hoechstNPixels - Array with the number of Hoechst-positive pixels per frame.
%       .hoechstSum     - Array with the total Hoechst intensity per frame.
%
% Example:
%   nd = getNuclearData("exports/20250328_Nocodazole/cropped_pos0", 0, 5, "CFP", "DAPI", true);

    arguments
        croppedImages (1,1) string
        dropletID (1,1) double
        tm  table
        tp table
        nucleiCount double
        automaticNucleiCount logical
    end

    % Construct file names (using your naming convention).
    nuclearMaskFile = fullfile(croppedImages, sprintf("nuclear_%03d.mat", dropletID));
    
    % Load mat files
    try
        nucData = load(nuclearMaskFile);
    catch ME
        error("Failed to load nuclear data: %s", ME.message);
    end
    
    nuclearData.nuclearArea = nucData.nuclearArea; % Vector of 1 x timepoints
    nuclearData.idxToFrameNuc = nucData.idxToFrame; % Vector of 1 x timepoints
    
    % sanity check, frame consistency
    assert(all(tm.FRAME == nuclearData.idxToFrameNuc'));

    % Append the obtained data to the tracking table.
    tm.NPIXEL_NUC = nuclearData.nuclearArea';

    % Judge whether there are positive nuclear signal pixels
    if max(nuclearData.nuclearArea) ~= 0
        nucleiCount = 1;
        fprintf('- Nuclear object detected')
        if automaticNucleiCount
            [nucleiCount, nucleiCountSeries] = postprocessing.detectMultiNuclei(nuclearMaskFile);
            tm.NUCLEI_COUNT = nucleiCountSeries';
        end
        %{
        if mn > 1
            spermCount = mn;
        else
            if size(tp, 1) >= 3 % more than three peaks
                if sum(power(tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3/2) > 0.0001) > size(tp, 1) - 2
                    
                    % single nuclear exists
                    if max(power(tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3/2)) > 0.05 && max(power(tp.DNA_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3/2)) > 0.005
                        % maximum n/c volume > 0.05 & maximum Hoechst volume > 0.005
                        spermCount = 1;
                        fprintf(" - Single nucleus\n");
                    else
                        spermCount = nan;
                        fprintf(" - Fail type 1 - bright pixel in DAPI but area is not enough large\n");
                    end
                    
                elseif sum(power(tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3/2) > 0.001) < 2
                    
                    spermCount = 0;
                    fprintf(" - No nucleus\n");
                    
                else
                    
                    spermCount = nan;
                    fprintf(" - Fail type 2 - uncategorized error\n");
                    
                end
                
            else
                spermCount = NaN;
                fprintf(" - Number of cycle too small -");
            end
        end
        %}    
    else
        fprintf('- No nuclear object detected')
        tm.NUCLEI_COUNT = NaN*ones(size(tm,1),1);
        nucleiCount = NaN;
    end    
    
end
