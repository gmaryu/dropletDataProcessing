function [tm, tp, spermCount] = getNuclearData(croppedImages, dropletID, tm, tp, spermCount, automaticSpermCount, hoechstoffset)
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
%   posId         - (1,1) numeric position identifier.
%   dropletID     - (1,1) numeric droplet identifier.
%   nucChannel    - (1,1) string for the nuclear channel name.
%   dnaChannel    - (1,1) string for the DNA/Hoechst channel name.
%   automaticSpermCount - (1,1) logical flag indicating whether to apply auto nuclei conunting.
%   hoechstoffset - (1,1) logical flag indicating whether to apply Hoechst offset correction.
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
        spermCount double
        automaticSpermCount logical
        hoechstoffset logical
    end

    % Construct file names (using your naming convention).
    nuclearMaskFile = fullfile(croppedImages, sprintf("nuclear_%03d.mat", dropletID));
    dnaMaskFile = fullfile(croppedImages, sprintf("dna_%03d.mat", dropletID));
    
    % Load mat files
    try
        nucData = load(nuclearMaskFile);
        dnaData = load(dnaMaskFile);
    catch ME
        error("Failed to load nuclear or DNA data: %s", ME.message);
    end
    
    nuclearData.nuclearArea = nucData.nuclearArea;
    nuclearData.idxToFrameNuc = nucData.idxToFrame;
    
    nuclearData.hoechstSum = dnaData.hoechstsum;
    nuclearData.hoechstNPixels = dnaData.npts;
    nuclearData.idxToFrameDNA = dnaData.idxToFrame;

    if ~hoechstoffset
        % If offset correction is disabled, adjust using smoothbg (if available).
        if isfield(dnaData, 'smoothbg')
            nuclearData.hoechstSum = nuclearData.hoechstSum + dnaData.smoothbg .* nuclearData.hoechstNPixels;
        end
    end
    
    % sanity check, frame consistency
    assert(all(tm.FRAME == nuclearData.idxToFrameNuc'));
    assert(all(tm.FRAME == nuclearData.idxToFrameDNA'));

    % Append the obtained data to the tracking table.
    tm.NPIXEL_NUC = nuclearData.nuclearArea';
    tm.NPIXEL_DNA = nuclearData.hoechstNPixels';
    tm.SUMINTENSITY_DNA = nuclearData.hoechstSum';


    % define logic to count sperm dna copies
    if automaticSpermCount
        mn = postprocessing.detectMultiNuclei(nuclearMaskFile);
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
                    %}
                elseif sum(power(tp.NUC_NPIXELS_Q90 ./ tp.AREA_NPIXELS_MEDIAN, 3/2) > 0.001) < 2
                    
                    spermCount = 0;
                    fprintf(" - No nucleus\n");
                    
                else
                    
                    spermCount = nan;
                    fprintf(" - Fail type 2 - uncategorized error\n");
                    
                end
                %}
            else
                spermCount = NaN;
                fprintf(" - Number of cycle too small -");
            end
        end

    end    

end
