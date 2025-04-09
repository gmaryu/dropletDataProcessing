function nuclearData = getNuclearData(croppedImages, posId, dropletID, nucChannel, dnaChannel, hoechstoffset)
% getNuclearData Load nuclear and DNA quantification data for a droplet.
%
%   nuclearData = getNuclearData(croppedImages, posId, dropletID, nucChannel, dnaChannel, hoechstoffset)
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
        posId (1,1) double
        dropletID (1,1) double
        nucChannel (1,1) string
        dnaChannel (1,1) string
        hoechstoffset logical
    end

    % Construct file names (using your naming convention).
    nuclearMaskFile = fullfile(croppedImages, sprintf("nuclear_%03d.mat", dropletID));
    dnaMaskFile = fullfile(croppedImages, sprintf("dna_%03d.mat", dropletID));
    
    try
        nucData = load(nuclearMaskFile);
        dnaData = load(dnaMaskFile);
    catch ME
        error("Failed to load nuclear or DNA data: %s", ME.message);
    end
    
    nuclearData.nuclearArea = nucData.nuclearArea;
    % Optional: you might want to include idxToFrame if needed.
    % nuclearData.idxToFrameNuc = nucData.idxToFrame;
    
    nuclearData.hoechstSum = dnaData.hoechstsum;
    nuclearData.hoechstNPixels = dnaData.npts;
    
    if ~hoechstoffset
        % If offset correction is disabled, adjust using smoothbg (if available).
        if isfield(dnaData, 'smoothbg')
            nuclearData.hoechstSum = nuclearData.hoechstSum + dnaData.smoothbg .* nuclearData.hoechstNPixels;
        end
    end
end
