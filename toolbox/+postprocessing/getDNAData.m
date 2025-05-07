function [tm, tp, spermCount, nucleiCount] = getDNAData(db, dropletID, posId, tm, tp, spermCount, hoechstoffset)

% getNuclearData Load nuclear and DNA quantification data for a droplet.
%
%   nuclearData = getNuclearData(db, posId, dropletID, nucChannel, dnaChannel, automaticSpermCount, hoechstoffset)
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
        db
        dropletID (1,1) double
        posId   (1,1) double
        tm  table
        tp table
        spermCount double
        hoechstoffset logical
    end

    %profile on
    
    %% fileIO
    % Construct file names (using your naming convention).
    croppedImages = db.croppedImages;
    maskMatFiles = db.maskMatFiles;
    dnaMaskFile = fullfile(maskMatFiles, sprintf("dna_%03d.mat", dropletID));
    nuclearMaskFile = fullfile(maskMatFiles, sprintf("nuclear_%03d.mat", dropletID));
    hoechstImgFiles = fullfile(croppedImages, sprintf("droplet_%03d/Pos%d_DAPI_???.tif", dropletID, posId));

    % Load mat files
    % All droplets have nuclear_xxx.mat file
    nucData = load(nuclearMaskFile);
    
    % Cytoplasm-only dropelts doesn't have dna_xxx.mat file
    if exist(dnaMaskFile, "file")
        dnaData = load(dnaMaskFile);
    else
        % no nuclear area, therefore no dna positive area
        fprintf(" - No DNA.mat data");
        tm.SPERM_COUNT = NaN*ones(size(tm,1),1);
        tm.SUM_SPERM_HOECHST_INT = NaN*ones(size(tm,1),1);
        tm.SUM_NUCLEUS_HOECHST_INT = NaN*ones(size(tm,1),1);
        tm.NPIXEL_DNA = NaN*ones(size(tm,1),1);
        tm.SUM_NUCLEUS_HORCHST_INT_MOD = NaN*ones(size(tm,1),1);
        tm.NPIXEL_NUC_MOD = NaN*ones(size(tm,1),1);
    
        spermCount = NaN;
        nucleiCount = NaN;
        return
    end
    
    % Load hoechst images
    fs = dir(hoechstImgFiles);
    imgroot = fs(1).folder;
    N = length(fs);
    rawHoechstImages = cell(1, N);
    for i = 1:N
        % raw image stack
        rawHoechstImages{i} = imread(fullfile(imgroot, fs(i).name));
    end
    rawImagesCat = cat(3, rawHoechstImages{:});
    
    % Area specific hoechst 
    nucHoechstInt = double(rawImagesCat).*nucData.nuclearMask;
    dnaHoechstInt = double(rawImagesCat).*dnaData.dnaMask;
    nucHoechstInt_mod = double(rawImagesCat).*dnaData.NucDNAMask;

    tm.SUM_SPERM_HOECHST_INT = squeeze(sum(sum(dnaHoechstInt, 1), 2));
    tm.SUM_NUCLEUS_HOECHST_INT = squeeze(sum(sum(nucHoechstInt, 1), 2));
    tm.SUM_NUCLEUS_HORCHST_INT_MOD = squeeze(sum(sum(nucHoechstInt_mod, 1), 2));
    tm.NPIXEL_DNA = squeeze(sum(sum(dnaData.dnaMask, 1), 2));
    tm.NPIXEL_NUC_MOD = squeeze(sum(sum(dnaData.NucDNAMask, 1), 2));

    %% Sperm Counting 
    % if any(dnaData.dnaMask(:) > 0)
    %     % positive nuclear area and positive dna area
    %     spermCount = 1;
    % else
    %     % positive nuclear area, but no positive dna area
    %     spermCount = 0;
    % end

    % Judge whether there are positive hoechst signal pixels
    if max(dnaData.hoechstArea) ~= 0
        %spermCount = 1;
        fprintf('- DNA object detected')
        [~, spermCountSeries] = postprocessing.detectMultiNuclei(dnaMaskFile);
        spermCount = spermCountSeries(1);% superm number at first time frame when tracking started
        tm.SPERM_COUNT = spermCountSeries';
        nucleiCount = max(tm.NUCLEI_COUNT); % maximum nuclei number in a time series

        % check co-localization
        hasPositive = any(nucData.nuclearMask.*dnaData.dnaMask > 0, 'all');
        if hasPositive == 0
            spermCount = 0;
            nucleiCount = 0;
            tm.SPERM_COUNT = NaN*ones(size(tm,1),1);
            tm.NUCLEI_COUNT = NaN*ones(size(tm,1),1);
            tm.NPIXEL_DNA = NaN*ones(size(tm,1),1);
            tm.NPIXEL_NUC = NaN*ones(size(tm,1),1);
            tm.SUM_NUCLEUS_HORCHST_INT_MOD = NaN*ones(size(tm,1),1);
            tm.NPIXEL_NUC_MOD = NaN*ones(size(tm,1),1);
            fprintf('- No colocalization');
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
        % nuclear positve are, but no DNA positve area. (Nuclear area artifact)
        fprintf('- No DNA object detected')
        tm.SPERM_COUNT = NaN*ones(size(tm,1),1);
        tm.SUM_SPERM_HOECHST_INT = NaN*ones(size(tm,1),1);
        tm.SUM_NUCLEUS_HOECHST_INT = NaN*ones(size(tm,1),1);
        tm.NPIXEL_DNA = NaN*ones(size(tm,1),1);
        tm.SUM_NUCLEUS_HORCHST_INT_MOD = NaN*ones(size(tm,1),1);
        tm.NPIXEL_NUC_MOD = NaN*ones(size(tm,1),1);
        spermCount = NaN;
        nucleiCount = NaN;
    end    
  
    %profile off
end
