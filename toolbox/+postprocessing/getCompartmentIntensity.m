function [tm, spermCount, nucleiCount] = getCompartmentIntensity(db, dropletID, posId, tm)

% getCompartmentIntensity measure nuclear, cytoplasmic, and whole area fluorescent intensity for a droplet.
%
%   [tm, spermCount, nucleiCount] = getCompartmentIntensity(db, dropletID, posId, tm)
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
    end
    fprintf('    Intensity Quantification: ')
    %% test params
    % db = database{2};
    % posId = 10;
    % dropletID = 33;
    % tm = control_nuc.timeSeries(control_nuc.timeSeries.POS_ID==posId & control_nuc.timeSeries.TRACK_ID == dropletID,:);

    %% fileIO
    % Construct file names (using your naming convention).
    croppedImages = fullfile(db.croppedImages,sprintf("droplet_%03d", dropletID));
    maskMatFiles = db.maskMatFiles;
    dnaMaskFile = fullfile(maskMatFiles, sprintf("dna_%03d.mat", dropletID));
    nuclearMaskFile = fullfile(maskMatFiles, sprintf("nuclear_%03d.mat", dropletID));
    
    %% Load mask files (single droplet mask stack (logical))
    % All droplets have nuclear_xxx.mat file
    if exist(nuclearMaskFile, "file")
        nucData = load(nuclearMaskFile);
        nucMask = nucData.nuclearMask;
    else
        fprintf('-- No Nuclear mask .mat file');
        spermCount = NaN;
        nucleiCount = NaN;
        return
    end

    if exist(dnaMaskFile, "file")
        dnaData = load(dnaMaskFile);
        dnaMask = dnaData.dnaMask;
        % nucMask = NucDNAMask; % use this line in future. Some NucDNAMask
        % calculated before 2025/09/18 is wrong. Use the line below  to
        % correct.
        nucMask_original = nucMask;
        nucMask = nucMask | dnaMask;
    else
        fprintf('-- No DNA mask .mat file');
        tm.SPERM_COUNT = NaN*ones(size(tm,1),1);
        tm.SUM_SPERM_HOECHST_INT = NaN*ones(size(tm,1),1);
        tm.SUM_NUCLEUS_HOECHST_INT = NaN*ones(size(tm,1),1);
        tm.NPIXEL_DNA = NaN*ones(size(tm,1),1);
        tm.SUM_NUCLEUS_HORCHST_INT_MOD = NaN*ones(size(tm,1),1);
        tm.NPIXEL_NUC_MOD = NaN*ones(size(tm,1),1);

        spermCount = NaN;
        nucleiCount = max(tm.NUCLEI_COUNT);
    end
    
    %% Prepare cropped images stack
    % --- Load fluorescent images and label ---
    fs = dir(fullfile(croppedImages,'*.tif'));
    names = {fs.name};

    % --- extract fluorescent channels information ---
    colors = extractBetween(names, "_", "_");  
    colors = colors(~ismissing(colors));
    colors = unique(cellstr(colors));
    if isfield(nucData,'labelImageCat')
        omitList = {'BF','label'};
    else
        omitList = {'BF'};
    end
    colors(ismember(colors, omitList)) = [];

    % --- creat H x W x T array for each channel ---
    data = struct();
    for c = 1:numel(colors)
        f = fs(contains(names,colors{c}));

        % sort by name
        [~,idx] = sort({f.name});
        f = f(idx);
        T = numel(f);

        stack = zeros(size(nucMask));
        Hmax = size(stack,1); Wmax = size(stack,2);
        for t = 1:T
    		img = imread(fullfile(f(t).folder, f(t).name));
    		h = size(img,1); w = size(img,2);
    		r0 = floor((Hmax - h)/2);
            c0 = floor((Wmax - w)/2);
            stack((1:h)+r0, (1:w)+c0, t) = img;
        end

        %{
        stack = zeros(h,w,T,'like',sample);
        for t = 1:T
            stack(:,:,t) = imread(fullfile(f(t).folder, f(t).name));
        end
        %}

        data.(colors{c}) = stack;
    end
  
    if isfield(nucData,'labelImageCat')
        labelstack = logical(nucData.labelImageCat);
    else
        labelstack = logical(data.label);
        data = rmfield(data,'label');
    end
    
    % --- Check if FRET pair exists in the color palette ---
    % Do NOT create pixel-wise ratio image here.
    % Ratio will be calculated later as:
    %   sum(FRET in mask) / sum(CFP in mask)
    hasRatio = false;
    ratioNumeratorChannel = "";

    if isfield(data, 'CFP')
        if isfield(data, 'FRET')
            hasRatio = true;
            ratioNumeratorChannel = "FRET";
        elseif isfield(data, 'Custom')
            hasRatio = true;
            ratioNumeratorChannel = "Custom";
        else
            warning('CFP exists, but no FRET/Custom channel');
        end
    else
        warning('No CFP channel');
    end

    % --- create cytoplasmic area mask ---
    cytMask = labelstack & ~nucMask;

    % --- calculate mean intensity in each compartment ---
    % For ordinary channels:
    %   mean intensity = sum(I in mask) / number of pixels in mask
    %
    % For FRET/CFP ratio:
    %   ratio = sum(FRET in mask) / sum(CFP in mask)
    % instead of mean(FRET ./ CFP).

    colors = fieldnames(data);

    % Ensure masks are binary logical masks
    nucMaskBin   = nucMask > 0;
    cytMaskBin   = cytMask > 0;
    wholeMaskBin = labelstack > 0;

    for c = 1:numel(colors)
        ch = colors{c};
        I = data.(ch);

        if ~isequal(size(I,1), size(nucMaskBin,1)) || ...
                ~isequal(size(I,2), size(nucMaskBin,2)) || ...
                ~isequal(size(I,3), size(nucMaskBin,3))
            error('Size mismatch: %s stack and mask3D must have the same H×W×T.', ch);
        end

        % Nuclear mean intensity
        y1 = maskedMean3D(I, nucMaskBin);
        varName1 = matlab.lang.makeValidName(['NucMean_' ch]);

        % Cytoplasmic mean intensity
        y2 = maskedMean3D(I, cytMaskBin);
        varName2 = matlab.lang.makeValidName(['CytMean_' ch]);

        % Whole droplet / whole compartment mean intensity
        y3 = maskedMean3D(I, wholeMaskBin);
        varName3 = matlab.lang.makeValidName(['WholeMean_' ch]);

        tm.(varName1) = y1;
        tm.(varName2) = y2;
        tm.(varName3) = y3;
    end

    % --- calculate FRET/CFP ratio as sum numerator / sum CFP ---
    if hasRatio
        Fnum = data.(ratioNumeratorChannel);
        CFP  = data.CFP;

        tm.NucMean_RATIO = maskedSumRatio3D(Fnum, CFP, nucMaskBin);
        tm.CytMean_RATIO = maskedSumRatio3D(Fnum, CFP, cytMaskBin);
        tm.WholeMean_RATIO = maskedSumRatio3D(Fnum, CFP, wholeMaskBin);
    end
    
    % --- DNA related information (previously in getDNAData)
    if exist('dnaMask','var') && ~isempty(dnaMask)
        nucHoechstInt = double(data.DAPI) .* nucMask_original;
        dnaHoechstInt = double(data.DAPI) .* dnaMask;
        nucHoechstInt_mod = double(data.DAPI) .* nucMask;

        tm.SUM_SPERM_HOECHST_INT = squeeze(sum(sum(dnaHoechstInt, 1), 2));
        tm.SUM_NUCLEUS_HOECHST_INT = squeeze(sum(sum(nucHoechstInt, 1), 2));
        tm.SUM_NUCLEUS_HORCHST_INT_MOD = squeeze(sum(sum(nucHoechstInt_mod, 1), 2));
        tm.NPIXEL_NUC = squeeze(sum(sum(nucMask_original, 1), 2));
        tm.NPIXEL_DNA = squeeze(sum(sum(dnaMask, 1), 2));
        tm.NPIXEL_NUC_MOD = squeeze(sum(sum(nucMask, 1), 2));

        % Judge whether there are positive hoechst signal pixels
        if max(dnaData.hoechstArea) ~= 0
            % Case for spermCount = 1;
            fprintf('- DNA object detected')
            [~, spermCountSeries] = postprocessing.detectMultiNuclei(dnaMaskFile);
            spermCount = spermCountSeries(1);% superm number at first time frame when tracking started
            tm.SPERM_COUNT = spermCountSeries';
            nucleiCount = max(tm.NUCLEI_COUNT); % maximum nuclei number in a time series

            % check co-localization
            hasPositive = any(nucData.nuclearMask.*dnaData.dnaMask > 0, 'all');
            if hasPositive == 0
                % overwrite to NaN
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
        else
            spermCount = 0;
            nucleiCount = 0;
            tm.SPERM_COUNT = NaN*ones(size(tm,1),1);
            tm.NUCLEI_COUNT = NaN*ones(size(tm,1),1);
            tm.NPIXEL_DNA = NaN*ones(size(tm,1),1);
            tm.NPIXEL_NUC = NaN*ones(size(tm,1),1);
            tm.SUM_NUCLEUS_HORCHST_INT_MOD = NaN*ones(size(tm,1),1);
            tm.NPIXEL_NUC_MOD = NaN*ones(size(tm,1),1);

            fprintf('- No DNA mask');
        end
    end

end

function y = maskedMean3D(I, mask)
% maskedMean3D
% Calculates frame-wise mean intensity inside a binary mask.
%
% I    : H x W x T image stack
% mask : H x W x T logical mask
%
% y(t) = sum(I(:,:,t) within mask(:,:,t)) / number of mask pixels

    mask = mask > 0;

    num = squeeze(sum(sum(double(I) .* double(mask), 1), 2));
    den = squeeze(sum(sum(mask, 1), 2));

    den(den == 0) = NaN;
    y = num ./ den;
end

function r = maskedSumRatio3D(numeratorStack, denominatorStack, mask)
% maskedSumRatio3D
% Calculates frame-wise ratio as:
%
%   sum(numeratorStack within mask) / sum(denominatorStack within mask)
%
% This is recommended for FRET/CFP ratio quantification because it avoids
% averaging pixel-wise ratios.

    mask = mask > 0;

    if ~isequal(size(numeratorStack), size(denominatorStack)) || ...
       ~isequal(size(numeratorStack), size(mask))
        error('Size mismatch: numerator, denominator, and mask must have the same H×W×T size.');
    end

    num = squeeze(sum(sum(double(numeratorStack)   .* double(mask), 1), 2));
    den = squeeze(sum(sum(double(denominatorStack) .* double(mask), 1), 2));

    den(den == 0) = NaN;
    r = num ./ den;
end