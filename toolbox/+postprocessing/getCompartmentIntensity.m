function [tm, spermCount, nucleiCount] = getCompartmentIntensity(db, dropletID, posId, tm)

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
        return
    end

    if exist(dnaMaskFile, "file")
        dnaData = load(dnaMaskFile);
        dnaMask = dnaData.dnaMask;
        % nucMask = NucDNAMask; % use this line in future. Some NucDNAMask
        % calculated before 2025/09/18 is wrong. Use the line below  to
        % correct.
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
        nucleiCount = NaN;
    end
    
    %% Prepare cropped images stack
    % Load fluorescent images and label
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

        stack = zeros(size(nucData.labelImageCat));
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
    
    % --- Check if FRET pair exist in the color palette ---
    if ismember('CFP', colors)
        if ismember('FRET', colors)
            data.RATIO = double(data.FRET) ./ double(data.CFP);
        elseif ismember('Custom', colors)
            data.RATIO = double(data.Custom) ./ double(data.CFP);
        else
            warning('CFP exits, but no FRET/Custom channel');
        end
    else
        warning('No CFP channel');
    end
  
    % --- create cytoplasmic area mask ---
    cytMask = labelstack & ~nucMask;

    % --- calculate mean intensity in each compartment ---
    colors=fieldnames(data);
    for c = 1:numel(colors)
        ch = colors{c};
        I = data.(ch);
        
        if ~isequal(size(I,1), size(nucMask,1)) || ~isequal(size(I,2), size(nucMask,2)) || ~isequal(size(I,3), size(nucMask,3))
            error('Size mismatch: %s stack and mask3D must have the same H×W×T.', ch);
        end

        num1 = squeeze(sum(sum(I .* cast(nucMask,'like',I), 1), 2));
        den1 = squeeze(sum(sum(nucMask, 1), 2));                    
        den1(den1==0) = NaN;                                     
        y1 = num1 ./ den1;   
        varName1 = matlab.lang.makeValidName(['NucMean_' ch]); 


        num2 = squeeze(sum(sum(I .* cast(cytMask,'like',I), 1), 2)); 
        den2 = squeeze(sum(sum(cytMask, 1), 2));                 
        den2(den2==0) = NaN;                                     
        y2 = num2 ./ den2;   
        varName2 = matlab.lang.makeValidName(['CytMean_' ch]); 

        tm.(varName1) = y1;
        tm.(varName2) = y2;
    end

    % --- DNA related information (previously in getDNAData)
    if exist('dnaMask','var') && ~isempty(dnaMask)
        nucHoechstInt = double(data.DAPI) .* nucMask;
        dnaHoechstInt = double(data.DAPI) .* dnaMask;
        nucHoechstInt_mod = double(data.DAPI) .* nucMask;

        tm.SUM_SPERM_HOECHST_INT = squeeze(sum(sum(dnaHoechstInt, 1), 2));
        tm.SUM_NUCLEUS_HOECHST_INT = squeeze(sum(sum(nucHoechstInt, 1), 2));
        tm.SUM_NUCLEUS_HORCHST_INT_MOD = squeeze(sum(sum(nucHoechstInt_mod, 1), 2));
        tm.NPIXEL_DNA = squeeze(sum(sum(dnaData.dnaMask, 1), 2));
        tm.NPIXEL_NUC_MOD = squeeze(sum(sum(nucMask, 1), 2));

        % Judge whether there are positive hoechst signal pixels
        if max(dnaData.hoechstArea) ~= 0
            % Case for spermCount = 1;
            fprintf('- DNA object detected \n')
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
        end
    end

    %% obtain droplet area in pixel
    tm.AREA_NPIXEL = squeeze(sum(labelstack,[1 2]));
end
