function createOverlayVideo(result_str, numFrames, Position, CFPChannel, DAPIChannel, maskDir, alphaValue, outputDir)
% createOverlayVideo  Build side‐by‐side overlay movies of CFP and DAPI plus masks.
%
%   createOverlayVideo(sample_ids, numFrames, Position, CFPChannel, DAPIChannel, maskDir, alphaValue, outputDir)
%
% Inputs:
%   sample_ids   – vector of droplet IDs (e.g. [0,1,2,…])
%   numFrames    – number of frames to include per droplet
%   Position     – string prefix before channel in filenames (e.g. 'Pos7')
%   CFPChannel   – string name of CFP channel (e.g. 'CFP')
%   DAPIChannel  – string name of DAPI channel (e.g. 'DAPI')
%   maskDir      – folder containing subfolders 'droplet_###' and mask .mat files
%   alphaValue   – scalar in [0,1] for mask transparency
%   outputDir    – folder in which to save the resulting MP4 videos
%
% For each sample_id, this function:
%   • Loads the per–frame binary masks from 'nuclear_###.mat' (CFP mask)
%     and 'dna_###.mat' (DAPI mask) in maskDir.
%   • Reads up to numFrames TIFFs named Position_CFP_*.tif and Position_DAPI_*.tif
%     from the subfolder 'droplet_###'.
%   • For each timepoint, builds a 1×2 montage:
%       – Left: CFP image in grayscale with blue mask overlay.
%       – Right: DAPI image in grayscale with yellow mask overlay.
%     A title (e.g. "Droplet 001 – Frame 15") is placed above.
%   • Captures each figure as a frame and writes out an MP4 video.
%
% Requires Image Processing Toolbox.
%
% Example:
%   createOverlayVideo(0:5, 50, "Pos7", "CFP", "DAPI", "./masks", 0.3, "./videos")

    arguments
        result_str       (1,1) struct
        numFrames        (1,1) double {mustBePositive}
        Position         (1,1) string
        CFPChannel       (1,1) string
        DAPIChannel      (1,1) string
        maskDir          (1,1) string
        alphaValue       (1,1) double %{mustBeGreaterOrEqual(alphaValue,0), mustBeLessOrEqual(alphaValue,1)}
        outputDir        (1,1) string
    end

    % Ensure output directory exists
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    n = sscanf(Position, 'Pos%d'); 
    sample_ids = result_str.info.TRACK_ID(result_str.info.NUCLEI_COUNT > 0 & result_str.info.POS_ID == n);
    % Loop through each droplet ID
    for id = sample_ids'
        % Subfolder with images
        folderName = fullfile(maskDir, sprintf('droplet_%03d', id));
        if ~isfolder(folderName)
            warning('Folder not found: %s', folderName);
            continue;
        end

        % Load CFP mask
        cfpMaskFile = fullfile(maskDir, sprintf('nuclear_%03d.mat', id));
        if ~isfile(cfpMaskFile)
            warning('CFP mask missing: %s', cfpMaskFile);
            continue;
        end
        cfpData = load(cfpMaskFile);  % expects variable 'nuclearMask'
        cfpMask = cfpData.nuclearMask; % size [H W T]

        % Load DAPI mask
        dapiMaskFile = fullfile(maskDir, sprintf('dna_%03d.mat', id));
        if ~isfile(dapiMaskFile)
            warning('DAPI mask missing: %s', dapiMaskFile);
            continue;
        end
        dapiData = load(dapiMaskFile);  % expects variable 'dnaMask'
        dapiMask = dapiData.dnaMask;     % size [H W T]

        % Collect CFP and DAPI image files
        cfpPattern  = sprintf('%s_%s_*.tif', Position, CFPChannel);
        dapiPattern = sprintf('%s_%s_*.tif', Position, DAPIChannel);
        cfpFiles  = dir(fullfile(folderName, cfpPattern));
        dapiFiles = dir(fullfile(folderName, dapiPattern));
        if isempty(cfpFiles) || isempty(dapiFiles)
            warning('Images missing in %s', folderName);
            continue;
        end

        % Sort them by numeric index
        cfpFiles  = sortByFrameIndex(cfpFiles);
        dapiFiles = sortByFrameIndex(dapiFiles);

        % Cap to numFrames
        nF = min([numel(cfpFiles), numel(dapiFiles), numFrames]);

        % Prepare video writer
        vidName = fullfile(outputDir, sprintf('overlay_%03d.mp4', id));
        vw = VideoWriter(vidName, 'MPEG-4');
        vw.FrameRate = 10;
        open(vw);

        % Create an invisible figure
        %fig = figure('Visible','off','Position',[100 100 800 400]);
        fig = figure('Position',[100 100 800 400]);
        sg = sgtitle(sprintf('Droplet %03d', id));
        sg.FontSize = 14;

        for t = 1:nF
            % Read grayscale images
            Icfp  = imread(fullfile(folderName, cfpFiles(t).name));
            Idapi = imread(fullfile(folderName, dapiFiles(t).name));

            % Normalize to [0 1]
            Icfp  = mat2gray(Icfp);
            Idapi = mat2gray(Idapi);

            % Masks for this frame
            Mcfp  = cfpMask(:,:,min(t,size(cfpMask,3)));
            Mdapi = dapiMask(:,:,min(t,size(dapiMask,3)));

            % Left panel: CFP + blue mask
            ax1 = subplot(1,2,1);
            imshow(Icfp,'InitialMagnification','fit'); hold on;
            h = imshow(cat(3, zeros(size(Mcfp)), zeros(size(Mcfp)), Mcfp));
            set(h,'AlphaData', alphaValue);
            title(sprintf('CFP – Frame %d', t));

            % Right panel: DAPI + yellow mask
            ax2 = subplot(1,2,2);
            imshow(Idapi,'InitialMagnification','fit'); hold on;
            h2 = imshow(cat(3, Mdapi, Mdapi, zeros(size(Mdapi))));
            set(h2,'AlphaData', alphaValue);
            title(sprintf('DAPI – Frame %d', t));

            drawnow;

            % Capture and write frame
            frame = getframe(fig);
            writeVideo(vw, frame);

            cla(ax1); cla(ax2);
        end

        close(vw);
        close(fig);
        fprintf('Saved video: %s\n', vidName);
    end
end

%% Helper function to sort files by the first numeric substring
function sortedFiles = sortByFrameIndex(fileStruct)
    nums = zeros(numel(fileStruct),1);
    for k = 1:numel(fileStruct)
        tokens = regexp(fileStruct(k).name, '(\d+)', 'tokens');
        nums(k) = str2double(tokens{1}{1});
    end
    [~, idx] = sort(nums);
    sortedFiles = fileStruct(idx);
end
