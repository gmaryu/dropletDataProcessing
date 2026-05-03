function tm = getDropletNpixels(db, dropletID, posId, tm)

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
    %fprintf('    Intensity Quantification: ')
    %% test params
    % db = database{2};
    % posId = 20;
    % dropletID = 3;
    % tm = data_SP_control.timeSeries(data_SP_control.timeSeries.POS_ID==posId & data_SP_control.timeSeries.TRACK_ID == dropletID,:);

    %% fileIO
    % Construct file names (using your naming convention).
    croppedImages = fullfile(db.croppedImages,sprintf("droplet_%03d", dropletID));

    %% Prepare cropped images stack
    % --- Load fluorescent images and label ---
    fs = dir(fullfile(croppedImages,'*.tif'));
    names = {fs.name};

    % --- creat H x W x T array for each channel ---
    f = fs(contains(names,'label'));

    % sort by name
    [~,idx] = sort({f.name});
    f = f(idx);
    T = numel(f);

    Hmax = 0;
    Wmax = 0;
    for i = 1:numel(f)
        fname = fullfile(f(i).folder, f(i).name);

        info = imfinfo(fname);
        H = info(1).Height;
        W = info(1).Width;

        Hmax = max(Hmax, H);
        Wmax = max(Wmax, W);
    end
    stack = zeros(Hmax,Wmax,T);
    for t = 1:T
        img = imread(fullfile(f(t).folder, f(t).name));
        h = size(img,1); w = size(img,2);
        r0 = floor((Hmax - h)/2);
        c0 = floor((Wmax - w)/2);
        stack((1:h)+r0, (1:w)+c0, t) = img;
    end

    %% obtain droplet area in pixel
    tm.NPIXEL_DROPLET = squeeze(sum(sum(logical(stack), 1), 2));
    
end
