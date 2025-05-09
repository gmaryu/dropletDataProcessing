function nuclearSegmentation(db, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo)
% nuclearQuantification  Process nuclear segmentation and DNA intensity quantification.
%
%   nuclearQuantification(database, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo)
%
% This function iterates over each entry in the given database (a cell array of
% structures). For each droplet (each subdirectory within the croppedImages folder),
% the function constructs file patterns for the nuclear and DNA images based on the specified
% channels. If the corresponding mask files are to be updated (based on the overwrite flags), it calls
% the segmentation routines cropBrightChunk and sumHoechstIntwNucMask to generate the .mat files.
%
% Inputs:
%   database         - Cell array of database structures. Each structure should include at least:
%                        .posId: numeric identifier for the position.
%                        .croppedImages: path to the folder containing cropped droplet images.
%   nucChannel       - (1,1) string specifying the nuclear channel (e.g., "CFP").
%   dnaChannel       - (1,1) string specifying the DNA (Hoechst) channel (e.g., "DAPI").
%   overwriteNucMask - Logical flag; if true, the nuclear segmentation (mask) is recalculated.
%   overwriteDNAInfo - Logical flag; if true, the DNA intensity quantification is recalculated.
%
% Example:
%   postprocessing.nuclearQuantification(database, "CFP", "DAPI", true, true);
%
% Note: This function assumes that the naming convention for the droplet directories
%   is such that the folder name includes "droplet" and that the nuclear and DNA files
%   are named by replacing "droplet" with "nuclear" or "dna", respectively.

arguments
    db {iscell(db)}
    nucChannel (1,1) string
    dnaChannel (1,1) string
    overwriteNucMask logical
    overwriteDNAInfo logical
end

% Prepare a folder for output masks
mkdir(db.maskMatFiles);

% Process only subdirectories (skip '.' and '..') List of cropped droplets
subdirs = dir(db.croppedImages);
subdirs = subdirs([subdirs(:).isdir]);
subdirs = subdirs(~ismember({subdirs(:).name},{'.','..'}));

for j = 1:length(subdirs)

    % Construct the file patterns for nuclear and DNA images.
    nuclearImages = fullfile(db.croppedImages, subdirs(j).name, ...
        sprintf("Pos%d_%s_???.tif", db.posId, nucChannel));
    dnaImages     = fullfile(db.croppedImages, subdirs(j).name, ...
        sprintf("Pos%d_%s_???.tif", db.posId, dnaChannel));
    labelImages = fullfile(db.croppedImages, subdirs(j).name, ...
        sprintf("Pos%d_label_???.tif", db.posId));

    % Construct the output file names for the segmentation results.
    nuclearMaskFile = fullfile(db.maskMatFiles, sprintf("%s.mat",strrep(subdirs(j).name, "droplet", "nuclear")));
    dnaMaskFile = fullfile(db.maskMatFiles, sprintf("%s.mat",strrep(subdirs(j).name, "droplet", "dna")));

    % Process nuclear segmentation if requested.
    if overwriteNucMask
        fprintf("Processing nuclear mask for %s of Pos %d...\n", subdirs(j).name, db.posId);
        try
            % Call the cropping/segmentation routine.
            [nuclearArea, idxToFrameNuc] = postprocessing.cropBrightChunk(nuclearImages, labelImages, nuclearMaskFile);
            fprintf("Nuclear mask obtained.\n");
        catch ME
            fprintf("Nuclear segmentation failed for %s: %s\n", subdirs(j).name, ME.message);
            continue;
        end
    
    end

    % Process DNA (Hoechst) intensity if requested.
    if overwriteDNAInfo
        fprintf("Processing DNA segmentation for %s of Pos %d...\n", subdirs(j).name, db.posId);
        try
            [dnaArea, idxToFrameDNA] = postprocessing.cropDNAMask(dnaImages, labelImages, nuclearMaskFile, dnaMaskFile);
            fprintf("DNA quantification completed.\n");
        catch ME
            fprintf("DNA quantification failed for %s: %s\n", subdirs(j).name, ME.message);
            continue;
        end
    end
end



end
