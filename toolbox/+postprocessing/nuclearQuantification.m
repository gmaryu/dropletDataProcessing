function nuclearQuantification(database, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo)
% nuclearQuantification Processes nuclear segmentation and DNA intensity quantification.
%
%   nuclearQuantification(database, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo)
%
% This function loops through each database entry’s cropped images folder and, for each droplet 
% (i.e. each subdirectory), constructs file names for the corresponding nuclear and DNA images.
% If the overwrite flags are true, it calls cropBrightChunk and sumHoechstIntwNucMask to generate 
% the .mat files used for subsequent quantification.
%
% Inputs:
%   database         - Cell array of database structures.
%   nucChannel       - String indicating the nuclear channel name.
%   dnaChannel       - String indicating the DNA channel name.
%   overwriteNucMask - Logical; if true, the nuclear mask will be recalculated.
%   overwriteDNAInfo - Logical; if true, DNA quantification is recalculated.
%
% Example:
%   postprocessing.nuclearQuantification(database, "CFP", "DAPI", true, true);
    
    arguments
        database {iscell(database)}
        nucChannel (1,1) string
        dnaChannel (1,1) string
        overwriteNucMask logical
        overwriteDNAInfo logical
    end

    for i = 1:length(database)
        db = database{i};
        subdirs = dir(db.croppedImages);
        for j = 1:length(subdirs)
            if ~subdirs(j).isdir || any(strcmp(subdirs(j).name, [".", ".."]))
                continue;
            end
            % Build file patterns for nuclear and DNA images.
            nuclearImages = sprintf("%s/%s/Pos%d_%s_???.tif", db.croppedImages, subdirs(j).name, db.posId, nucChannel);
            spermImages   = sprintf("%s/%s/Pos%d_%s_???.tif", db.croppedImages, subdirs(j).name, db.posId, dnaChannel);
            
            nuclearMask = sprintf("%s/%s.mat", db.croppedImages, strrep(subdirs(j).name, "droplet", "nuclear"));
            spermMask   = sprintf("%s/%s.mat", db.croppedImages, strrep(subdirs(j).name, "droplet", "dna"));
            
            if overwriteNucMask
                fprintf(" - Processing %s of Pos %d\n", subdirs(j).name, db.posId);
                try
                    [nuclearArea, idxToFrameNuc] = postprocessing.cropBrightChunk(nuclearImages, nuclearMask);
                    fprintf(" - Nuclear mask obtained.\n");
                catch ME
                    fprintf(" - Ignored due to error: %s\n", ME.message);
                    continue;
                end
            end
            
            if overwriteDNAInfo
                fprintf(" - Processing %s for DNA quantification of Pos %d\n", subdirs(j).name, db.posId);
                [~, ~, ~, ~, idxToFrameDNA] = postprocessing.sumHoechstIntwNucMask(spermImages, nuclearMask, spermMask);
                fprintf(" - Hoechst intensity quantified.\n");
            end
        end
    end
end
