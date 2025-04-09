function database = genDatabase(exportsPath)
% genDatabase  Generates a database structure from exported CSV file names.
%
%   database = genDatabase(exportsPath)
%
% This function searches the specified exportsPath for files ending with 
% "_segmented_tracks.csv". For each matching file, it constructs a structure
% with the following fields:
%   - posId              : Numeric identifier for the position.
%   - trackMateSpotsCsv  : Full file path to the corresponding "PosX_segmented_spots.csv".
%   - trackMateTracksCsv : Full file path to the corresponding "PosX_segmented_tracks.csv".
%   - croppedImages      : Full folder path for the corresponding cropped images.
%   - spermCountCsv      : Full file path to the corresponding "spermcount_PosX.csv".
%   - forceIgnoreCsv     : (Optional) Full file path to "force_ignore.csv", if present.
%
% Input:
%   exportsPath - (1x1 string) Path to the folder containing exported CSV files.
%
% Output:
%   database - (1xN cell array) Each cell contains a structure for each position.
%
% Example:
%   db = genDatabase("exports/20250328_Nocodazole");

    arguments
        exportsPath (1,1) string
    end

    % Search for all files ending with '_segmented_tracks.csv' in the exportsPath.
    fileListStruct = dir(fullfile(exportsPath, '*_segmented_tracks.csv'));
    fileNames = {fileListStruct.name};

    % Initialize the database cell array.
    database = cell(1, numel(fileNames));

    for i = 1:numel(fileNames)
        name = fileNames{i};
        % Use regular expression to extract the position number from the filename "Pos<number>_segmented_tracks.csv".
        tokens = regexp(name, '^Pos(\d+)_segmented_tracks\.csv$', 'tokens');

        if ~isempty(tokens)
            posId = str2double(tokens{1}{1});
            posName = sprintf("Pos%d", posId);

            % Build the struct entry with the appropriate file paths.
            entry = struct();
            entry.posId = posId;
            entry.trackMateSpotsCsv  = fullfile(exportsPath, sprintf('%s_segmented_spots.csv', posName));
            entry.trackMateTracksCsv = fullfile(exportsPath, sprintf('%s_segmented_tracks.csv', posName));
            entry.croppedImages      = fullfile(exportsPath, sprintf('cropped_%s', posName));
            entry.spermCountCsv      = fullfile(exportsPath, sprintf('spermcount_%s.csv', posName));
            entry.forceIgnoreCsv     = fullfile(exportsPath, 'force_ignore.csv');

            database{i} = entry;
        else
            warning('Filename "%s" did not match the expected pattern.', name);
        end
    end

end
