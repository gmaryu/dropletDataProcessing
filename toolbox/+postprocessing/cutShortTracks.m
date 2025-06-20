function fi_table = cutShortTracks(database, fi_table, thres)
% cutShortTracks  Append entries for short tracks to the fault indication table.
%
%   fi_table = cutShortTracks(database, fi_table, thres)
%
% This function processes a cell array "database" where each element is a structure 
% containing tracking information. For each position, it reads the tracking data from 
% the provided CSV file and identifies tracks whose length (as indicated by the 
% NUMBER_SPOTS column) is less than a specified threshold fraction (thres) of the 
% maximum track length at that position. For each such short track, an entry is added 
% to the input fault indication table (fi_table) with the position ID, track ID, and a 
% reason ("short trace").
%
% Inputs:
%   database - (cell array) Each cell is a structure with at least the fields:
%              .posId             : Numeric identifier for the position.
%              .trackMateTracksCsv: Full file path to the "segmented_tracks" CSV file.
%
%   fi_table - (table) An existing table with columns "PosID", "DropID", and "Reason"
%              to which new entries will be appended.
%
%   thres    - (scalar double) A fraction (0 < thres < 1) that defines the threshold 
%              for a track to be considered short relative to the longest track.
%
% Output:
%   fi_table - (table) The updated table including entries for all tracks that meet 
%              the "short trace" condition.
%
% Example:
%   fi_table = table('Size',[0 3],'VariableTypes',{'double','double','cell'},...
%                    'VariableNames',{'PosID','DropID','Reason'});
%   database = {... % populate database
%               };
%   fi_table = cutShortTracks(database, fi_table, 0.7);

    arguments
        database {iscell(database)}
        fi_table table
        thres (1,1) double {mustBeGreaterThan(thres, 0), mustBeLessThan(thres, 1)}
    end

    for i = 1:length(database)
        db = database{i};

        % Read the tracking data from CSV
        trackMateTracks = readtable(db.trackMateTracksCsv);

        % Extract tracking lengths from the NUMBER_SPOTS column.
        tracking_lengths = trackMateTracks.NUMBER_SPOTS;

        % Find the maximum track length ignoring NaN values.
        max_len = nanmax(tracking_lengths);

        % Identify the tracks shorter than max_len * thres.
        mask = tracking_lengths < (max_len * thres);

        % If any short tracks exist, record their details.
        if any(mask)
            numShort = sum(mask);
            tmpTable = table(...
                repmat(db.posId, numShort, 1), ...                     % PosID repeated for each short track
                trackMateTracks.TRACK_ID(mask), ...                     % Short track IDs
                repmat({"short trace"}, numShort, 1), ...                 % Reason for flagging
                'VariableNames', {'PosID', 'DropID', 'Reason'});

            % Append the temporary table to the input fi_table.
            fi_table = [fi_table; tmpTable];
        end
    end

end