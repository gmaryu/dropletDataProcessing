function fi_table = radiusDiffDetection(database, fi_table, thres)
% radiusDiffDetection  Detects tracks with significant radius changes.
%
%   fi_table = radiusDiffDetection(database, fi_table, thres)
%
% This function iterates over the entries in "database" (a cell array of structures)
% and processes the corresponding trackMateSpots CSV file for each position. For every 
% unique track (identified by TRACK_ID), it computes the differences between consecutive
% radius measurements. If the maximum absolute difference exceeds a fraction (thres) of the 
% initial radius, the track is flagged as having a "radius change". The fault information 
% is then appended to the fi_table.
%
% Inputs:
%   database - (cell array) Each cell is a structure with at least:
%              .trackMateSpotsCsv - Full file path to the trackMateSpots CSV file.
%              .posId             - Numeric identifier for the position.
%
%   fi_table - (table) Existing fault indication table with columns "PosID", "DropID", and
%              "Reason" to which new entries will be appended.
%
%   thres    - (scalar double) A threshold fraction for the initial radius value to determine 
%              significant changes.
%
% Output:
%   fi_table - (table) The updated table with flagged tracks appended.
%
% Example:
%   fi_table = table('Size',[0 3],'VariableTypes',{'double','double','cell'},...
%                    'VariableNames',{'PosID','DropID','Reason'});
%   fi_table = radiusDiffDetection(database, fi_table, 0.5);

    arguments
        database {iscell(database)}
        fi_table table
        thres (1,1) double {mustBeGreaterThan(thres, 0)}
    end

    for i = 1:length(database)
        fprintf("Radius Diff Detection: %d / %d\n", i, length(database));

        db = database{i};
        trackMateSpots = readtable(db.trackMateSpotsCsv);
        % Skip the first row if needed, then sort by FRAME and TRACK_ID.
        trackMate = sortrows(sortrows(trackMateSpots(2:end, :), "FRAME"), "TRACK_ID");

        % Get the unique track IDs.
        ids = unique(trackMate.TRACK_ID);

        fi_posid = [];
        fi_dropid = [];
        fi_reason = {};

        for j = 1:length(ids)
            % Extract and sort the data for the current track.
            track = sortrows(trackMate(trackMate.TRACK_ID == ids(j), :), "FRAME");

            radius = track.RADIUS;
            if numel(radius) < 2
                % Skip if insufficient data to compute differences.
                continue;
            end
            
            % Compute differences between consecutive radius values.
            radiusDif = radius(1:end-1) - radius(2:end);

            % Flag the track if the maximum absolute difference exceeds the threshold fraction.
            if max(abs(radiusDif)) > radius(1) * thres
                fi_posid(end+1) = db.posId;
                fi_dropid(end+1) = ids(j);
                fi_reason{end+1} = 'radius change';
            end
        end

        % Create a temporary table for the current position.
        tmp_fi_table = table(fi_posid', fi_dropid', fi_reason', ...
                             'VariableNames', {'PosID', 'DropID', 'Reason'});
        % Append temporary table to the overall fault indication table.
        fi_table = [fi_table; tmp_fi_table];
    end

end
