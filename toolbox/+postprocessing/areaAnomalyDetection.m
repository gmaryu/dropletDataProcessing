function fi_table = areaAnomalyDetection(database, fi_table, thres)

    arguments
        database {iscell(database)}
        fi_table table
        thres (1,1) double {mustBeGreaterThan(thres, 0)}
    end

    for i = 1:length(database)
        fprintf("Area Anomaly Detection: %d / %d\n", i, length(database));

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

            areaAnomaly = abs((track.AREA - median(track.AREA)) / median(track.AREA));
            if sum(areaAnomaly > 0.05) > 0.01 * numel(areaAnomaly)
                validIdx = areaAnomaly < 0.05;
                origLen = height(track);
                track = track(validIdx, :);
                if height(track) < thres * origLen
                    %fprintf("ID:%d - Area anomaly (%.3f, %.3f, %d)\n", ids(i), median(areaAnomaly), std(areaAnomaly), sum(areaAnomaly > 0.1));
                    fi_posid(end+1) = db.posId;
                    fi_dropid(end+1) = ids(j);
                    fi_reason{end+1} = 'Area Anomaly';
                    continue;
                end
            end

        end

        % Create a temporary table for the current position.
        tmp_fi_table = table(fi_posid', fi_dropid', fi_reason', ...
                             'VariableNames', {'PosID', 'DropID', 'Reason'});
        % Append temporary table to the overall fault indication table.
        fi_table = [fi_table; tmp_fi_table];
    end
end