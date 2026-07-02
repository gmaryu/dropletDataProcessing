% Private helper shared within +postprocessing (see extractShapeFeatures.m).
function dropletID = makeDropletID(positionID, trackID)
    dropletID = strcat("P", string(positionID), "_D", string(trackID));
end

