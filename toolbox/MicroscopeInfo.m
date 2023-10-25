classdef MicroscopeInfo
    % TODO: Add description
    properties (Constant)
        newEpiChannelID = dictionary(["BF", "BFP", "CFP", "FRET", "YFP"], ...
                                     [0, 1, 2, 4, 5]);
        oldEpiChannelID = dictionary(["BF", "CFP", "mCherry", "Custom", "YFP"], ...
                                     [4, 5, 3, 8, 6]);
        micronToPixel = 1645.0 / 2000.0; % 1645 px = 2000 um at resolution 1 for 20x objective
    end
end