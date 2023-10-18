classdef MicroscopeInfo
    % TODO: Add description
    properties (Constant)
        newEpiChannelID = dictionary(["BF", "BFP", "CFP", "FRET", "YFP"], ...
                                     [0, 1, 2, 4, 5]);
        oldEpiChannelID = dictionary(["BF", "CFP", "mCherry", "Custom", "YFP"], ...
                                     [4, 5, 3, 8, 6]);
        % TODO: Add pixel to micron conversion
        % pixelToMicron = 1;
    end
end