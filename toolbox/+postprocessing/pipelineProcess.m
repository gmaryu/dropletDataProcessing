function data = pipelineProcess(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
                                spermCondition, hoechstCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
                                automaticNucleiCount, hoechstoffset, FRETNumerator, FRETDenominator)
% pipelineProcess  Run the full analysis pipeline across multiple positions.
%
%   data = pipelineProcess(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
%                          spermCondition, hoechstCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
%                          automaticNucleiCount, hoechstoffset, FRETNumerator, FRETDenominator)
%
% Description:
%   Executes the end‑to‑end image analysis workflow on a collection of positions.  
%   For each entry in the input `database` cell array whose `posId` is listed in
%   `totalPositions`, it calls `processPosition` to perform:
%     • TrackMate data parsing and oscillation peak detection  
%     • Optional nuclear segmentation and DNA quantification  
%     • Per-droplet oscillation and morphology analysis  
%   After processing all valid positions, it merges the individual results into a
%   single output structure via `mergeDatabase`.
%
% Inputs:
%   database            – Cell array of structs; each struct must contain at least:
%                           .posId           (numeric) position identifier  
%                           .trackMateSpotsCsv, 
%                           .segmented_tracks path, etc.  
%   totalPositions      – Vector of numeric position IDs to process  
%   frameToMin          – Scalar double; conversion factor from frame index to minutes  
%   pixelToUm           – Scalar double; conversion factor from pixels to micrometers  
%   initialPeakTimeBound– Scalar double; time (min) beyond which oscillations are ignored  
%   forceIgnore         – Table with columns PosID and DropID listing droplets to skip  
%   spermCondition      – Logical; if true, perform nuclear segmentation & DNA mask steps  
%   hoechstCondition    – Logical; if true, include Hoechst mask in downstream analysis  
%   nucChannel          – String; name of the nuclear channel (e.g. "CFP")  
%   dnaChannel          – String; name of the DNA channel (e.g. "DAPI")  
%   overwriteNucMask    – Logical; if true, regenerate nuclear mask .mat files  
%   overwriteDNAInfo    – Logical; if true, regenerate DNA mask .mat files  
%   automaticNucleiCount– Logical; if true, perform automatic nuclei counting  
%   hoechstoffset       – Logical; if false, apply background offset correction to Hoechst signal  
%   FRETNumerator       – String; name of field for FRET numerator  ("FRET" channel)
%   FRETDenominator     – String; name of field for FRET denominator  ("CFP" channel)
%
% Output:
%   data – Struct with fields:
%       • data.info        : Merged droplet‐level summary table  
%       • data.timeSeries  : Merged tracking time series table  
%       • data.cycle       : Merged oscillation cycle table  
%       • data.FrameToMin  : frameToMin (echoed)  
%       • data.PixelToUm   : pixelToUm (echoed)  
%       • data.InitialPeakTimeBound : initialPeakTimeBound (echoed)  
%
% Workflow:
%   1. Iterate through each element of `database`.  
%   2. If `db.posId` ∈ `totalPositions`, call:
%        db = postprocessing.processPosition(db, ...parameters...);  
%      to process that position.  
%   3. Store the updated `db` back into `database{i}`.  
%   4. After the loop, merge all positions via:
%        data = postprocessing.mergeDatabase(database, ...parameters...);  
%
% Example:
%   dbAll = loadDatabase(...);  
%   forceIgnore = readtable('force_ignore.csv');  
%   data = pipelineProcess(dbAll, 0:10, 6, 2.649, 100, forceIgnore, ...
%                          true, true, "CFP", "DAPI", true, true, true, true, ...
%                          "MEAN_INTENSITY_CH5", "MEAN_INTENSITY_CH3");
%
% See also:
%   postprocessing.processPosition, postprocessing.mergeDatabase


arguments
    database            cell
    totalPositions      {mustBeNumericOrLogical}
    frameToMin          double
    pixelToUm           double
    initialPeakTimeBound double
    forceIgnore         table
    spermCondition logical
    hoechstCondition logical
    nucChannel string
    dnaChannel string
    overwriteNucMask logical
    overwriteDNAInfo logical
    automaticNucleiCount logical
    hoechstoffset logical
    FRETNumerator string
    FRETDenominator string
end

% Loop through each database entry and process only the selected positions.
%%parfor i = 1:length(database)
for i = 1:length(database)

    db = database{i};
    if ismember(db.posId, totalPositions)
        db = postprocessing.processPosition(db, frameToMin, pixelToUm, initialPeakTimeBound, forceIgnore, ...
            spermCondition, hoechstCondition, nucChannel, dnaChannel, overwriteNucMask, overwriteDNAInfo, ...
            automaticNucleiCount, hoechstoffset, FRETNumerator, FRETDenominator);
        database{i} = db; % Update the database entry
    end
end

% Merge results across positions into final data structure.
data = postprocessing.mergeDatabase(database, totalPositions, frameToMin, pixelToUm, initialPeakTimeBound);
end