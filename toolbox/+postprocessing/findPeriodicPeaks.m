function retv = findPeriodicPeaks(signal, frameToMin)
% findPeriodicPeaks Detects periodic peaks in a signal.
%
%   retv = findPeriodicPeaks(signal, frameToMin)
%
% The function uses findpeaks on both the signal and its negative to determine
% peak and trough positions. It returns a matrix with each row as
% [start_index, end_index, trough_index] if the peaks and troughs match the
% expected pattern; otherwise, it returns NaN.

    signal = signal(:);

    p = 0.08;
    maxw = 180 / frameToMin;  % Maximum expected peak width in frames

    % Main peak/trough detection
    [~, ip] = findpeaks(signal, "MinPeakProminence", p, "MaxPeakWidth", maxw);
    [~, it] = findpeaks(-signal, "MinPeakProminence", p);

    % Rescue a missed first peak WITHOUT mirroring/reflection.
    % Only do this if:
    %   - there is at least one trough,
    %   - no detected peak exists before the first trough,
    %   - and the original signal itself supports a rise-then-fall pattern.
    if ~isempty(it)
        firstTrough = it(1);

        hasPeakBeforeFirstTrough = any(ip < firstTrough);

        if ~hasPeakBeforeFirstTrough && firstTrough >= 3
            searchIdx = 1:firstTrough;
            [candVal, candRelIdx] = max(signal(searchIdx));
            candIdx = searchIdx(candRelIdx);

            % Never rescue frame 1 as a peak.
            if candIdx > 1 && candIdx < firstTrough
                % Require local maximum shape in the original signal
                isLocalMax = signal(candIdx) >= signal(candIdx-1) && ...
                             signal(candIdx) >  signal(candIdx+1);

                % Require evidence that the signal rose before the candidate
                leftRise = candVal - signal(1);

                % Require evidence that it fell afterward into the first trough
                rightDrop = candVal - signal(firstTrough);

                % Optional: require candidate width not to be absurdly broad/narrow
                spanToTrough = firstTrough - candIdx;

                if isLocalMax && ...
                   leftRise >= p/2 && ...
                   rightDrop >= p && ...
                   spanToTrough <= maxw

                    ip = sort([ip(:); candIdx]);
                end
            end
        end
    end

    % Match peaks and troughs into cycles
    if numel(ip) == numel(it)
        if all(it - ip > 0)
            retv = [ip(1:end-1), ip(2:end), it(1:end-1)];
        elseif all(it - ip < 0)
            retv = [ip(1:end-1), ip(2:end), it(2:end)];
        else
            retv = nan;
        end
    elseif numel(ip) == numel(it) + 1
        if all(it - ip(1:end-1) > 0)
            retv = [ip(1:end-1), ip(2:end), it(1:end)];
        else
            retv = nan;
        end
    elseif numel(ip) + 1 == numel(it)
        % discard first trough
        if all(it(2:end) - ip(1:end) > 0)
            retv = [ip(1:end-1), ip(2:end), it(2:end-1)];
        else
            retv = nan;
        end
    else
        retv = nan;
    end
end