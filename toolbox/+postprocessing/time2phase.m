function cycles = time2phase(cycles, peakMatrix)
%% INPUTS
%{
DESCRIPTION
    Convert time to phase.
INPUTS
    position: scalar
    oscillation_table: 
    sumTrackTable: 
%}

%Frame2min = cycles.START_MINUTE(1)/cycles.START_FRAME(1);
%cycles.END_MINUTE = cycles.END_FRAME .* Frame2min;

%% collect values for phase conversion
% empty vector to put phase data and converted values
phaseConversion.cycle=0;
phaseConversion.phase=0;
phaseConversion.frame=0;

for cycle=1:size(peakMatrix,2)
    target_peak=peakMatrix(:,cycle);

    phaseConversion.cycle(cycle) = cycle;
    phaseConversion.phase(cycle) = cycle*2*pi;
    phaseConversion.frame(cycle) = median(target_peak,'omitnan');
end

tmp_frame=linspace(max(cycles.END_FRAME),1,max(cycles.END_FRAME));
tmp_phase=NaN*ones(length(tmp_frame),1); % vector for phase values (max tracking length)
for tp=1:size(tmp_frame(:),1)
    t=tmp_frame(tp);
    tmp_phase(tp) = interp1(phaseConversion.frame, phaseConversion.phase, t, 'linear','extrap');
end

droplet_POS_ID = unique(cycles(:, {'POS_ID','TRACK_ID'}), 'rows');
for i = 1:size(droplet_POS_ID,1)
    POS = droplet_POS_ID.POS_ID(i);
    ID  = droplet_POS_ID.TRACK_ID(i);
    mask = cycles.POS_ID == POS & cycles.TRACK_ID == ID;
    tmpPeaks = cycles(mask,:);
    if any(~tmpPeaks.IGNORED)
        %disp('do')
        targets = tmpPeaks.START_FRAME;
        [tf, idx] = ismember(targets, tmp_frame);

        out = nan(numel(targets), 1);
        out(tf) = tmp_phase(idx(tf));

        assert(all(tf), 'targets の中に tmp_frame に存在しない値があります。');
        
        cycles.START_PHASE(mask,:) = out./(2*pi);
    end
    
end

%disp('hoge')
end
 
%{
%% extract position information from oscillation table and sumTrackTable
positions_oscillation_table= extractfield(oscillation_table,'position');
positions_sumTrackTable=[];
for i=1:size(sumTrackTable,2)
    tmpTable=sumTrackTable{i};
    tmpPos=extractfield(tmpTable, 'position');
    positions_sumTrackTable=[positions_sumTrackTable,tmpPos(1)];
end

%% check designated position is included in these oscillatino and tracking result
if ismember(position, positions_oscillation_table) && ismember(position, positions_sumTrackTable)
    % extract tables 
    tmp_Otable=oscillation_table(positions_oscillation_table==position).data;
    tmp_Ttable=sumTrackTable{positions_sumTrackTable==position};
    
    %
    tmp_Otable.peak_phase=NaN*zeros(size(tmp_Otable,1),1);

    % get a list of unique dropletID from oscillation table. Sorted
    % oscillation table has new droplet ID, so use prev_IDs in that case
    IDs=tmp_Otable.dropID;
    unique_IDs=unique(IDs);

    %% collect values for phase conversion 
    % empty vector to put phase data and converted values
    phaseConversion.cycle=0;
    phaseConversion.phase=0;
    phaseConversion.frame=0;

    for cycle=1:max(tmp_Otable.cycleID)
        row_index=tmp_Otable.cycleID==cycle;
        target_peak=tmp_Otable(row_index,:);

        phaseConversion.cycle(cycle) = cycle;
        phaseConversion.phase(cycle) = cycle*2*pi;
        phaseConversion.frame(cycle) = median(target_peak.absPeakTime);
    end



    %% conversion
    phaseConvertedTable=cell(1,length(unique_IDs));
    
    tmp_phase=NaN*ones(max(extractfield(tmp_Ttable,'Len')),1); % vector for phase values (max tracking length)
    tmp_frame=linspace(1,max(extractfield(tmp_Ttable,'Len')),max(extractfield(tmp_Ttable,'Len')));
    for tp=1:size(tmp_frame,2)
        t=tmp_frame(tp);
        tmp_phase(tp) = interp1(phaseConversion.frame, phaseConversion.phase, t, 'linear','extrap');

    end


    % loop with droplet IDs from oscillation table
    for uid=1:length(unique_IDs)
        % add tmp_phase to taracking result as a new column
        tmp_id=unique_IDs(uid); %droplet ID


        % signal data extraction
        tracks_IDs=extractfield(tmp_Ttable,'id');
        tmp_track=tmp_Ttable(tracks_IDs==tmp_id).Feat;

        firstTrackFrame=min(tmp_track.t);
        lastTrackFrame=max(tmp_track.t);
        tmp_track.est_phase=tmp_phase(firstTrackFrame:lastTrackFrame);

        phaseconv_1st=find(~isnan(tmp_phase), 1 );
        phaseconv_last=find(~isnan(tmp_phase), 1, 'last' );
        
        raw_corr_signal=tmp_track.FC_ratio_whole(phaseconv_1st:phaseconv_last);
        raw_corr_phase=tmp_track.est_phase(phaseconv_1st:phaseconv_last);
        
        % signal decay correction
        troughsFitResult = linearBaseLineCorrection(tmp_Otable, tmp_id);
        yfit = troughsFitResult(1)*(firstTrackFrame:lastTrackFrame)+troughsFitResult(2);
        signal_corrected = tmp_track.FC_ratio_whole./yfit';



        phaseConvertedTable{tmp_id}.position=position;
        phaseConvertedTable{tmp_id}.dropID=tmp_id;
        phaseConvertedTable{tmp_id}.t=tmp_track.t(find(~isnan(tmp_phase), 1 ):find(~isnan(tmp_phase), 1, 'last' ));
        phaseConvertedTable{tmp_id}.signal=raw_corr_signal;
        phaseConvertedTable{tmp_id}.phase=raw_corr_phase;
        phaseConvertedTable{tmp_id}.crrectedSignal=signal_corrected;
        
        % add phase of each peak to oscillation table
        row_index=IDs==tmp_id;
        target_peak=tmp_Otable(row_index,:);
        target_peak.peak_phase=tmp_phase(target_peak.absPeakTime);
        tmp_Otable(row_index,:)=target_peak;
        phaseConvertedTable{tmp_id}.peaks=target_peak.absPeakTime;
        
    end
    
    sumTrackTable{positions_sumTrackTable==position} = tmp_Ttable;
    oscillation_table(positions_oscillation_table==position).data=tmp_Otable;
   
else
    disp('Error')
end
%}