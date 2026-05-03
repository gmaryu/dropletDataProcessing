function plotLineScreening(dataSet, database, varargin)

% ----- parse options -----
p = inputParser;
addParameter(p,'LineVar1',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'LineVar2',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'LineVar3',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'ImgCh1',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'ImgCh2',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'AdditionalExportFolderName',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'IgnoreVar',"IGNORED",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'Colors',[],@(c)isnumeric(c)&&size(c,2)==3);
addParameter(p,'MarkerSize',8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TileLayout','auto',@(v)(ischar(v)||isstring(v)) || (isnumeric(v)&&numel(v)==2));
parse(p,varargin{:});
opt = p.Results;

LineVar1   = string(opt.LineVar1);
LineVar2   = string(opt.LineVar2);
LineVar3   = string(opt.LineVar3); 
ImgCh1   = string(opt.ImgCh1);
ImgCh2 = string(opt.ImgCh2);
fldName = string(opt.AdditionalExportFolderName);
ignName = string(opt.IgnoreVar);


if strlength(LineVar1)==0
    error('LineVar is required (char or string scalar).');
end

%% ----- data extraction ------
info = dataSet.info;
cycle = dataSet.cycle;
ts = dataSet.timeSeries;
posIDs = cellfun(@(s) s.posId, database);


%% ----- loop by individual droplet ------
for i = 1:size(info,1)
    if ~info.IGNORED(i)
        tmpPos = info.POS_ID(i);
        tmpTrackID = info.TRACK_ID(i);
        tmpMedDia = info.ORIGINAL_MED_DIAMETER(i);
        fprintf("Pos%d TrackID%d", tmpPos, tmpTrackID);

        % index in database
        idx = find(posIDs == tmpPos);
        if ~isempty(idx)
            tmpStruct = database{idx};
            croppedImages = tmpStruct.croppedImages;
            
            tmpCyc = cycle(cycle.POS_ID == tmpPos & cycle.TRACK_ID== tmpTrackID, :);
            tmpTS = ts(ts.POS_ID == tmpPos & ts.TRACK_ID== tmpTrackID, :);
            

            % line plot data
            tmpSignal = tmpTS.(LineVar1);
            tmpSignal2 = tmpTS.(LineVar2);
            tmpSignal3 = tmpTS.(LineVar3);
            tmpTime = tmpTS.('FRAME')*dataSet.FrameToMin;
           
            % images
            numCyc = size(tmpCyc,1);

            % replication limit
            rlimit = info.rmax(i);
            rlimit_DNAVR = info.rmax_DNACR(i);
            
            % firstFrame = min(tmpTS.('FRAME'));
            % tmpImg1 = imread(fullfile(croppedImages, sprintf('droplet_%03d',tmpTrackID),sprintf('Pos%d_%s_%03d.tif',tmpPos,ImgCh1,firstFrame)));
            % tmpImg2 = imread(fullfile(croppedImages, sprintf('droplet_%03d',tmpTrackID),sprintf('Pos%d_%s_%03d.tif',tmpPos,ImgCh2,firstFrame)));


            % make figure panel
            f = figure('Visible','off');
            f.Units = 'centimeters';
            f.Position = [0,0,10,20];

            % 全体を 1行2列 に分ける
            t = tiledlayout(4,1,"TileSpacing","compact","Padding","compact");

            % 左側にラインプロット
            nexttile(1);
            plot(tmpTime,tmpSignal);
            title(sprintf('Diameter: %d', ceil(tmpMedDia)));

            nexttile(2);
            plot(tmpTime,tmpSignal2);
            xline(unique([tmpCyc.START_INDEX*dataSet.FrameToMin;tmpCyc.END_INDEX*dataSet.FrameToMin]),':');
            ylabel('NCVR')

            nexttile(3);
            plot(tmpTime,tmpSignal3);
            yline(rlimit_DNAVR,'r');
            yline(2.^(0:fix(rlimit))*tmpCyc.DNACR_NUC(1),'--');
            xline(unique([tmpCyc.START_INDEX*dataSet.FrameToMin;tmpCyc.END_INDEX*dataSet.FrameToMin]),':');
            %ylim([0,rlimit_DNAVR+0.1]);
            title(sprintf('ReplicationLimit: %d',rlimit));
            ylabel('DNACR')


            % image panel
            imPanel1 = tiledlayout(t,2,numCyc+1);
            imPanel1.Layout.Tile = 4; % 
            imPanel1.TileSpacing = "tight";
            imPanel1.Padding = "tight";

            for j = 1:size(tmpCyc,1)+1
                if j == 1
                   Frame = min(tmpTS.('FRAME'));
                else
                   Frame = tmpCyc.END_FRAME(j-1) - 3;
                end
    
                img1Path = fullfile(croppedImages, sprintf('droplet_%03d',tmpTrackID),sprintf('Pos%d_%s_%03d.tif',tmpPos,ImgCh1,Frame));
                img2Path = fullfile(croppedImages, sprintf('droplet_%03d',tmpTrackID),sprintf('Pos%d_%s_%03d.tif',tmpPos,ImgCh2,Frame));
                
                if isfile(img1Path)
                    tmpImg1 = imread(img1Path);
                end
                if isfile(img2Path)
                    tmpImg2 = imread(img2Path);
                end

                % % 右上（実際には左側のタイル）
                nexttile(imPanel1,j);
                imshow(tmpImg1,[]);
                title(sprintf('%d min',Frame*dataSet.FrameToMin));

                % 右下（実際には右側のタイル）
                nexttile(imPanel1,numCyc+1+j);
                imshow(tmpImg2,[]);
                %title(ImgCh2);

            end

            [exportfolder, name, ext] = fileparts(tmpStruct.forceIgnoreCsv);
            outfolder = fullfile(exportfolder,'LineScreening',fldName);
            
            if ~exist(outfolder, 'dir')
                mkdir(outfolder);
            end

            outfn = fullfile(outfolder,sprintf('Pos%d_Track%03d.png',tmpPos,tmpTrackID));            
            exportgraphics(f, outfn);
            close gcf

        else
            warning('PosID=%d was not found', tmpPos);
            break
        end
    end

end


end



%%
