function plotLinePlusImages(dataSet, database, varargin)

% ----- parse options -----
p = inputParser;
addParameter(p,'LineVar1',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(p,'LineVar2',"",@(x)ischar(x)||(isstring(x)&&isscalar(x)));
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
ImgCh1   = string(opt.ImgCh1);
ImgCh2 = string(opt.ImgCh2);
fldName = string(opt.AdditionalExportFolderName);
ignName = string(opt.IgnoreVar);


if strlength(LineVar1)==0
    error('LineVar is required (char or string scalar).');
end

% ----- data extraction ------
info = dataSet.info;
ts = dataSet.timeSeries;
posIDs = cellfun(@(s) s.posId, database);


%% ----- loop by individual droplet ------
for i = 1:size(info,1)
    if ~info.IGNORED(i)
        tmpPos = info.POS_ID(i);
        tmpTrackID = info.TRACK_ID(i);
        tmpMedDia = info.ORIGINAL_MED_DIAMETER(i);

        % 一致するセルのインデックスを探す
        idx = find(posIDs == tmpPos);
        if ~isempty(idx)
            tmpStruct = database{idx};
            croppedImages = tmpStruct.croppedImages;
            
            tmpTS = ts(ts.POS_ID == tmpPos & ts.TRACK_ID== tmpTrackID, :);
            
            % line plot data
            tmpSignal = tmpTS.(LineVar1);
            tmpTime = tmpTS.('FRAME')*dataSet.FrameToMin;
           
            % images
            firstFrame = min(tmpTS.('FRAME'));
            tmpImg1 = imread(fullfile(croppedImages, sprintf('droplet_%03d',tmpTrackID),sprintf('Pos%d_%s_%03d.tif',tmpPos,ImgCh1,firstFrame)));
            tmpImg2 = imread(fullfile(croppedImages, sprintf('droplet_%03d',tmpTrackID),sprintf('Pos%d_%s_%03d.tif',tmpPos,ImgCh2,firstFrame)));


            % make figure panel
            f = figure('Visible','off');
            f.Units = 'centimeters';
            f.Position = [0,0,10,3];

            % 全体を 1行2列 に分ける
            t = tiledlayout(1,2,"TileSpacing","compact","Padding","compact");

            % 左側にラインプロット
            nexttile(1);
            plot(tmpTime,tmpSignal);
            title(sprintf('Diameter: %d', tmpMedDia));
            if ~isempty(LineVar2)
                hold on
                yyaxis right
                tmpSignal2 = tmpTS.(LineVar2);
                plot(tmpTime,tmpSignal2);
                hold off
            end

            % 右半分をさらに 1行2列 に分割
            tRight = tiledlayout(t,1,2);
            tRight.Layout.Tile = 2; % 右側に配置
            tRight.TileSpacing = "compact";
            tRight.Padding = "compact";

            % 右上（実際には左側のタイル）
            nexttile(tRight,1);
            imshow(tmpImg1,[]);
            title(ImgCh1);

            % 右下（実際には右側のタイル）
            nexttile(tRight,2);
            imshow(tmpImg2,[]);
            title(ImgCh2);

            [exportfolder, name, ext] = fileparts(tmpStruct.forceIgnoreCsv);
            outfolder = fullfile(exportfolder,'LinePlusImages',fldName);
            
            if ~exist(outfolder, 'dir')
                mkdir(outfolder);
            end

            outfn = fullfile(outfolder,sprintf('Pos%d_Track%03d.png',tmpPos,tmpTrackID));            
            exportgraphics(f, outfn);
            close gcf

        else
            warning('PosID=%d was not found', targetID);
            break
        end
    end

end


end



%%
