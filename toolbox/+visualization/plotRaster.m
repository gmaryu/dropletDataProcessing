function plotRaster(data, varargin)
% plotRaster  簡易ラスター図を描く
%
% Usage
%   plotRaster(data)
%   plotRaster(data, Name,Value, ...)
%
% Optional Name–Value pairs
%   visibleLine    "on" | "off"   既定:"off"
%   Color          1×3 RGB        既定:[0 0 0]
%   LineWidth      正の数          既定:0.5

%% ---- 引数定義とデフォルト ----
arguments
    data struct                    % 必須
    visibleLine (1,1) string   {mustBeMember(visibleLine,["on","off"])} = "off"
    Color (1,3) double       {mustBeGreaterThanOrEqual(Color,0),mustBeLessThanOrEqual(Color,1)} = [0 0 0]
    LineWidth (1,1) double   {mustBePositive} = 0.5
end

%% ---- プロット本体 ----

t = data.info.POS_ID;                      % 例: 1列目が時刻
y = data.info.TRACK_ID;                  % 2列目以降がスパイク(0/1)
hold on
for k = 1:size(y,2)
    tk = t(logical(y(:,k)));
    plot(tk, k*ones(size(tk)), '.', ...
         'Color',Color, 'MarkerSize',4);
end
if visibleLine == "on"
    yline(0,'-', 'Color',Color, 'LineWidth',LineWidth);
end
hold off
end
