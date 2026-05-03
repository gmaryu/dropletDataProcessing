customFigSpec.LINEWIDTH = 72 / 25.4 * 0.2; % 0.2 mm line weight
customFigSpec.axSpec = {...
    %'LineWidth', customFigSpec.LINEWIDTH, ...
    'FontSize', 6, ... % 6 pt for tick labels
    'FontSizeMode', 'manual', ...
    'LabelFontSizeMultiplier', 7 / 6, ... % 6 * 7 / 6 = 7 pt for axis labels
    'TickLength', [0.02, 0.05], ...
    'XColor', 'k', ...
    'YColor', 'k', ...
    'ZColor', 'k', ...
    'Units', 'centimeters', ...
    };
customFigSpec.colBlue = [50 50 200] / 255;
customFigSpec.colRed = [200 50 50] / 255;
customFigSpec.colGreen = [30 120 30] / 255;
customFigSpec.colPurple = [120 30 200] / 255;

customFigSpec.colRedDim = [200, 120, 120] / 255;
customFigSpec.colBlueDim = [120, 120, 200] / 255;
customFigSpec.colGreenDim = [100 180 100] / 255;
customFigSpec.colPurpleDim = [180 100 200] / 255;

%colcbrewer = [228,26,28; 55,126,184; 77,175,74; 152,78,163; 255,127,0; 255,255,51; 166,86,40] / 255;
%colcbrewer = [27,158,119; 217,95,2; 117,112,179; 231,41,138; 102,166,30; 230,171,2; 102,102,102] / 255;

% 20-color qualitative colormap (0–1 normalized)/High visibility and CVD
% safety considered
customFigSpec.cmap20 = [
    217,  95,   2;   %  1 orange2   - deep warm orange
    230, 159,   0;   %  2 orange1   - lighter orange
    204, 187,  68;   %  3 yellow1   - muted yellow-ochre
    240, 228,  66;   %  4 yellow2   - bright yellow
    153, 153,  51;   %  5 olive1    - olive-toned yellow-green
    102, 166,  30;   %  6 green2    - yellowish green
     34, 136,  51;   %  7 green1    - deep green
      0, 158, 115;   %  8 teal1     - teal (blue-green)
     86, 180, 233;   %  9 blue3     - light sky blue
    102, 194, 255;   % 10 cyan1     - bright cyan-blue
     68, 119, 170;   % 11 blue2     - standard mid-blue
     51,  82, 130;   % 12 blue1     - darker subdued blue
    117, 112, 179;   % 13 purple1   - soft purple
    102,  45, 145;   % 14 purple2   - deep purple
    204, 121, 167;   % 15 red1      - bright pinkish red
    231,  41, 138;   % 16 red2      - magenta-leaning red-purple
    136,  34,  85;   % 17 red3      - dark wine red
    221, 221, 221;   % 18 gray3     - light gray
    153, 153, 153;   % 19 gray1     - medium gray
    102, 102, 102;   % 20 gray2     - dark gray
] ./ 255;


% 30-color qualitative colormap (no gray, 0–1 normalized)
customFigSpec.cmap30 = [
    214,  85,   0;
    230, 120,  20;
    240, 150,  60;
    245, 180,  90;
    220, 190,  40;
    200, 170,  20;
    170, 160,  30;
    140, 150,  40;
    110, 150,  50;
     80, 150,  60;
     50, 140,  70;
     20, 130,  90;
      0, 140, 120;
     40, 160, 160;
     70, 170, 200;
     90, 180, 230;
     70, 150, 220;
     50, 120, 200;
     40,  95, 170;
     60,  80, 150;
     90,  70, 160;
    120,  70, 170;
    150,  70, 170;
    180,  80, 160;
    200,  90, 150;
    220,  90, 130;
    230,  70, 110;
    210,  60,  90;
    180,  50,  70;
    150,  40,  60;
] ./ 255;
