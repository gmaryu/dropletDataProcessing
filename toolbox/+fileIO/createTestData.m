function createTestData(srcFolder, destRoot, posVec, frameVec, colorSet)
    % createTestData  Copy test data files from a source folder to a specified destination root folder,
    %                organizing files into subfolders based on position, while keeping all color images together.
    %
    %   createTestData(srcFolder, destRoot, posVec, frameVec, colorSet)
    %
    % Inputs:
    %   srcFolder - (char or string) Path to the source data folder.
    %   destRoot  - (char or string) Path to the destination root folder where test data folders will be created.
    %   posVec    - (numeric vector) Vector of position numbers.
    %   frameVec  - (numeric vector) Vector of frame numbers.
    %   colorSet  - (cell array of strings) Set of color identifiers (e.g., {'4-BF', '5-CFP'}).
    %
    % For each position in posVec, this function:
    %   1. Creates a subfolder within destRoot named "PosX".
    %   2. Searches for files in srcFolder whose names match the position, frame, and color.
    %   3. Copies the matching files into the corresponding position folder, keeping the original file names.
    %
    % Example:
    % imageNameFormats = [
    %     "img_%09d_4-BF_000.tif", ...
    %     "img_%09d_1-DAPI_000.tif", ...
    %     "img_%09d_5-CFP_000.tif", ...
    %     "img_%09d_6-YFP_000.tif", ...
    %     "img_%09d_8-Custom_000.tif", ...
    %     ];
    % fileIO.createTestData("C:\Data\Original", "C:\Data\Test", [3,5], 1:150,imageNameFormats);
    
        srcFolder = string(srcFolder);
        destRoot  = string(destRoot);
        
        for p = posVec
            % Create destination folder for this position
            destFolder = fullfile(destRoot, sprintf('Pos%d', p));
          
            if ~exist(destFolder, 'dir')
                mkdir(destFolder);
                fprintf('Created folder: %s\n', destFolder);
            end
            
            for f = frameVec
                for color = colorSet
                    % Construct the search pattern
                    pattern = sprintf(color{1},f);
                    files = dir(fullfile(srcFolder, sprintf('Pos%d', p), sprintf('*%s*',pattern)));
                    
                    if isempty(files)
                        fprintf('No files found for Pos%d  with pattern: %s\n', ...
                            p, pattern);
                        continue;
                    end
                    srcFile = fullfile(srcFolder, sprintf('Pos%d', p), files(1).name);
                    img = imread(srcFile);
                    img_size = size(img);
   
                    if img_size(1) > 2000 % bin size 1x1
                        %disp(img_size)
                        for k = 1:length(files)
                            srcFile = fullfile(srcFolder, sprintf('Pos%d', p), files(k).name);
                            img = imread(srcFile);
                            scale = 0.5;
                            scaled_img = imresize(img,scale);
                            destFile = fullfile(destFolder, files(k).name);
                            imwrite(scaled_img,destFile)
                            %fprintf('Copied file: %s\n', files(k).name);
                        end
                    else % bin size 2x2
                        % Copy matching files to the destination folder
                        for k = 1:length(files)
                            srcFile = fullfile(srcFolder, sprintf('Pos%d', p), files(k).name);
                            destFile = fullfile(destFolder, files(k).name);
                            copyfile(srcFile, destFile);
                            %fprintf('Copied file: %s\n', files(k).name);
                        end
                    end
                end
            end
            fprintf('Copied files for Position: %d\n', p);
        end
    end