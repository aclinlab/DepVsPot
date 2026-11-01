%% --- Fly-level batch extraction of response traces from 5D ScanImage data ---
%
% NB: Currently only extracts data from 2 z-slices (4 ROIs) but can be adapted to be more
% flexible.
%
% Dependencies:
% - requires Matlab's Image Processing Toolbox
% - needs 'readScanImageTiffLinLab' for data loading (add to path)
%   - needs ScanImageTiffReader.m (download from https://vidriotech.gitlab.io/scanimage-tiff-reader/)
% - needs stimTimes.mat file, with stimulus start times for experiment (in
% seconds)
%
% Output:
% DetailedData.mat file for each fly, which contains all data needed for subsequent
% analysis and plotting (MultNoise_StatsAndPlotting+Deming_ErrorVariance.m)
%
% This script:
% - Loops through the selected folder and loads all trials per fly (in ascending order)
% - Uses ScanImageTiffReader for fast loading of .tif files
% - Saves all SI metadata in a cell array
% - Opens GUI for user to specify:
%       - which channel to analyse
%       - smoothing window
%       - response window 
%       - fixed and rolling F0 baseline window
% - Creates two non-modal 'implay' windows (original and smoothed)
% - Creates window with smoothed average projection of the 2 z-slices
%       - Opens two separate contrast/brightness adjustment windows for each
%         projection
%       - Click once to confirm contrast, zoom etc, 
%       - Then manually draw the ROIs
%           - 2 per window (1  = background, 2 = signal)
% - Saves all ROI masks and corresponding z-slice indices in cell array, to
% correctly map each ROI back to z-plane during trace extraction
%       - allROIs and all_roiSlices 
% - Opens 'live' dF/F feedback window 
%   (using fixed F0, to check if ROIs and esp. background, are selected correctly)
% - Runs a Monte Carlo Split-Half Resampling (Repeated Random Sampling) 
%   to detect photon shot noise (not physical poisson noise - because in DF/F space)
%       - Shot noise resampling estimates photon shot noise by randomly splitting each
%         ROI into two halves and computing peak responses independently. 
%         The resulting distributions are used in subsequent statistical analysis
%         (MultNoise_StatsAndPlotting_DEming_ErrorVariance.m) to build a 
%         sample-size matched shot noise null for the TLS-based multiplicative
%         noise test.
% - Extracts peak stimulus responses for each ROI (all_amplitudes_ROI1/2)
%       - subtracts background from raw traces
%       - smoothes and calculates DF/F with rolling baseline
%       - extracts peak amplitude during specified response window
%       - saves results in 5x10 matrix (adds for every trial) for each ROI
%         
% - Saves progress in original fly folder to allow resuming analysis (if
%   exists, prompts user to 'start fresh' or 'resume')
% - Saves key output variables in DetailedData.mat for each fly 
%   (in master output folder - hardcoded), to be used for subsequent statistical 
%   analysis and plotting
%
% NB: MATLAB's image display functions (like implay and imshow)
% expect matrices to be in (row, column) or (Y, X) format, which corresponds 
% to [Height, Width]. ScanImage TIF files save data in the order it was acquired: 
% (X, Y), or [Width, Height], because X axis (pixels per line) is the "fast" scanning axis.
%
% Avoid using 'squeeze' to remove singelton dimension. Keep dimension order rigid, otherwise
% 1 z-slice and 2 channels, vs 2 z-slices and 1 channel would be the same.

% Written by Katharina Bergmann

clear; clc; close all hidden;

%% --------------------------------------------------------------------------------- 
%                           1. SETUP & PARAMETERS 
% ----------------------------------------------------------------------------------

% --- 1) FOLDER SELECTION --- 

if ~exist('readScanImageTiffLinLab', 'file') % needed to load .tif file efficiently
    error('"readScanImageTiffLinLab.m" not found. Please add it to your MATLAB path.');
end

fprintf('\nSelect the folder containing the TIF files for this fly.\n');
folderPath = uigetdir('', 'Select the Folder Containing Your TIF Trials');
if isequal(folderPath, 0), disp('User selected Cancel');
    return;
end
% Get the files
fileList = dir(fullfile(folderPath, '*.tif'));
if isempty(fileList), error('No .tif files were found in the selected folder.');
end


% --- 2) SORTING .TIF FILES FOR: "YYYYMMDD_FlyX_TrialY_..." ---
% - Iterates directly over fileList 
% - Extracts digits specifically after "Trial"

numTrials = numel(fileList);
trialNumbers = zeros(1, numTrials);

for i = 1:numTrials
    currentName = fileList(i).name;
    
    % Regex: Find "Trial" followed immediately by digits (\d+)
    tokens = regexp(currentName, 'Trial(\d+)', 'tokens');
    
    if ~isempty(tokens)
        % Unpack the nested cell: {{'10'}} -> '10' -> 10
        trialNumbers(i) = str2double(tokens{1}{1});
    else
        % Fallback: keep original index if pattern not found
        trialNumbers(i) = i; 
        fprintf('Warning: Could not extract Trial number from: %s\n', currentName);
    end
end

% Sort fileList based on the extracted numbers
[~, sortedIndices] = sort(trialNumbers);
fileList = fileList(sortedIndices);
actualTrialNumbers = trialNumbers(sortedIndices);

% CHECK ORDER: Print the sorted order
fprintf('\n--- File Sorting Check ---\n');
if numel(fileList) >= 1, fprintf(' First: %s (Trial %d)\n', fileList(1).name, trialNumbers(sortedIndices(1))); end
if numel(fileList) > 1,  fprintf(' Last:  %s (Trial %d)\n', fileList(end).name, trialNumbers(sortedIndices(end))); end
fprintf('--------------------------\n');

% --- 3) SETUP OUTPUT STORAGE ---
% Defines where you want to keep your analyzed data
% NB: This path is hard-coded as it refers to one master folder. 
% Progress file saved in Fly Raw Imaging Data folder (allows trial-level
% resume per fly)
% DetailedData saved in master folder (accumulates across flies for
% population analysis)

choice = questdlg('Save output to existing master folder or create new?', ...
                  'Output Folder', ...
                  'Existing','New','Existing');

switch choice
    case 'Existing'
        masterOutputFolder = uigetdir(pwd, 'Select output folder');
        
    case 'New'
        baseFolder = uigetdir(pwd, 'Select parent folder');
        folderName = inputdlg('Enter new folder name:');
        
        if isempty(baseFolder) || isempty(folderName)
            error('Folder selection cancelled.');
        end
        
        masterOutputFolder = fullfile(baseFolder, folderName{1});
        mkdir(masterOutputFolder);
end

if masterOutputFolder == 0
    error('No folder selected.');
end

% Create unique identifier (flyID) from the folder name
[~, flyID, ~] = fileparts(folderPath);
fprintf('Unique Fly Identifier: %s\n', flyID);

% Define Filenames for this specific fly
detailedDataFile = fullfile(masterOutputFolder, [flyID '_DetailedData.mat']); % Contains all relevant information for that fly

% --- 4) PARAMETER INPUT ---
% Create dialogue box to set analysis parameters (answer = inputdlg(prompt,dlgtitle,fieldsize,definput,opts))
prompt = {
    'Which channel to analyze? (e.g., 1 or 2)', ...
    'Enter smoothing window (seconds):', ...
    'Response Window ROI 1 (seconds):', ...
    'Response Window ROI 2 (seconds):', ...
    'Fixed F0 period (seconds, e.g., "0 5"):',...
    'Rolling F0 period(seconds)'...
    };
definput = {'1', '0.5', '2', '2', '0 5', '5'}; % changed to 5s response window for ROI1
answer = inputdlg(prompt, 'Analysis Parameters for All Trials', [1 60], definput);
channelToAnalyze = str2double(answer{1});
smoothing_window_sec = str2double(answer{2});
response_window_sec_roi1 = str2double(answer{3});
response_window_sec_roi2 = str2double(answer{4});
fixed_F0_sec = str2num(answer{5}); % Use str2num to get an array
if isempty(fixed_F0_sec) || numel(fixed_F0_sec) ~= 2
    error('Fixed F0 period must be two numbers (e.g., "0 5").');
end

rolling_F0_sec = str2double(answer{6});

% Select and load stimTimes file
fprintf('\nSelect your stimTimes.mat file\n');
[stimFile, stimPath] = uigetfile('*.mat', 'Select your stimTimes.mat file'); % You need to have a .mat file with your stimulus onset times
if isequal(stimFile, 0)
    disp('User selected Cancel');
    return;
end

stimData = load(fullfile(stimPath, stimFile));

if ~isfield(stimData, 'stimTimes')
    error('The .mat file must contain "stimTimes".');
end

stimTimes_sec = stimData.stimTimes;
numStim = length(stimData.stimTimes);

% Initialize accumulators
all_amplitudes_ROI1 = []; 
all_amplitudes_ROI2 = [];
all_traces_ROI1 = []; 
all_traces_ROI2 = [];
all_SI = cell(numTrials, 1); % To store all SI metadata
all_ROIs = cell(numTrials, 1); % To store all ROI masks
all_roiSlices = cell(numTrials,1);
all_shotnoise_A1 = [];
all_shotnoise_B1 = [];
all_shotnoise_A2 = [];
all_shotnoise_B2 = [];

startTrial = 1;
saveFilePath = fullfile(folderPath, 'analysis_progress.mat'); % in raw data folder containing .tifs
auto_accept_saved = false; % Reset flag at start

% --- 5) PRE-LOAD ROIS IF AVAILABLE ---
% Avoid having to re-draw ROIs while developing the analysis script 
% Loads ROIs for all trials

if exist(detailedDataFile, 'file')
    % Check if the file actually has ROIs
    tempInfo = whos('-file', detailedDataFile);
    if ismember('all_ROIs', {tempInfo.name})
        
        answer = questdlg('Found previous analysis file. Load saved ROIs to skip drawing?', ...
            'Load ROIs', 'Yes, Load ROIs', 'No, Keep Empty', 'Yes, Load ROIs');
        
        if strcmp(answer, 'Yes, Load ROIs')
            fprintf('Loading saved ROIs from: %s\n', detailedDataFile);
            tempData = load(detailedDataFile, 'all_ROIs', 'all_roiSlices');
            
            % Fill the ROI array 
            n_load = min(length(tempData.all_ROIs), numTrials);
            all_ROIs(1:n_load) = tempData.all_ROIs(1:n_load);
            if isfield(tempData, 'all_roiSlices')
                all_roiSlices(1:n_load) = tempData.all_roiSlices(1:n_load);
            end
            fprintf(' -> ROIs loaded for %d trials.\n', n_load);
        end
    end
end

% --- 6) RESUME PROGRESS
% Allows resuming analysis after interruption, e.g., half way through
% trials for one fly
%
% NB: On resume, finished trials keep the ROIs they were analysed with
if exist(saveFilePath, 'file')
    answer = questdlg('Found existing progress file. Resume analysis?', ...
        'Resume Analysis', 'Yes, Resume', 'No, Start Over', 'Yes, Resume'); % last inout to questdlg sets default button

    if strcmp(answer, 'Yes, Resume')
        fprintf('Loading progress from: %s\n', saveFilePath);
        progressData = load(saveFilePath);
        
        % Load all data variables to restore 
        all_amplitudes_ROI1 = progressData.all_amplitudes_ROI1;
        all_amplitudes_ROI2 = progressData.all_amplitudes_ROI2;
   
        all_traces_ROI1     = progressData.all_traces_ROI1;
        all_traces_ROI2     = progressData.all_traces_ROI2;
        all_SI              = progressData.all_SI;
        all_shotnoise_A1    = progressData.all_shotnoise_A1;
        all_shotnoise_B1    = progressData.all_shotnoise_B1;
        all_shotnoise_A2    = progressData.all_shotnoise_A2;
        all_shotnoise_B2    = progressData.all_shotnoise_B2;
        stimFrames          = progressData.stimFrames;
        
        % ROI merge on resume:
        % Finished trials (1:lastTrialDone) -> ROIs from the progress file (the ROIs these results were computed with)
        % Remaining trials                  -> keep any ROIs pre-loaded in step 5 (still confirmed or redrawn per trial)
        
        % Find where you left off: each completed trial adds one column to all_amplitudes_ROI1
        lastTrialDone = size(progressData.all_amplitudes_ROI1, 2);
        
        if lastTrialDone == 0
            startTrial = 1;
        else
            startTrial = lastTrialDone + 1;
            % Overwrite the "done" trials with the specific ROIs used during that run
            all_ROIs(1:lastTrialDone) = progressData.all_ROIs(1:lastTrialDone);
            if isfield(progressData, 'all_roiSlices')
                all_roiSlices(1:lastTrialDone) = progressData.all_roiSlices(1:lastTrialDone);
            end
        end
        
        fprintf(' -> Resuming at Trial %d.\n', startTrial);
    else
        fprintf(' -> Starting fresh (Existing progress ignored).\n');
        % Note: If you loaded ROIs in Step 5, they are still preserved here

    end
end




%% -------------------------------------------------------------------------------
%               2. MAIN LOOP: DATA EXTRACTION & SHOT NOISE RESAMPLING 

% ---------------------------------------------------------------------------------

% --- LOOP THROUGH ALL AVAILABLE TRIALS PER FLY ---
for trial_idx = startTrial:numTrials

    % --- 7) LOAD & DISPLAY DATA ---
    currentFileName = fileList(trial_idx).name;
    fullFilePath = fullfile(folderPath, currentFileName);
    fprintf('\n--- Processing Trial %d of %d: %s ---\n', trial_idx, numTrials, currentFileName);

    fprintf('Loading and reshaping ScanImage data with your function...');
    [result, SI] = readScanImageTiffLinLab(fullFilePath);

    all_SI{trial_idx} = SI;
    fprintf(' Done.\n');

    height = SI.hRoiManager.linesPerFrame;
    width = SI.hRoiManager.pixelsPerLine;
    numZ = SI.hStackManager.numSlices;

    % --- HANDLING 3-SLICE FILES ---
    % Script handles exactly 2 z-slices. If 3 slices present (due to
    % acquisition settings, user selects which 2 to keep)

    if numZ == 3
        fprintf('  *** Alert: Trial has 3 Z-slices. ***\n');
        
        % Ask user which slices to keep 
        quest = sprintf('Trial "%s" has 3 slices. Pick 2 to analyze:', currentFileName);
        choice = questdlg(quest, '3-Slice Detected', ...
            'Slices 1 & 2', 'Slices 1 & 3', 'Slices 2 & 3', 'Slices 1 & 2');
        
        switch choice
            case 'Slices 1 & 2', keepIdx = [1 2];
            case 'Slices 1 & 3', keepIdx = [1 3];
            case 'Slices 2 & 3', keepIdx = [2 3];
            otherwise, error('Operation cancelled by user.');
        end
        
        % Remove the unwanted slice from the 5D matrix
        % Dimensions: [Rows, Cols, Z, Channels, Time]
        result = result(:, :, keepIdx, :, :);
        
        % Reset numZ to 2 so the rest of the script runs normally
        numZ = 2; 
        fprintf('  -> Retained slices %d & %d. Data reshaped to 2 slices.\n', keepIdx(1), keepIdx(2));
        
    elseif numZ ~= 2
        error('Script is designed for 2 z-slices. Found %d.', numZ);
    end
    numChannels = numel(SI.hChannels.channelSave);

    scanFrameRate = SI.hRoiManager.scanFrameRate;
    volumeRate = SI.hRoiManager.scanVolumeRate;
    time_per_volume_sec = 1 / volumeRate;
    fprintf('ScanImage Frame Rate: %.2f Hz\n', scanFrameRate);
    fprintf('Calculated Volume Rate (%d slices): %.2f Hz (%.4f s/vol)\n', numZ, volumeRate, time_per_volume_sec);

    numTime = size(result, 5); % number of frames in time, 5th dimension = time

    % --- MOVIE PLAYER ---
    % Opens two movie player windows for manual inspection of time series, 
    % one with raw data movie, the other smoothed
    disp('Preparing movies... this may take a moment.');
    smoothing_window_frames = round(smoothing_window_sec / time_per_volume_sec);
    if mod(smoothing_window_frames, 2) == 0, smoothing_window_frames = smoothing_window_frames + 1; 
    end
    fprintf('Smoothing movie data for player...');

    % Smooth the 5D data along the 5th (time) dimension (no 'squeeze', smoothdata(A,dim,method,window))
    smoothedResultData = smoothdata(result(:,:,:,channelToAnalyze,:), 5, 'movmean', smoothing_window_frames);
    fprintf(' Done.\n');

    movieData1 = zeros(height, width * 2, 1, numTime, 'like', result); % Original
    movieData2 = zeros(height, width * 2, 1, numTime, 'like', result); % Smoothed

    for t = 1:numTime
        % 1. Get original data (uint16)
        slice1_orig = result(:, :, 1, channelToAnalyze, t)';
        slice2_orig = result(:, :, 2, channelToAnalyze, t)';

        limits1 = stretchlim(slice1_orig);
        limits2 = stretchlim(slice2_orig);

        movieData1(:, :, 1, t) = [imadjust(slice1_orig, limits1), ...
            imadjust(slice2_orig, limits2)];

        % 2. Get smoothed data (double)
        slice1_smooth_double = smoothedResultData(:, :, 1, 1, t)';
        slice2_smooth_double = smoothedResultData(:, :, 2, 1, t)';

        % Convert to uint16 for manual contrast in player
        movieData2(:, :, 1, t) = [uint16(slice1_smooth_double), ...
            uint16(slice2_smooth_double)];
    end
    fprintf('Launching movie players for current trial. \n');

    hPlayer1 = implay(movieData1, volumeRate);
    hPlayer1.Parent.Name = sprintf('TRIAL %d: %s - ORIGINAL', trial_idx, currentFileName);

    hPlayer2 = implay(movieData2, volumeRate);
    hPlayer2.Parent.Name = sprintf('TRIAL %d: %s - SMOOTHED (%.1fs)', trial_idx, currentFileName, smoothing_window_sec);


    % Get the 5D data for the analysis channel.
    channelData = result(:, :, :, channelToAnalyze, :);

    % --- Define both baseline types F0 ---
    stimFrames = round(stimTimes_sec * volumeRate);
    rolling_F0_frames = round(rolling_F0_sec * volumeRate);
    fixed_F0_frames = [max(1, round(fixed_F0_sec(1) * volumeRate)), ...
        min(numTime, round(fixed_F0_sec(2) * volumeRate))];

    fprintf('Using fixed F0 for feedback plots: Frames %d to %d (%.2f to %.2f s)\n', ...
        fixed_F0_frames(1), fixed_F0_frames(2), ...
        fixed_F0_frames(1)/volumeRate, fixed_F0_frames(2)/volumeRate);


% --- 8) ROI DRAWING: CHECK SAVED -> CONFIRM -> OR DRAW NEW ---

    % Setup flags - Default state = manual drawing
    skip_manual_drawing = false;
    rois_are_correct = false; % controls manual loop


    % Average Projection for Visualisation
    avgTimeProjection = mean(double(smoothedResultData), 5);
    slice1_data = avgTimeProjection(:, :, 1, 1)';
    slice2_data = avgTimeProjection(:, :, 2, 1)';

    
    % --- Check if ROIs exist for this trial ---
    if trial_idx <= length(all_ROIs) && ~isempty(all_ROIs{trial_idx})
        savedMasks = all_ROIs{trial_idx};
        
        % Retrieve slice info BEFORE plotting
        if trial_idx <= length(all_roiSlices) && ~isempty(all_roiSlices{trial_idx})
            roiSlices = all_roiSlices{trial_idx};
        else
            % Fallback only if data is missing
            roiSlices = [1; 1; 2; 2]; 
            warning('No slice info found for saved ROIs. Defaulting to standard.');
        end

        % Validate dimensions
        if size(savedMasks{1}, 1) == height && size(savedMasks{1}, 2) == width
            
            % Ask User or Auto-Accept
            if auto_accept_saved
                answer = 'Yes';
            else
                % Show saved ROIs
                hCheck = figure('Name', sprintf('Verify Saved ROIs: %s', currentFileName), ...
                    'Position', [100 100 1000 500], 'Color', 'w');
                
                % Setup Axes
                ax_check1 = subplot(1,2,1); 
                imagesc(slice1_data); axis image; colormap parula; hold on; 
                title('Slice 1 (Saved)');
                
                ax_check2 = subplot(1,2,2); 
                imagesc(slice2_data); axis image; colormap parula; hold on; 
                title('Slice 2 (Saved)');
                
                % Dynamic Plotting based on roiSlices
                roiColors = {'m', 'c', 'm', 'c'};
                for r = 1:4
                    thisSlice = roiSlices(r);
                    if thisSlice == 1
                        plot_boundary(savedMasks{r}, roiColors{r}, ax_check1); 
                    elseif thisSlice == 2
                        plot_boundary(savedMasks{r}, roiColors{r}, ax_check2);
                    end
                end

                % Ask User (last option = repeat = default choice)
                answer = questdlg('Use these saved ROIs?', 'Found Existing Data', 'Yes', 'Yes to All', 'No', 'Yes');
                close(hCheck);
            end
            
            % Process answer
            if strcmp(answer, 'Yes') || strcmp(answer, 'Yes to All')

                if strcmp(answer, 'Yes to All')
                auto_accept_saved = true; % set to true when user selects 'Yes to All', bypasses confirmation dialogue for remaining trials
                                          % resets to false at start of each fly
                fprintf(' -> Auto-accept enabled for remaining trials.\n');
                end

                roiMasks = savedMasks;

                all_roiSlices{trial_idx} = roiSlices;

                % Set flags to skip manual loop
                skip_manual_drawing = true;
                rois_are_correct = true; % Skips the manual loop
               
                
                % --- Calculate feedback traces ---
                
                channelData = result(:,:,:, channelToAnalyze,:);
                raw_F_traces_feedback = zeros(numTime, 4);
                for i = 1:4
                    zSliceForROI = roiSlices(i);
                    currentMask = roiMasks{i};
                    pixelIndices = find(currentMask);
                    for t = 1:numTime
                        currentSlice = channelData(:, :, zSliceForROI, 1, t)'; % Transpose here before indexing data with mask in next line: 
                        % ROI masks are drawn on transposed images (height, width)
                        raw_F_traces_feedback(t, i) = mean(currentSlice(pixelIndices));
                    end
                end
                fprintf(' -> Accepted saved ROIs.\n');
            else
                fprintf('  -> User chose to redraw.\n');
            end
        end
    end
    

    % --- Manual drawing loop --- 
    % Only runs if skip_manual_drawing is FALSE (if no saved ROIs or user said "No, Redraw")

    while ~skip_manual_drawing && ~rois_are_correct       
        
        % Create figure with subplots - try to mimic FIJI image quality by
        % using custom color map, smart contrast and bilinear interpolation
        % (Matlab default is nearest-neighbour -> looks more grainy)

        % Define FIJI-style Fire Colormap
        % (Black -> Blue -> Magenta -> Orange -> Yellow -> White)
        keyColors = [0 0 0 ; 0 0 1; 1 0 1; 1 0.5 0; 1 1 0; 1 1 1];
        newMap = interp1(linspace(0, 1, 6), keyColors, linspace(0,1,256));
        
        % Specify window handle
        hFig = figure('Name', sprintf('TRIAL %d: %s', trial_idx, currentFileName), 'Position', [50 50 1200 600]);

        ax1 = subplot(1, 2, 1);
        % Smart contrast (clip top 0.2% outliers, flatten first to find percentiles)
        temp_vec1 = double(slice1_data(:));
        min_val1 = min(temp_vec1);
        max_val1 = prctile(temp_vec1, 99.8);
        
        % Specify actual image handle
        hIm1 = imagesc(slice1_data);
        hIm1.Interpolation = 'bilinear'; 

        axis image;
        colormap(ax1, newMap);
        clim(ax1, [min_val1, max_val1]); % set colormap limits - clip top 0.2%
        title('Slice 1');

        ax2 = subplot(1, 2, 2);
        % Smart contrast
        temp_vec2 = double(slice2_data(:));
        min_val2 = min(temp_vec2);
        max_val2 = prctile(temp_vec2, 99.8);
        
        hIm2 = imagesc(slice2_data);
        hIm2.Interpolation = 'bilinear';

        axis image; 
        colormap(ax2, newMap);
        clim(ax2, [min_val2, max_val2]);
        title('Slice 2');
        
        % For manually adjusting contrast
        imcontrast(ax1);
        imcontrast(ax2);

        roiMasks = cell(4,1);
        roiSlices = zeros(4,1);

        % Labels 
        roiLabels = {'1: BKG (Pair 1)', '2: SIG (Pair 1)', '3: BKG (Pair 2)', '4: SIG (Pair 2)'};
        roiColors = {'magenta', 'cyan', 'magenta', 'cyan'};

        graphic_handles = [];

        for i = 1:4
            % Ask for the slice once per pair (dialogue box appears on drawing ROIs 1 & 3)
            if mod(i, 2) == 1
                sliceChoice = questdlg(sprintf('Draw pair %d (ROIs %d & %d) on which slice?', (i+1)/2, i, i+1), ...
                    'Choose slice', 'Slice 1', 'Slice 2', 'Slice 1');
                switch sliceChoice
                    case 'Slice 1', targetAx = ax1; sliceNum = 1;
                    case 'Slice 2', targetAx = ax2; sliceNum = 2;
                    otherwise,      error('ROI drawing cancelled.');
                end
                figure(hFig);   % bring the drawing window back to the front
            end

            sgtitle(sprintf('ZOOM/PAN/CONTRAST, then press any key to draw ROI %d (%s)', i, roiLabels{i}));
            waitforbuttonpress;

            % Switch off any zoom/pan tool left on, so drawpolygon can take over mouse and keyboard
            zoom(hFig, 'off');
            pan(hFig, 'off');
            
            sgtitle(['DRAW ROI ' roiLabels{i}]);
            h_roi = drawpolygon(targetAx, 'Color', roiColors{i}, ...
                'LineWidth', 1, ...
                'FaceAlpha', 0.0); % No fill
            if isempty(h_roi.Position), error('ROI drawing cancelled.'); end

            roiSlices(i) = sliceNum;
            roiMasks{i}  = h_roi.createMask();

            pos = mean(h_roi.Position, 1);
            h_text = text(h_roi.Parent, pos(1), pos(2), num2str(i), 'Color', 'w', 'FontSize', 16, 'FontWeight', 'bold');
            graphic_handles = [graphic_handles, h_roi, h_text];
        end
        
        % Safety check to ensure background and signal are drawn on same
        % slice
        if roiSlices(1) ~= roiSlices(2)
             warndlg(sprintf('WARNING: Background 1 (Slice %d) and Signal 1 (Slice %d) are on different slices!', ...
                 roiSlices(1), roiSlices(2)), 'Slice Mismatch Warning');
        end
        if roiSlices(3) ~= roiSlices(4)
             warndlg(sprintf('WARNING: Background 2 (Slice %d) and Signal 2 (Slice %d) are on different slices!', ...
                 roiSlices(3), roiSlices(4)), 'Slice Mismatch Warning');
        end
        
        
   % --- Immediate Feedback Loop ---
        % For immediate ROI trace inspection. Either accept ROIs or redraw

        sgtitle('All ROIs drawn. Calculating feedback traces...');

        % Extract raw traces
        raw_F_traces_feedback = zeros(numTime, 4);
        for i = 1:4
            zSliceForROI = roiSlices(i); % saves which slice ROI was drawn on
            currentMask = roiMasks{i};
            pixelIndices = find(currentMask);
            for t = 1:numTime
                currentSlice = channelData(:, :, zSliceForROI, 1, t)';
                raw_F_traces_feedback(t, i) = mean(currentSlice(pixelIndices));
            end
        end

        % Save individual raw ROI traces
        F_raw_bkg1 = raw_F_traces_feedback(:, 1);
        F_raw_sig1 = raw_F_traces_feedback(:, 2);
        F_raw_bkg2 = raw_F_traces_feedback(:, 3);
        F_raw_sig2 = raw_F_traces_feedback(:, 4);

        % Subtract the raw background from the raw signal
        F_raw_corrected_s1 = F_raw_sig1 - F_raw_bkg1;
        F_raw_corrected_s2 = F_raw_sig2 - F_raw_bkg2;

        % Calculate all 4 dF/F traces for plotting - using fixed F0 set in
        % GUI at beginning - Pass the raw traces and smooth_window to the helper function (calculate_fixed_dff)
        dff_bkg1 = calculate_fixed_dff(F_raw_bkg1, smoothing_window_frames, fixed_F0_frames);
        dff_sig1_corrected = calculate_fixed_dff(F_raw_corrected_s1, smoothing_window_frames, fixed_F0_frames);
        dff_bkg2 = calculate_fixed_dff(F_raw_bkg2, smoothing_window_frames, fixed_F0_frames);
        dff_sig2_corrected = calculate_fixed_dff(F_raw_corrected_s2, smoothing_window_frames, fixed_F0_frames);

        % Plot the final traces
        dF_F0_feedback = nan(numTime, 4); % Pre-allocate

        dF_F0_feedback(:, 1) = dff_bkg1; % Plot Background 1
        dF_F0_feedback(:, 2) = dff_sig1_corrected; % Plot Signal 1 (Corrected)
        dF_F0_feedback(:, 3) = dff_bkg2; % Plot Background 2
        dF_F0_feedback(:, 4) = dff_sig2_corrected; % Plot Signal 2 (Corrected)

        % Plot feedback
        hFeedbackFig = figure('Name', 'ROI Feedback Traces');
        trace_titles = {'Background 1', 'Signal 1 (Corrected)', 'Background 2', 'Signal 2 (Corrected)'};

        for i = 1:4
            subplot(2, 2, i);
            plot(dF_F0_feedback(:, i), 'Color', roiColors{i}, 'LineWidth', 1);
            hold on;
            for j = 1:numel(stimFrames)
                xline(stimFrames(j), '--k');
            end
            hold off;
            title(trace_titles{i});
            ylabel('\DeltaF/F_0');
            xlabel('Time (frames)');
            axis tight; grid on;
        end

        % Prompt: Are you happy with traces?
        sgtitle('Feedback: dF/F Traces (Fixed F0)');
        figure(hFig);
        answer = questdlg('Are you happy with these ROIs?', 'Confirm ROIs', 'Yes', 'No, Redraw', 'Yes');
        
        if strcmp(answer, 'Yes')
            rois_are_correct = true;
            all_ROIs{trial_idx} = roiMasks;
            close(hFeedbackFig);
            close(hFig); % Close ROI drawing window
        else
            delete(graphic_handles);
            close(hFig);
            close(hFeedbackFig);
            disp('ROIs cleared. Please redraw.');
        end

        all_roiSlices{trial_idx} = roiSlices;
    end % End of ROI Drawing loop


    % --- 9) SHOT NOISE VALIDATION  ---
    % Runs split-half Monte Carlo Resampling (No replacement, unlike bootstrapping)
    % (Uses Rolling Baseline, 1000 resamples per trial) on transposed image
    
    % Resampling is done on RAW traces
    % Then uses rolling dF/F per half to mimic the actual data processing
    
    slice_idx_1 = roiSlices(2); % Signal 1
    slice_idx_2 = roiSlices(4); % Signal 2
    
    % Extract slice data
    temp_S1 = permute(reshape(channelData(:, :, slice_idx_1, 1, :), width, height, numTime), [2 1 3]); 
    temp_S2 = permute(reshape(channelData(:, :, slice_idx_2, 1, :), width, height,  numTime), [2 1 3]);
    
    data_S1_2D = reshape(temp_S1, [], numTime); % pixels x time
    data_S2_2D = reshape(temp_S2, [], numTime);
    
    % ROI pixels
    pixels_S1 = find(roiMasks{2});
    pixels_S2 = find(roiMasks{4});
    n_pix_1 = numel(pixels_S1);
    n_pix_2 = numel(pixels_S2);
    
    % Background traces (raw)
    F_bkg_s1_raw = raw_F_traces_feedback(:, 1);
    F_bkg_s2_raw = raw_F_traces_feedback(:, 3);
    
    response_window_frames_roi1 = round(response_window_sec_roi1 * volumeRate);
    response_window_frames_roi2 = round(response_window_sec_roi2 * volumeRate);
    
    n_iters = 1000;
    hWait = waitbar(0,sprintf('Running %d resamples...', n_iters));
    
    % Preallocate for peaks
    stim_maxes_A1 = nan(n_iters, numel(stimFrames));
    stim_maxes_B1 = nan(n_iters, numel(stimFrames));
    stim_maxes_A2 = nan(n_iters, numel(stimFrames));
    stim_maxes_B2 = nan(n_iters, numel(stimFrames));
    
    for k = 1:n_iters
        if mod(k,50) == 0
            waitbar(k/n_iters, hWait);
        end
    
        % --- Signal 1 ROI split ---
        idx_perm1 = randperm(n_pix_1);
        halfA1 = pixels_S1(idx_perm1(1:floor(n_pix_1/2)));
        halfB1 = pixels_S1(idx_perm1(floor(n_pix_1/2)+1:end));
    
        trace_A1_raw = mean(data_S1_2D(halfA1, :), 1)' - F_bkg_s1_raw;
        trace_B1_raw = mean(data_S1_2D(halfB1, :), 1)' - F_bkg_s1_raw;
    
        dF_A1 = calculate_rolling_dff(trace_A1_raw, smoothing_window_frames, stimFrames, rolling_F0_frames, numTime);
        dF_B1 = calculate_rolling_dff(trace_B1_raw, smoothing_window_frames, stimFrames, rolling_F0_frames, numTime);
    
        for i = 1:numel(stimFrames)
            startFrame = stimFrames(i);
            endFrame = min(numTime, startFrame + response_window_frames_roi1 - 1);
            stim_maxes_A1(k,i) = max(dF_A1(startFrame:endFrame));
            stim_maxes_B1(k,i) = max(dF_B1(startFrame:endFrame));
        end
    
        % --- Signal 2 ROI split ---
        idx_perm2 = randperm(n_pix_2);
        halfA2 = pixels_S2(idx_perm2(1:floor(n_pix_2/2)));
        halfB2 = pixels_S2(idx_perm2(floor(n_pix_2/2)+1:end));
    
        trace_A2_raw = mean(data_S2_2D(halfA2, :), 1)' - F_bkg_s2_raw;
        trace_B2_raw = mean(data_S2_2D(halfB2, :), 1)' - F_bkg_s2_raw;
    
        dF_A2 = calculate_rolling_dff(trace_A2_raw, smoothing_window_frames, stimFrames, rolling_F0_frames, numTime);
        dF_B2 = calculate_rolling_dff(trace_B2_raw, smoothing_window_frames, stimFrames, rolling_F0_frames, numTime);
    
        for i = 1:numel(stimFrames)
            startFrame = stimFrames(i);
            endFrame = min(numTime, startFrame + response_window_frames_roi2 - 1);
            stim_maxes_A2(k,i) = max(dF_A2(startFrame:endFrame));
            stim_maxes_B2(k,i) = max(dF_B2(startFrame:endFrame));
        end
    end
    
    close(hWait);
    fprintf(' Shot noise resampling done.\n');
    
    % --- POOL SHOT-NOISE SAMPLES ---
    % Accumulate into 3D arrays [5 stimuli x 1000 iterations x nTrials]
    % cat(3, ...) accumulates trials along 3rd dimension
    % Flattening to 1D happens later in the analysis script
    all_shotnoise_A1 = cat(3, all_shotnoise_A1, permute(stim_maxes_A1, [2 1]));
    all_shotnoise_B1 = cat(3, all_shotnoise_B1, permute(stim_maxes_B1, [2 1]));
    all_shotnoise_A2 = cat(3, all_shotnoise_A2, permute(stim_maxes_A2, [2 1]));
    all_shotnoise_B2 = cat(3, all_shotnoise_B2, permute(stim_maxes_B2, [2 1]));

   
    % --- 10) EXTRACT MBON RESPONSES ----
    % Uses rolling dF/F and smoothing window, as set in GUI
    % Extracts peak responses per stimulus window and saves them in 5x10 matrix for
    % each fly and ROI

    % Get the raw traces
    raw_bkg1 = raw_F_traces_feedback(:, 1);
    raw_sig1 = raw_F_traces_feedback(:, 2);
    raw_bkg2 = raw_F_traces_feedback(:, 3);
    raw_sig2 = raw_F_traces_feedback(:, 4);

    % Subtract background
    raw_corrected_s1 = raw_sig1 - raw_bkg1;
    raw_corrected_s2 = raw_sig2 - raw_bkg2;

    % Calculate dF/F - rolling F0
    dF_F0_traces = nan(numTime, 2); % Pre-allocate
    % ROI1
    dF_F0_traces(:, 1) = calculate_rolling_dff(raw_corrected_s1, smoothing_window_frames, stimFrames, rolling_F0_frames, numTime);
    % ROI2
    dF_F0_traces(:, 2) = calculate_rolling_dff(raw_corrected_s2, smoothing_window_frames, stimFrames, rolling_F0_frames, numTime);

    % Pre-allocate
    responseAmplitudes_ROI1 = nan(numel(stimFrames), 1);
    responseAmplitudes_ROI2 = nan(numel(stimFrames), 1);

    % Find peak response within user-defined response window
    for i = 1:numel(stimFrames)
        % ROI1
        startFrame1 = stimFrames(i); 
        endFrame1 = min(numTime, startFrame1 + response_window_frames_roi1 - 1); % min(numTime... prevents indexing beyond data
        [max_val, ~] = max(dF_F0_traces(startFrame1:endFrame1, 1));
        responseAmplitudes_ROI1(i) = max_val;
        % ROI2
        startFrame2 = stimFrames(i); 
        endFrame2 = min(numTime, startFrame2 + response_window_frames_roi2 - 1);
        [max_val, ~] = max(dF_F0_traces(startFrame2:endFrame2, 2));
        responseAmplitudes_ROI2(i) = max_val;
    end


    % Save peak MBON responses
    all_amplitudes_ROI1 = [all_amplitudes_ROI1, responseAmplitudes_ROI1]; % concatenate
    all_amplitudes_ROI2 = [all_amplitudes_ROI2, responseAmplitudes_ROI2];


    % Added padding as for one fly, some trials were 5 frames shorter
    % No safety limit of #frames added yet, as this was isolated case and limited to
    % 5 frames. 
    % NB substantial differences in length  will distort mean traces and should be handled differently
    if isempty(all_traces_ROI1)
        all_traces_ROI1 = dF_F0_traces(:, 1);
        all_traces_ROI2 = dF_F0_traces(:, 2);
    else
        len_stored = size(all_traces_ROI1, 1);
        len_new    = size(dF_F0_traces, 1);
        
        if len_new < len_stored
            % New trial is SHORTER: Pad the NEW trial with NaNs at the end
            pad_len = len_stored - len_new;
            fprintf('  Warning: Trial %d is shorter (-%d frames). Padding with NaNs.\n', ...
                trial_idx, pad_len);
            
            nan_pad = nan(pad_len, 1);
            
            % Add NaNs to bottom of new vectors
            col1_padded = [dF_F0_traces(:, 1); nan_pad];
            col2_padded = [dF_F0_traces(:, 2); nan_pad];
            
            % Append
            all_traces_ROI1 = [all_traces_ROI1, col1_padded];
            all_traces_ROI2 = [all_traces_ROI2, col2_padded];
            
        elseif len_new > len_stored
            % New trial is LONGER: Pad the EXISTING matrix with NaNs
            pad_len = len_new - len_stored;
            n_stored_trials = size(all_traces_ROI1, 2);
            fprintf('  Warning: Trial %d is longer (+%d frames). Padding previous trials with NaNs.\n', ...
                trial_idx, pad_len);
            
            nan_mat = nan(pad_len, n_stored_trials);
            
            % Add NaNs to bottom of existing matrices
            all_traces_ROI1 = [all_traces_ROI1; nan_mat];
            all_traces_ROI2 = [all_traces_ROI2; nan_mat];
            
            % Append the new long trial
            all_traces_ROI1 = [all_traces_ROI1, dF_F0_traces(:, 1)];
            all_traces_ROI2 = [all_traces_ROI2, dF_F0_traces(:, 2)];
            
        else
            % Exact match
            all_traces_ROI1 = [all_traces_ROI1, dF_F0_traces(:, 1)];
            all_traces_ROI2 = [all_traces_ROI2, dF_F0_traces(:, 2)];
        end
    end


    %% -------------------------------------------------------------------------------
    %               3. SAVE PROGRESS & OUTPUT
    % ---------------------------------------------------------------------------------
    
    % ---  11) SAVE PROGRESS

    % Allows to resume data analysis at later point
    fprintf('\n Saving Progress to %s \n', saveFilePath);
    save(saveFilePath, 'all_amplitudes_ROI1',...
        'all_amplitudes_ROI2', 'all_traces_ROI1', 'all_traces_ROI2', 'all_SI', ...
        'all_ROIs', 'all_roiSlices', 'all_shotnoise_A1', 'all_shotnoise_B1', ...
        'all_shotnoise_A2', 'all_shotnoise_B2', 'stimFrames');

    
    % ---  12) SAVE OUTPUT TO MASTER FOLDER
    
    fprintf(' Saving Output to Master Folder: %s \n', detailedDataFile);
    save(detailedDataFile, 'flyID', 'numTrials', 'actualTrialNumbers','stimFrames', ...
        'all_amplitudes_ROI1', 'all_amplitudes_ROI2','all_traces_ROI1', ...
        'all_traces_ROI2', 'all_SI', 'all_ROIs', 'all_roiSlices', 'all_shotnoise_A1',...
        'all_shotnoise_B1','all_shotnoise_A2', 'all_shotnoise_B2', 'stimFrames',...
        'n_pix_1', 'n_pix_2');

end % END MAIN LOOP THROUGH TRIALS

fprintf('\n All trials extracted. Move on to next fly or analyse data with MultNoise_StatsAndPlotting_Deming_ErrorVariance.m script.\n');


%% -------------------------------------------------------------------------------
%              4. HELPER FUNCTIONS
% ---------------------------------------------------------------------------------

% --- DRAW ROIS OUTLINES FROM SAVED ---
function plot_boundary(mask, color, axHandle)
    if nargin < 3, axHandle = gca; 
    end % Default to current axes
    [B, ~] = bwboundaries(mask);
    for k = 1:length(B)
        boundary = B{k};
        plot(axHandle, boundary(:,2), boundary(:,1), 'Color', color, 'LineWidth', 2);
    end
end


% --- ROLLING DF/F ---
% This function takes a RAW trace, smoothes it, and calculates
% its rolling-baseline dF/F. It does NO background subtraction

function dF_F_trace = calculate_rolling_dff(trace_raw, smooth_window, stimFrames, rolling_F0_frames, numTime)

    % Smooth the raw trace
    trace_smooth = movmean(trace_raw, smooth_window);
    
    % Calculate rolling baseline dF/F
    % Baseline window immediately precedes each stimulus onset.
    % Each epoch runs from stimulus i onset to stimulus i+1 onset, or end
    % of recording
    dF_F_trace = nan(numTime, 1);
    for i = 1:numel(stimFrames)
        baselineEndFrame = stimFrames(i) - 1;
        baselineStartFrame = max(1, baselineEndFrame - rolling_F0_frames + 1);

        % To include first 5 seconds, in all_traces_ROI. This is needed for plot_mean_traces function in 
        % MultNoise_StatsAndPlotting_Deming_ErrorVariance.m analysis script
    
        if i == 1 % only for first stimulus
            epochStartFrame = 1;
        else
            epochStartFrame = stimFrames(i);  % previoulsy only this, start populating trace at frame of first stimulus (58) 1-57 was NaNs
        end

        if i < numel(stimFrames)
            epochEndFrame = stimFrames(i+1) - 1;
        else
            epochEndFrame = numTime;
        end
    
        if baselineStartFrame > baselineEndFrame
            F0 = trace_smooth(baselineStartFrame);
        else
            F0 = mean(trace_smooth(baselineStartFrame:baselineEndFrame));
        end
    
        F_epoch = trace_smooth(epochStartFrame:epochEndFrame);
    
        % Check for unstable baseline. Near-zero F0 can occur. To avoid
        % division by near-zero, fall back to subtraction
        if F0 ~= 0 && abs(F0) > 1e-5 
            dF_F_trace(epochStartFrame:epochEndFrame) = (F_epoch - F0) / F0;
        else
            % If baseline is zero, just subtract it (F-F0)
            dF_F_trace(epochStartFrame:epochEndFrame) = F_epoch - F0;
        end
    end
end


% --- FIXED DF/F ---
% This helper function first smoothes the trace, then calculates dF/F
% using a single fixed baseline as set by user in GUI

function dF_F_trace = calculate_fixed_dff(trace_raw, smooth_window, fixed_F0_frames)

    trace_smooth = movmean(trace_raw, smooth_window);
    
    F0_val = mean(trace_smooth(fixed_F0_frames(1):fixed_F0_frames(2)));
    if F0_val ~= 0 && abs(F0_val) > 1e-5
        dF_F_trace = (trace_smooth - F0_val) / F0_val;
    else
        dF_F_trace = trace_smooth - F0_val;
    end
end
