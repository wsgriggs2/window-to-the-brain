%% human_glm_task.m
% This script performs the basic GLM analysis for the human experiments.
% The human experiments consisted of blocks of rest and drawing or blocks
% of playing guitar and rest.

%% User settings
clear; close all; clc

% Define p value threshold for parametric statistics
p_threshold = 1e-3;

% Load some colormaps
lateralizationMapComplete = load('LateralizationP.mat');
PositiveOnlyColormap = lateralizationMapComplete.LateralizationP;
PositiveOnlyColormap = flipud(PositiveOnlyColormap(1:floor(length(PositiveOnlyColormap(:,1))/2),:));

PositiveNegativeColormap = load('PercentChangeColormap.mat');
PositiveNegativeColormap = PositiveNegativeColormap.PercentChangeColormap;

%% Load data
[file_to_load, path_to_file] = uigetfile(w2b.io.getUserDataPath, 'Choose data file to load');
data = load(fullfile(path_to_file, file_to_load), ...
    'dop', 'angiogram', 'task', 'timestamps', 'UF', 'run_label');

% Define pixel ratio - they are non isotropic
x_pixelsize = 0.3;
z_pixelsize = data.UF.Lambda;
data_aspect_ratio = [ 1 x_pixelsize/z_pixelsize 1];


%% Perform motion correction
Dop_corrected = w2b.preproc.correctMotion(data.dop);
angiogram = w2b.util.makeAngiogram(Dop_corrected);

% Generate axe labels for future plots
X_img_mm = x_pixelsize/2 + (0:size(Dop_corrected,2)-1)*x_pixelsize;
Z_img_mm = z_pixelsize/2 + (0:size(Dop_corrected,1)-1)*z_pixelsize;

%% Preprocessing
% Get size of the doppler data
[nDepth, nWidth, nTimepoints] = size(Dop_corrected);

% Apply spatial filter to entire session.
Dop_blur = w2b.preproc.spatialFilter(Dop_corrected, {'gaussian', 5, 1});

% View the spatially smoothed data
figure;
imagesc(w2b.util.makeAngiogram(Dop_blur));
daspect(data_aspect_ratio);

% Apply temporal smoothing to reduce noise in the data
% Moving average filter
windowSize = 5;
if ~isnan(windowSize)
    Dop_proc = movmean(Dop_blur, [windowSize 0], 3);
else
    Dop_proc = Dop_blur;
end

% Detrend
window_length = 300;
Dop_proc = w2b.preproc.detrendSlidingWindow(Dop_proc, window_length);

% Two ways to calculate percent change
method_to_use = 'baseline';
if strcmp(method_to_use, 'baseline')
    % Scale doppler to percent change from rest blocks
    Dop_scaled = w2b.preproc.scaleDoppler(Dop_proc, 'baseline_index', ~logical(data.task));

elseif strcmp(method_to_use, 'mean activity')
    % Scale doppler to percent change from mean voxel activity
    Dop_scaled = w2b.preproc.scaleDoppler(Dop_proc);
end


%% Behavioral variables
task_start_ind = find(diff(data.task) > 0) + 1;
task_end_ind = [find(diff(data.task) < 0) size(Dop_corrected,3)];

%% Create GLM regressors
% Convert into ms
fUS_timestamps = data.timestamps(:, 1);
fUS_startTime = min(fUS_timestamps(:, 1)) * 1000;
fUS_endTime = max(fUS_timestamps(:, 1)) * 1000;

% Task will be modeled as boxcar
boxcar_regressor_names = {'task on'};
boxcar_event_timing = {fUS_timestamps(task_start_ind)};
boxcar_event_duration = {fUS_timestamps(task_end_ind) - fUS_timestamps(task_start_ind)};

% If there were impulse (delta functions), could add those here.
% No impulses, so just boxcar.
combined_event_labels = boxcar_regressor_names;
combined_event_timing = boxcar_event_timing;
combined_event_duration = boxcar_event_duration;
combined_event_modulation = cell(1, 1);


% Remove all NaNs from the regressors
nRegressors = length(combined_event_labels);
for regressor = 1:nRegressors
    % Pull out info about current regressor
    regressor_timing = combined_event_timing{regressor};
    regressor_duration = combined_event_duration{regressor};
    nTrials = length(regressor_timing);

    % Find NaNs
    nan_ind = isnan(regressor_timing);
    if ~isempty(regressor_duration)
        nan_ind = nan_ind | isnan(regressor_duration);
    end

    % Remove NaNs
    combined_event_timing{regressor} = regressor_timing(~nan_ind);
    if ~isempty(regressor_duration)
        combined_event_duration{regressor} = regressor_duration(~nan_ind);
    end
end

% Combine into table
regressor_table = table(combined_event_labels', ...
    combined_event_timing, ...
    combined_event_duration, ...
    combined_event_modulation,...
    'VariableNames', {'regressorName', 'timestamps', 'duration', 'modulation'});

% Create HRF
n_secs = 16; % In secs; length of desired HRF
hrf = w2b.stats.generateHRF(n_secs, 'tau', 0.7, 'delta', 3, 'n', 3, 'Fs', 1000);

%% Create regressors
% Create regressors of interest
regressors = w2b.stats.createGLMRegressors(fUS_timestamps(:, 1), regressor_table, hrf.shape);

% Add constant regressor
regressors = [regressors ones(size(regressors, 1), 1)];

% Select regressors of interest
regressors_of_interest = regressors;
regressors_of_interest_ind = 1:2;


%% Build GLM
y = Dop_scaled;
patch_size = 1; % Size of patches in pixels/voxels
[pvalue, tscore, beta_hat] = w2b.stats.voxelwiseGLM(...
    y, regressors_of_interest, patch_size);


%% Apply specific contrast
% For each regressor, apply the T contrast
[pvalue_cont, tscore_cont] = deal(zeros(size(y, 1), size(y, 2), 2));
for i = 1:2
    cont = zeros(1, height(regressor_table)+1);
    cont(i) = 1;

    if size(cont, 2) ~= size(regressors, 2)
        error('The contrast must have same number of elements as there are regressors');
    end
    [pvalue_cont(:, :, i), tscore_cont(:, :, i)] = w2b.stats.Tcontrast(...
        y, regressors, beta_hat, cont, patch_size, ...
        'multiple_comparison_correction', 'FDR');
end


%% Threshold image
p_threshold = 1e-10;
beta_hat_plot = beta_hat;

pvalue_thresholded = pvalue_cont;
tscore_thresholded = tscore_cont;
beta_hat_thresholded = beta_hat_plot;
threshold_mask = pvalue_cont>p_threshold;
pvalue_thresholded(threshold_mask) = NaN;
tscore_thresholded(threshold_mask) = NaN;
beta_hat_thresholded(threshold_mask) = NaN;


%% Plot T-scored map
tValuesToDisplay = [-14 14];

figure(1);
% Generate the colormap of the t-scored activation to each direction
cmap = w2b.plot.plotDuplexImage(X_img_mm, Z_img_mm, tscore_thresholded(:, :, 1),angiogram,...
    'colormap2use', PositiveNegativeColormap, 'nonlinear_bg', 2, ...
    'AutoColorBarLimits', tValuesToDisplay, ...
    'showColorbar', true, ...
    'ColorBarTitle', 't-score', ...
    'data_aspect_ratio', data_aspect_ratio);
title(sprintf('T-contrast map - %s\n', data.run_label));


%% Plot voxel across time
window_size = 3; % In voxels

% Since the voxels are not isotropic, we need to define separate sizes for
% x and z
window_size_x = window_size;
window_size_z = window_size*.3/data.UF.Lambda;
window_size_mm_x = window_size_x * x_pixelsize;
window_size_mm_z = window_size_z * z_pixelsize;

figure(1);
if exist('roi_point', 'var')
    delete(roi_point);
end
roi_point = drawpoint;
if exist('rec', 'var')
    delete(rec)
end

voxel_of_interest = roi_point.Position;
ROI_x_mm = voxel_of_interest(1)-window_size_mm_x:voxel_of_interest(1)+window_size_mm_x;
ROI_z_mm = voxel_of_interest(2)-window_size_mm_z:voxel_of_interest(2)+window_size_mm_z;
% Display ROI on plot
rec = rectangle('Position', [min(ROI_x_mm), min(ROI_z_mm), 2*window_size_mm_x, 2*window_size_mm_z], 'EdgeColor', 'w', 'LineWidth', 3);

% Convert from mm to pixel space
[~, Xidx] = min(abs(X_img_mm-voxel_of_interest(1)));
[~, Zidx] = min(abs(Z_img_mm-voxel_of_interest(2)));

% Define center of the patch of interest in pixel coordinates
pixelOfInterest = [Zidx Xidx];
ROI_x = round(pixelOfInterest(2)-window_size_x):round(pixelOfInterest(2)+window_size_x);
ROI_z = round(pixelOfInterest(1)-window_size_z):round(pixelOfInterest(1)+window_size_z);

% If selected ROI extends beyond plot, then confine to boundaries of
% existing plot.
if any(ROI_x <= 0)
    ROI_x(ROI_x <= 0) = 1;
end
if any(ROI_z <= 0)
    ROI_z(ROI_z <= 0) = 1;
end
if any(ROI_x > size(y, 2))
    ROI_x(ROI_x > size(y,2)) = size(y, 2);
end
if any(ROI_z > size(y, 1))
    ROI_z(ROI_z > size(y, 1)) = size(y, 1);
end


mean_activity = squeeze(mean(mean(y(ROI_z, ROI_x, :), 2), 1));

figure;
tiledlayout(1,2, 'TileSpacing', 'tight');
nexttile();
max_abs_value = max(abs(tscore_thresholded(:, :, 1)), [], 'all');
cmap = w2b.plot.plotDuplexImage(X_img_mm, Z_img_mm, tscore_thresholded(:, :, 1), angiogram,...
    'colormap2use', PositiveNegativeColormap, 'nonlinear_bg', 2, ...
    'showColorbar', true, ...
    'AutoColorBarLimits', [-max_abs_value max_abs_value], ...
    'ColorBarTitle', 't-score', ...
    'data_aspect_ratio', data_aspect_ratio);
title(sprintf('%s', combined_event_labels{regressors_of_interest_ind(1)}));

%Put rectangle onto the plot
ROI_x_mm = (voxel_of_interest(1)-window_size_mm_x:voxel_of_interest(1)+window_size_mm_x);
ROI_z_mm = (voxel_of_interest(2)-window_size_mm_z:voxel_of_interest(2)+window_size_mm_z);
% Display ROI on plot
rec2 = rectangle('Position', [min(ROI_x_mm), min(ROI_z_mm), 2*window_size_mm_x, 2*window_size_mm_z], 'EdgeColor', 'w', 'LineWidth', 3);
drawnow;

% Plot mean activity from ROI across time
nexttile()
task2_scaled = data.task * (max(mean_activity) - min(mean_activity));
plot(fUS_timestamps(:, 1) - fUS_startTime/1000, [mean_activity, 100+task2_scaled']-100);
xlabel('Time (s)');
ylabel('Scaled fUSI signal');
title('Timecourse of fUSI ROI with GLM fit overlay');
legend('mean Dop activity', 'task');


%% Average activity during rest and drawing blocks for the ROI

% Calculate roi statistics
activity_task = cell(size(voxel_of_interest, 1), 1);
activity_rest = cell(size(voxel_of_interest, 1), 1);
fprintf('\n--------Statistical results --------\n');
for roi = 1:size(voxel_of_interest, 1)

    % Convert from mm to pixel space
    [~, Xidx] = min(abs(X_img_mm-voxel_of_interest(roi, 1)));
    [~, Zidx] = min(abs(Z_img_mm-voxel_of_interest(roi, 2)));

    % Define center of the patch of interest in pixel coordinates
    pixelOfInterest = [Zidx Xidx];
    ROI_x = round(pixelOfInterest(2)-window_size_x):round(pixelOfInterest(2)+window_size_x);
    ROI_z = round(pixelOfInterest(1)-window_size_z):round(pixelOfInterest(1)+window_size_z);

    % Confine to plot boundaries
    if any(ROI_x <= 0)
        ROI_x(ROI_x <= 0) = 1;
    end
    if any(ROI_z <= 0)
        ROI_z(ROI_z <= 0) = 1;
    end
    if any(ROI_x > size(y, 2))
        ROI_x(ROI_x > size(y,2)) = size(y, 2);
    end
    if any(ROI_z > size(y, 1))
        ROI_z(ROI_z > size(y, 1)) = size(y, 1);
    end

    task_ind = logical(data.task);
    rest_ind = ~logical(data.task);
    activity_task{roi} = squeeze(mean(mean(y(ROI_z, ROI_x, task_ind), 2), 1));
    activity_rest{roi} = squeeze(mean(mean(y(ROI_z, ROI_x, rest_ind), 2), 1));

    mean_difference = mean(activity_task{roi}) - mean(activity_rest{roi});

    fprintf('ROI %d \nmean difference - %0.3f%% \n', roi,mean_difference);

end



