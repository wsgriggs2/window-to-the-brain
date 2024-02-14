%% rodent_glm_task.m
% This script performs the basic GLM analysis for the rodent experiments.
% The rodent experiments consisted of blocks of visual stimulation on/off.

%% Initialize workspace
clear; close all; clc;

% Load some colormaps
lateralizationMapComplete = load('LateralizationP.mat');
PositiveOnlyColormap = lateralizationMapComplete.LateralizationP;
PositiveOnlyColormap = flipud(PositiveOnlyColormap(1:floor(length(PositiveOnlyColormap(:,1))/2),:));

PositiveNegativeColormap = load('PercentChangeColormap.mat');
PositiveNegativeColormap = PositiveNegativeColormap.PercentChangeColormap;


%% Load data and extract relevant variables
% This is the metadata about the rodent experimental sessions and runs.
project_record_filename = 'rodent_session_record.json';

% Load rodent data
data = w2b.io.loadRodentData('project_record_filename', project_record_filename);

% Extract variables from data struct
Dop = data.dop;
UF = data.UF;
timestamps = data.timestamps;
session_name = data.session;
run_name = data.run;

% The fUS sequence acquires non-isotropic voxels, so for visualization,
% let's take that into account.
x_pixelsize = 0.3;
z_pixelsize = UF.Lambda;
data_aspect_ratio = [1 x_pixelsize/z_pixelsize 1];


%% Extract metadata about the requested session/run combination
% Load ProjectRecord which has all the metadata.
ProjectRecord = w2b.io.loadJSONAsTable(project_record_filename);

% Find block structure for the loaded data run
indx = strcmp(ProjectRecord.run_folder, run_name);
on_blocks = ProjectRecord.on_block_cell{indx};
off_blocks = ProjectRecord.off_block_cell{indx};

% Remove the first index so that the lengths of on and off start inds match
if off_blocks(1) == 1
    off_blocks(1) = [];
elseif on_blocks(1) == 1
    on_blocks(1) = [];
end

% Cover possibility that we end on different task state than we started
add_entry_to_task_end_ind = false;
if length(on_blocks) ~= length(off_blocks)
    % Add entry at end since this means that we ended on different state
    % than we started
    if length(on_blocks) > length(off_blocks)
        off_blocks = [off_blocks; size(Dop, 3)];
        add_entry_to_task_end_ind = true;
    else
        on_blocks = [on_blocks; size(Dop, 3)];
    end
end

% Specify task structure
task = zeros(1, size(Dop, 3));
for block = 1:size(on_blocks)
    task(on_blocks(block):off_blocks(block)) = 1;
end

% Find indices of when the task starts and ends
task_start_ind = find(diff(task) > 0) + 1;
task_end_ind = find(diff(task) < 0);
if add_entry_to_task_end_ind
    task_end_ind = [task_end_ind size(Dop, 3)];
end


%% Perform motion correction
Dop_corrected = w2b.preproc.correctMotion(Dop);


%% Crop image size to just brain
% Speeds up computations

% Display initial angiogram
figure;
angiogram = w2b.util.makeAngiogram(Dop_corrected);
w2b.plot.plotAngiogram(angiogram, 1:size(angiogram,2), 1:size(angiogram, 1),...
    'colormap', 'inferno', ...
    'colorbar_title', 'Intensity (AU)', ...
    'data_aspect_ratio', data_aspect_ratio);
title('Draw a rectangle around the brain');

% Draw rectangle on angiogram to crop to just the brain
fprintf('Draw a rectangle around the brain.\n');
rect_roi = drawrectangle;
fprintf('Press enter to accept the drawn rectangular ROI.\n');
pause;

% Restrict the image to the defined rectangular ROI
rect_roi = round(rect_roi.Position);
Dop_cropped = Dop(rect_roi(2):rect_roi(2)+rect_roi(4)-1, rect_roi(1):rect_roi(1)+rect_roi(3)-1, :);

% Display the new cropped angiogram
angiogram = w2b.util.makeAngiogram(Dop_cropped);
figure;
X_img_mm = x_pixelsize/2 + (0:size(angiogram,2)-1)*x_pixelsize;
Z_img_mm = z_pixelsize/2 + (0:size(angiogram,1)-1)*z_pixelsize;
w2b.plot.plotAngiogram(angiogram, X_img_mm, Z_img_mm,...
    'colormap', 'inferno', ...
    'colorbar_title', 'Intensity (AU)', ...
    'data_aspect_ratio', data_aspect_ratio);


%% Preprocess the Doppler data
[nDepth, nWidth, nTimepoints] = size(Dop_cropped);

% Define desired spatial filter
FWHM = 2; % In voxels; 
sigma = FWHM/(sqrt(8*log(2)));
filter_size = 2*ceil(2*sigma)+1; % This is the default filter size used in `imgaussfilt`.

% Apply spatial filter to entire session.
Dop_blur = w2b.preproc.spatialFilter(Dop_cropped, {'gaussian', filter_size, sigma});

% Apply smoothing to reduce noise in the data
% Moving average filter
windowSize = 1;
if ~isnan(windowSize)
    Dop_proc = movmean(Dop_blur, [windowSize 0], 3);
else
    Dop_proc = Dop_blur;
end

% Detrend the data
window_length = 300;
Dop_proc = w2b.preproc.detrendSlidingWindow(Dop_proc, window_length);

% Scale doppler to percent change
% Two ways to calculate percent change
method_to_use = 'baseline';
if strcmp(method_to_use, 'baseline')
    % Scale doppler to percent change from rest blocks
    Dop_scaled = w2b.preproc.scaleDoppler(Dop_proc, 'scaling_type', 'voxelwise_mean', ...
        'baseline_index', ~logical(task));

elseif strcmp(method_to_use, 'mean activity')
    % Scale doppler to percent change from mean voxel activity
    Dop_scaled = w2b.preproc.scaleDoppler(Dop_proc);
end


%% Create GLM regressors
% Convert into ms
fUS_startTime = min(timestamps(:, 1)) * 1000;
fUS_endTime = max(timestamps(:, 1)) * 1000;

% Task will be modeled as boxcar
boxcar_regressor_names = {'task on'};
boxcar_event_timing = {timestamps(task_start_ind)};
boxcar_event_duration = {timestamps(task_end_ind) - timestamps(task_start_ind)};

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
hrf = w2b.stats.generateHRF(n_secs, 'tau', 0.7, 'delta', 1, 'n', 3, 'Fs', 1000);


%% Create regressors
% Create regressors of interest
regressors = w2b.stats.createGLMRegressors(timestamps(:, 1), regressor_table, hrf.shape);

% Add constant regressor
regressors = [regressors ones(size(regressors, 1), 1)];

% Select regressors of interest
regressors_of_interest = regressors;
regressors_of_interest_ind = 1:2;


%% Build GLM
y = Dop_scaled;
patch_size = 1; % Size of patches in pixels/voxels
[~, ~, beta_hat] = w2b.stats.voxelwiseGLM(...
    y, regressors_of_interest, patch_size);


%% Apply specific T contrast
[pvalue_cont, tscore_cont] = deal(zeros(size(y, 1), size(y, 2), 8));
% For each regressor, apply the T contrast
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


%% Threshold contrasts for desired pvalue
p_threshold = 1e-5;

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

fig1 = figure;
tld = tiledlayout(1, 2);

nexttile();
% Generate the colormap of the t-scored activation to each direction
w2b.plot.plotDuplexImage(X_img_mm, Z_img_mm, tscore_thresholded(:, :, 1),angiogram,...
    'colormap2use', PositiveNegativeColormap, 'nonlinear_bg', 2, ...
    'AutoColorBarLimits', tValuesToDisplay, ...
    'showColorbar', true, ...
    'ColorBarTitle', 't-score', ...
    'data_aspect_ratio', data_aspect_ratio);
title_str = sprintf('%s', combined_event_labels{regressors_of_interest_ind(1)});
title(title_str, ...
    'Interpreter', 'none');


%% Plot mean activity from specified ROI across time
% How much of an ROI do you want?
window_size = 1.5; %In voxels

% Since the voxels are not isotropic, we need to define separate sizes for
% x and z
window_size_x = window_size;
window_size_z = window_size*.3/UF.Lambda;
window_size_mm_x = window_size_x * x_pixelsize;
window_size_mm_z = window_size_z * z_pixelsize;

% Reset T-scored map if ROI already drawn
figure(fig1);
nexttile(1);
if exist('roi_point', 'var')
    delete(roi_point);
end
roi_point = drawpoint;
if exist('rec', 'var')
    delete(rec)
end

% Draw new ROI onto T-scored map.
voxel_of_interest = roi_point.Position;
ROI_x_mm = voxel_of_interest(1)-window_size_mm_x:voxel_of_interest(1)+window_size_mm_x;
ROI_z_mm = voxel_of_interest(2)-window_size_mm_z:voxel_of_interest(2)+window_size_mm_z;
% Display ROI on plot
rec = rectangle('Position', [min(ROI_x_mm), min(ROI_z_mm), 2*window_size_mm_x, 2*window_size_mm_z], 'EdgeColor', 'w', 'LineWidth', 3);
drawnow;

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

% Scaled signal is centered around 100 by convention, so shift back to 0.
offset = -100;
y_offset = y + offset;

% Plot percent change across time
nexttile;
cla;
ylimits = [-15 25]; % In percent change
mean_activity = squeeze(mean(mean(y_offset(ROI_z, ROI_x, :), 2), 1));
task_scaled = task * (max(mean_activity) - min(mean_activity));
task_shifted_scaled = task_scaled + min(mean_activity);
plot(timestamps(:, 1) - fUS_startTime/1000, [mean_activity, task_shifted_scaled']);
xlabel('Time (s)');
ylabel('Voxelwise scaled fUS signal');
title('Timecourse of single fUS voxel with GLM fit overlay');
legend('mean Dop activity', 'task');
ylim(ylimits);
title(tld, sprintf('Session %s - Run %s', session_name, run_name), 'Interpreter', 'none');


%% How many voxels within ROI are visually modulated?
% Typically examine left and right LGN, but can pick any ROI.
fig2 = figure;

% Find range of T-score values.
max_abs_value = max(abs(tscore_thresholded(:, :, 1)), [], 'all');

% Generate the T-scored map
cmap = w2b.plot.plotDuplexImage(1:size(tscore_thresholded, 2), 1:size(tscore_thresholded, 1), tscore_thresholded(:, :, 1), angiogram,...
    'colormap2use', PositiveNegativeColormap, 'nonlinear_bg', 2, ...
    'showColorbar', true, ...
    'AutoColorBarLimits', tValuesToDisplay, ...
    'ColorBarTitle', 't-score', ...
    'data_aspect_ratio', data_aspect_ratio);
title('Select rectangular ROI to analyze how many activated voxels', ...
    'Interpreter', 'none');

% Draw rectangular
rect_roi = drawrectangle;
fprintf('Press enter to accept the drawn rectangular ROI.\n');
pause;

% Find how many voxels within desired ROI are statistically modulated
rect_roi = round(rect_roi.Position);
threshold_mask_cropped = threshold_mask(rect_roi(2):rect_roi(2)+rect_roi(4)-1, rect_roi(1):rect_roi(1)+rect_roi(3)-1, 1);
fprintf('Number of voxels significant within selected region is %d\n', nnz(~threshold_mask_cropped(:, :, 1)));


%% Compare activity during on/off blocks for the specified ROI
type_of_roi = 'rectangle';

% Convert Dop scaled to 1D so we can extract the drawn ROI
Dop_scaled_1D = reshape(Dop_scaled, [], size(Dop_scaled, 3));

% Initialize figure
figure;

% Create a button to stop the drawpoint feature
tb = uicontrol(gcf, 'Style', 'togglebutton',...
    'String', 'Stop updating graph');
% Plot angiogram
w2b.plot.plotAngiogram(angiogram, X_img_mm, Z_img_mm,...
    'colormap', 'inferno', ...
    'colorbar_title', 'Intensity (AU)', ...
    'data_aspect_ratio', data_aspect_ratio);
title('Draw a ROI around the desired region to compare');

% select point
roi_point = drawpoint;

[z_pixels, x_pixels] = size(angiogram);
while true
    switch type_of_roi
        case 'polygon'
            [~, x, y] = roipoly; %Select polygon anatomical polygon

            % Define conversion from world coordinates (mm) into image
            % coordinates(px)
            x_limits = xlim;
            x_spacing = (x_limits(2)-x_limits(1))/size(angiogram,2);
            y_limits = ylim;
            y_spacing = (y_limits(2)-y_limits(1))/size(angiogram,1);

            % Use coordinate conversion to convert vertices into image coordinates
            % y_limits/y_spacing used for initial offset (e.g. starting at 3
            % mm instead of 0 mm).
            x_px = x/x_spacing;
            y_px = y/y_spacing - y_limits(1)/y_spacing;

            %If polygon goes out of figure bounds, define edge to be figure bound.
            x(x>x_limits(2))=x_limits(2); x(x<x_limits(1))=x_limits(1);
            y(y>y_limits(2))=y_limits(2); y(y<y_limits(1))=y_limits(1);

            x_px(x_px>x_pixels)=x_pixels; x_px(x_px<0)=0;
            y_px(y_px>z_pixels)=z_pixels; y_px(y_px<0)=0;

            %Draw anatomical polygon onto subplot
            hold on;
            hline = plot(x, y,...
                ['g' '.-'],...
                'MarkerSize', 15);
            hold off;

            % Store the anatomical polygons
            AnatomicalPolygon_info = [x_px, y_px];

            % Smooth ROI
            smoothed_points = fnplt(cscvn([AnatomicalPolygon_info(:, 1)' AnatomicalPolygon_info(1, 1);
                AnatomicalPolygon_info(:, 2)'  AnatomicalPolygon_info(1, 2)]));

            AnatomicalPolygon_info = smoothed_points';

            % Create mask
            ROI_poly = polyshape(AnatomicalPolygon_info);
            ROI_mask = poly2mask(ROI_poly.Vertices(:, 1), ROI_poly.Vertices(:, 2), ...
                z_pixels, x_pixels);

        case 'rectangle'
            window_size = 1.5; %In voxels
            window_size_x = window_size;
            window_size_z = window_size*.3/UF.Lambda;
            window_size_mm_x = window_size_x * x_pixelsize;
            window_size_mm_z = window_size_z * z_pixelsize;

            if exist('rec2', 'var')
                delete(rec2);
            end

            voxel_of_interest = roi_point.Position;
            ROI_x_mm = voxel_of_interest(1)-window_size_mm_x:voxel_of_interest(1)+window_size_mm_x;
            ROI_z_mm = voxel_of_interest(2)-window_size_mm_z:voxel_of_interest(2)+window_size_mm_z;
            % Display ROI on plot
            rec2 = rectangle('Position', [min(ROI_x_mm), min(ROI_z_mm), 2*window_size_mm_x, 2*window_size_mm_z], 'EdgeColor', 'w', 'LineWidth', 3);
            drawnow;

            % Convert from mm to pixel space
            [~, Xidx] = min(abs(X_img_mm-voxel_of_interest(1)));
            [~, Zidx] = min(abs(Z_img_mm-voxel_of_interest(2)));

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
            if any(ROI_x > size(angiogram, 2))
                ROI_x(ROI_x > size(angiogram,2)) = size(angiogram, 2);
            end
            if any(ROI_z > size(angiogram, 1))
                ROI_z(ROI_z > size(angiogram, 1)) = size(angiogram, 1);
            end

            % Define ROI mask
            ROI_mask = false(size(angiogram));
            ROI_mask(ROI_z, ROI_x) = true;
    end
    if (get(tb,'Value')==1)    % If stop button is pressed, stop updating the graphs
        delete(roi_point);
        delete(tb);
        break;
    end
end

% Reshape to 1D
ROI_mask_1D = reshape(ROI_mask, [], 1);

% Pull out voxels in the mask only
ROI_activity = Dop_scaled_1D(ROI_mask_1D, :);

% Calculate roi statistics
fprintf('\n--------Statistical results --------\n');
task_ind = logical(task);
rest_ind = ~logical(task);
activity_task = mean(ROI_activity(:, task_ind), 2);
activity_rest = mean(ROI_activity(:, rest_ind), 2);
mean_difference = mean(activity_task) - mean(activity_rest);
sem_difference = std(activity_task - activity_rest)/sqrt(size(activity_task, 1));

% Display ROI statistics
fprintf('mean difference - %0.3f +/- %0.3f %%, \n', mean_difference, sem_difference);
