function regressors = createGLMRegressors(fUS_timestamps, events, hrf, varargin)
%% Generate GLM regressors for fUS data
% This function creates GLM regressors. Will generate regressor for each
% provided column in event_timestamp Inputs:
%    *  fUS_timestamps - vector (timepoints x 1); The acquisition time of
%                        each fUS image in seconds.
%   *   events - Table array where each row contains at least three
%                variables (timestamps, duration, 'modulation'). `timestamps`
%                indicates when the event begins. 'duration' indicates the
%                length of each event. Duration should have a duration for
%                each event, i.e. same size as timestamps. If duration is
%                empty for a row, then assumed that we should use a stick
%                function to represent event times. Modulation implicitly
%                used to demarcate main vs parametric modulators. If empty,
%                then main regressor. If has values, should be same size as
%                timestamps.
%   *   hrf - Matrix (n x T) representing HRF at millisecond time
%             resolution. n is number of basis vectors
% Outputs:
%   *   regressors - cell array where each entry is a variable length
%                    vector representing the HRF-convolved signal at the
%                    fUS timestamps.
%
% Written by Whitney Griggs. Sept. 1, 2021

%% Main function body

% Parse variable input
p = inputParser;
addOptional(p, 'verbose', false);
parse(p, varargin{:});
result = p.Results;

% Calculate the start and end time of the fUS acquisition Round should not
% be necessary, but including to ensure that no decimals included. Decimals
% introduce round-off error that causes problems with constructing the
% index of when events/fUS occurs later in this script.
fUS_startTime = round(min(fUS_timestamps) * 1000);
fUS_endTime = round(max(fUS_timestamps) * 1000);

% Create a high-resolution timelog between start and end time.
highRes_time = fUS_startTime:fUS_endTime;

% How many regressors do we have?
nRegressors = height(events);

% Find indices of the fUS slices
fUS_timeInd = ismember(highRes_time, round(fUS_timestamps*1000));

% Create indices to represent events
event_ind = zeros(length(highRes_time), nRegressors);

for regressor = 1:nRegressors  
    % Extract the event onset timestamps for the specified events (in
    % second initially)
    event_times = round(cell2mat(events{regressor, 'timestamps'}) * 1000);
    
    % Extract the event duration for the specified events (in seconds
    % initially)
    event_durations = round(cell2mat(events{regressor, 'duration'}) * 1000);
    
    % Extract modulation (if any)
    event_modulations = cell2mat(events{regressor, 'modulation'});
    
    % Check whether the modulation variable is already mean-centered and
    % scaled between [-1 1]. If not, then perform mean-centering and
    % scaling.
    if ~isempty(event_modulations)
        if mean(event_modulations) ~= 0  || max(event_modulations) ~= 1
            % Mean center the parameter
            event_modulations = event_modulations - mean(event_modulations);
            
            % Rescale so max = 1
            event_modulations = event_modulations/max(event_modulations);
        end
    end
    
    % Find indices of event onset in the high-resolution time record
    timeInd = ismember(highRes_time, event_times);
    events_to_keep = ismember(event_times, highRes_time);

    
    if isempty(event_durations)
        % If no duration, then it is modeled as a stick function
        if isempty(event_modulations)
            % If no modulation, then it is a main effect regressor.
            event_ind(timeInd, regressor) = 1;
        else
            event_modulations = event_modulations(events_to_keep);
            event_ind(timeInd, regressor) = event_modulations;
        end
    else
        
        % If there is a duration, then the duration can be different for
        % each event and needs to be handled separately.
        timeInd_compact = find(timeInd);
        for event = 1:length(timeInd_compact)
            % Extract duration of each specific event
            duration = round(event_durations(event));
            
            % Create new index representing the total duration of the
            % event.

            timeInd_withDuration = timeInd_compact(event):(timeInd_compact(event) + duration);

            if isempty(event_modulations)
                % If no modulation, then it is a main effect regressor.
                event_ind(timeInd_withDuration, regressor) = 1;
            else
               event_modulations = event_modulations(events_to_keep);
               event_ind(timeInd_withDuration, regressor) = event_modulations(event); 
            end
        end
    end
end % End for loop over regressors

% Convolve HRF with the regressor.
% This takes some time to run.
fprintf('Convolving %d regressors with the HRF. This may take some time. \n', nRegressors);

% Pre-allocate space
events_convolved = NaN(size(event_ind, 1), size(hrf, 1)*nRegressors);

% Iterate over each hrf basis vector
for basis_num = 1:size(hrf, 1)
    fprintf('Convolving HRF basis #%d with regressors. \n', basis_num);
    events_convolved(:, basis_num:size(hrf, 1):end) = filter(hrf(basis_num, :), 1, event_ind);
end
    
% Pull out the timepoints of interest
regressors = events_convolved(fUS_timeInd, :);

% Debug plot
if result.verbose
    figure;
    % This will create the stereotypical design matrix.
    imagesc(regressors); 
    colormap('gray');
    ylabel('fUS acquisition number');
    xlabel('Regressors');
    xticks(1+(size(hrf, 1)-1)/2:size(hrf, 1):(size(hrf, 1)*nRegressors));
    xticklabels(events{:, 'regressorName'});
    title('Design matrix');
    
    figure;
    % This will create correlation coeffiients for the regressor to see how
    % correlated they are. 
    corr_matrix = corrcoef(regressors);
    imagesc(corr_matrix);
    colorbar;
    ylabel('Regressor');
    xlabel('Regressor');
    xticks(1+(size(hrf, 1)-1)/2:size(hrf, 1):(size(hrf, 1)*nRegressors));
    yticks(1+(size(hrf, 1)-1)/2:size(hrf, 1):(size(hrf, 1)*nRegressors));
    xticklabels(events{:, 'regressorName'});
    yticklabels(events{:, 'regressorName'});
    title('Correlation between regressors');
end
