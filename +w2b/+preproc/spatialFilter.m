 function dataOut = spatialFilter(dataIn, filter_options)  
% Apply specified spatial filter
% First entry in `filter_options` cell defines the spatial filter
% Subsequent entries define the parameters for that specified
% filter.
%
% For 'disk' - specify filter size
% For 'gaussian - specify filter size and filter sigma

%% Define the spatial filter
[~, ~, n_windows, n_trials] = size(dataIn); 
dataOut = NaN(size(dataIn));
switch filter_options{1}
    case 'disk'
        h = fspecial(filter_options{1}, filter_options{2}); % disk size defined here
    case 'gaussian'
        h = fspecial(filter_options{1}, filter_options{2}, filter_options{3}); % disk size defined here
    otherwise
        error('This filter type has not been implemented');
end


% can't filter2 n-d arrays -> ugly nested for-loop (sorry)
% note: NOT faster with parfor overhead (I tried)
for window = 1:n_windows
    for trial = 1:n_trials
        dataOut(:, :, window, trial) = ...
            filter2(h, squeeze(dataIn(:, :, window, trial)));     
    end
end  
