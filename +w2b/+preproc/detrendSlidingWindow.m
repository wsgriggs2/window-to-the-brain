function dop  = detrendSlidingWindow(dop, window_length)
% dop  = detrendSlidingWindow(dop, window_length)
%
% detrendDrift taks in dop, which is a 4-D array of dimensions
% yPixels x xPixels x nWindows x nTrials
%
% dop can also be 3-D, i.e. yPixels x xPixels x time
%
% window_length is the window length used in the sliding window, def: 30
%
% This function computes the moving mean (see movmean) for each pixel in
% the image in nested loops. 

[yPixels, xPixels, nWindows, nTrials] = size(dop);

%% handle varargin
if ndims(dop)==4
    % reshape iDop (note: I confirmed that this reshape operation is correct by
    % building nested for loops of trials/windows 
    dop = reshape(dop,yPixels,xPixels,nTrials*nWindows);
    reshape_flag = true;
else
    reshape_flag = false;
end
    
if ndims(dop)~=3
    warning('doppler data is badly sized!')
    return
end

if ~exist('window_length','var')
    window_length = 30;
end

%% the business
gcp;
fprintf('Detrending now (may take a second)...');
for y = 1:yPixels    
    parfor x = 1:xPixels
        this_vector = squeeze(dop(y,x,:));
        % note that I use a directional window length [kb kf] with kf=0
        % because, in a real time scenario, we won't be able to see into
        % the future. <<< movmean usage below >>>
        window_avg = movmean(this_vector,[window_length 0]);
        avg_avg = mean(squeeze(dop(y,x,:))); 
        dop(y,x,:) = permute(dop(y,x,:),[3 2 1]) - window_avg;
        % here I add the avg activation of this pixel (across the whole 
        % session) back in so we don't have negative values everywhere & we
        % preserve the raw doppler magnitude (means).
        dop(y,x,:) = dop(y,x,:) + avg_avg;
    end    
end

%% reshaping when necessary
if reshape_flag
    dop = reshape(dop,[yPixels xPixels nWindows nTrials]);
end
end
        