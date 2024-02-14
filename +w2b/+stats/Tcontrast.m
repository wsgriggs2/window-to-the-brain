function [pvalue, tscore] = Tcontrast(...
    Dop, regressors, beta_hat, t_contrast, patch_size, varargin)

% Calculate p-values associated with specific GLM T-contrasts


%% Parse variable inputs
p = inputParser;
validDopMat = @(x) length(size(x)) == 3;
validPatchSize = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p, 'Dop',  validDopMat);
addRequired(p, 'patch_size', validPatchSize);
addRequired(p, 'beta_hat');
addRequired(p, 't_contrast');
addOptional(p, 'multiple_comparison_correction', 'FDR');
addOptional(p, 'verbose', true);
parse(p, Dop, patch_size, beta_hat, t_contrast, varargin{:});
result = p.Results;

%% Reassign some of the variables for convenience
MC_correction = result.multiple_comparison_correction;

%% Extract size dimensions from Dop
[nDepth, nWidth, ~] = size(Dop);

%% Allocate space for the output images
[pvalue, tscore] = deal(NaN(nDepth, nWidth));

%% Determine tiling of the patches
if mod(nDepth, patch_size) == 0 && mod(nWidth, patch_size) == 0
    Zmat = (patch_size + 1) / 2 : patch_size : nDepth;
    Xmat = (patch_size + 1) / 2 : patch_size : nWidth;
elseif mod(nDepth, patch_size) == 0
    Zmat = (patch_size + 1) / 2 : patch_size : nDepth;
    Xmat = get_optimal_tiling(patch_size, nWidth);
elseif mod(nWidth, patch_size) == 0
    Zmat = get_optimal_tiling(patch_size, nDepth);
    Xmat = (patch_size + 1) / 2 : patch_size : nWidth;
else
    Zmat = get_optimal_tiling(patch_size, nDepth);
    Xmat = get_optimal_tiling(patch_size, nWidth);
end
if result.verbose
    h = waitbar(0,'Please wait ... computing stat mask');
end

%% Preallocate space for 1D pvalues
if strcmp(MC_correction, 'FDR')
    pvalue_1D = NaN(length(Zmat)*length(Xmat), 1);
end

counter = 0;
%% For speed, likely want to parallelize
for Z = Zmat
    if result.verbose
        waitbar(Z / length(Zmat))
    end
    for X = Xmat
        iWidth = X - (patch_size - 1)/2 : X + (patch_size - 1)/2;
        iDepth = Z - (patch_size - 1)/2 : Z + (patch_size - 1)/2;
        
        patch_activity = squeeze(mean(mean(Dop(iDepth, iWidth, :), 1), 2));
        
        patch_beta_hat = squeeze(beta_hat(iDepth, iWidth, :));
        
        error = (patch_activity - regressors*patch_beta_hat);
        
        var_e = error'*error/(length(patch_activity) - 1 - size(regressors, 2));
        
        t_stat = t_contrast*patch_beta_hat/sqrt(var_e*t_contrast*pinv(regressors'*regressors)*t_contrast');
        
        p_stat = 1-tcdf(t_stat, length(patch_activity) - size(regressors, 2));
        
        tscore(iDepth, iWidth) = t_stat;
        pvalue(iDepth, iWidth) = p_stat;
        
        if strcmp(MC_correction, 'FDR')
            counter = counter + 1;
            pvalue_1D(counter) = p_stat;
        end
    end
end

switch MC_correction
    case 'FDR'
        % Apply FDR correction. Since this is not necessarily done at the
        % voxel-level, we need to calculate this based upon the # of unique
        % tests ran, not the number of voxels.
        [Q] = mafdr(pvalue_1D, 'BHFDR', true);
        pvalues_corrected = Q;
        counter = 0;
        for Z = Zmat
            for X = Xmat
                iWidth = X - (patch_size - 1)/2 : X + (patch_size - 1)/2;
                iDepth = Z - (patch_size - 1)/2 : Z + (patch_size - 1)/2;
                
                counter = counter + 1;
                pvalue(iDepth, iWidth) = pvalues_corrected(counter);
            end
        end
    case 'Bonferroni'
        %Multiplying by number of tests done to create Bonferroni corrected
        %p-value.
        pvalues_corrected = pvalue * (length(Zmat)*length(Xmat)); 
        pvalues_corrected(pvalues_corrected > 1) = 1; % Probability cannot exceed 1.
        pvalue = pvalues_corrected;
    otherwise
        % Currently unnecessary, since no MC correction is applied
end
if result.verbose
    close(h)
    fprintf('Computed new statistics. patch_size=%d \n', patch_size);
end

% Turn warnings back on
warning('on','stats:glmfit:IterationLimit');
end

function patch_centers = get_optimal_tiling(patch_size, image_dimension)
% Generate the optimal tiling in 1D. If mod(image_dimension, patch_size) is
% even, then leave equal space on both sides. If this remainder is odd,
% then leave an extra pixel spacing)relative to the bottom or right spacing
% at the top or left of the image.
% 
% Inputs:
%   patch_size - The scalar representing patch size. Must be an int.
%   image_dimension - The 1D image size. Either depth or width.
% Outputs:
%   patch_centers - The center of the patch.

% Get remainder size
unusable_space = mod(image_dimension, patch_size);
spannable_distance = image_dimension - unusable_space;

% Determine optimal spacing
if mod(unusable_space,2)
    % if remainder is odd
   patch_centers = 1 + (unusable_space - 1) / 2 + (patch_size + 1) / 2 : patch_size : spannable_distance + 1;
else
    % if remainder is even
   patch_centers = unusable_space / 2 + (patch_size + 1) / 2 : patch_size : spannable_distance;
   
end
end