function [pvalue, tscore, beta] = voxelwiseGLM(...
    Dop, regressors, patch_size, varargin)

% Returns the beta value, p-value and t-score for each non-overlapping 
% voxel patch for the GLM using `regressors` and the fUS data
% 
% Inputs:
%   * Dop: 3D Doppler structure; (x, z, time)
%   * regressors: (n, p); n is number of observations (timepoints) and p is
%                  number of predictor variables
%   * patch_size: Size of the patch to be used to calculate the statistics
% Optional inputs:
%   * 'multiple_comparison_method' - 'FDR', 'bonferroni', or 'none';
%                                    Default is 'FDR'
%   * `verbose` - Display output to the console/command window.
% Outputs:
%   * pvalue - The p-value associated with each voxel and regressor under the GLM
%   * tscore - The t-score associated with each voxel and regressor under the GLM
%   * beta - The beta value(s) associated with each voxel and regressor under the GLM
% 
% Written by Whitney Griggs on August 31, 2021

% Suppress some warnings
warning('off','stats:glmfit:IterationLimit');

% Parse variable inputs
p = inputParser;
validDopMat = @(x) length(size(x)) == 3;
validPatchSize = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addRequired(p, 'Dop',  validDopMat);
addRequired(p, 'patch_size', validPatchSize);
addOptional(p, 'multiple_comparison_correction', 'FDR');
addOptional(p, 'verbose', true);
parse(p, Dop, patch_size, varargin{:});
result = p.Results;

% Reassign some of the variables for convenience
MC_correction = result.multiple_comparison_correction;

% Extract size dimensions from Dop1 and Dop2
[nDepth, nWidth, ~] = size(Dop);

% Allocate space for the output images
[pvalue, tscore, beta] = deal(NaN(nDepth, nWidth, size(regressors, 2)));

% Determine tiling of the patches
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

% Preallocate space for 1D pvalues
if strcmp(MC_correction, 'FDR')
    pvalue_2D = NaN(length(Zmat)*length(Xmat), size(regressors, 2));
end

counter = 0;

% For speed, likely want to parallelize
for Z = Zmat
    if result.verbose
        waitbar(Z / length(Zmat))
    end
    for X = Xmat
        iWidth = X - (patch_size - 1)/2 : X + (patch_size - 1)/2;
        iDepth = Z - (patch_size - 1)/2 : Z + (patch_size - 1)/2;
        
        patch_activity = squeeze(mean(mean(Dop(iDepth, iWidth, :), 1), 2));
        [B, ~, stats]=glmfit(regressors, patch_activity', 'normal', 'constant', 'off'); % uses matlab glmfit function
        beta(iDepth, iWidth, :) = B;
        tscore(iDepth, iWidth, :) = stats.t;
        pvalue(iDepth, iWidth, :) = stats.p;
        
                
        if strcmp(MC_correction, 'FDR')
            counter = counter + 1;
            pvalue_2D(counter, :) = stats.p;
        end
    end
end

switch MC_correction
    case 'FDR'
         % See above for meaning of following warning message: "The estimated
        % PI0 is greater than 1. Please check the p-values are valid or 
        % try a different lambda method. PI0 is set to 1."
        
        % Apply FDR correction. Since this is not necessarily done at the
        % voxel-level, we need to calculate this based upon the # of unique
        % tests ran, not the number of voxels.
        pvalue_1D = reshape(pvalue_2D, [], 1);
        [Q] = mafdr(pvalue_1D, 'BHFDR', true);
        pvalue_2D_corrected = reshape(Q, size(pvalue_2D));
        counter = 0;
        for Z = Zmat
            for X = Xmat
                iWidth = X - (patch_size - 1)/2 : X + (patch_size - 1)/2;
                iDepth = Z - (patch_size - 1)/2 : Z + (patch_size - 1)/2;
                
                counter = counter + 1;
                pvalue(iDepth, iWidth, :) = pvalue_2D_corrected(counter, :);
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