function [Dop_scaled, scaling_transform] = scaleDoppler(Dop, varargin)
% [Dop_scaled, scaling_transform] = scaleDoppler(Dop, varargin)
% Function to scale Doppler data.
% Two possible methods
% 1. 'Voxelwise_mean' - Used by AFNI
%     This function uses the same method as AFNI for scaling of data to a
%     common scale. See Chen, Taylor, and Cox 2017 under Units section.
% 
% 2. Scales to mean of the voxel.
%    'Grand_mean_scaling' - Used by SPM and FSL
%     Normalize entire image to have mean of 1. 
%
% INPUTS
%   Dop                     3D doppler (Z x X x time)
%   varargin
%       scaling_type        'voxelwise_mean' or 'grand_mean_scaling'
%       baseline_index      Specify a custom baseline, otherwise
%                           this defaults to mean.
%
% OUTPUTS
%   Dop_scaled              Scaled 3D doppler
%   scaling_transform       Transform used to scale the data

%% Input parser
p = inputParser;
addOptional(p, 'scaling_type', 'voxelwise_mean');
addOptional(p, 'baseline_index', []);
parse(p, varargin{:});
result = p.Results;


%% Z by X by timepoints
[n_depth, n_width, n_timepoints] = size(Dop);


%% flatten the data into a voxel-column
dataOut = permute(Dop, [2, 1, 3]);% [n_width, n_depth, n_timepoints);
dataOut = reshape(dataOut, [n_depth*n_width, n_timepoints]);


%% perform the scaling
switch result.scaling_type
    case 'voxelwise_mean'
        if isempty(result.baseline_index)
            voxel_means = mean(dataOut, 2);
        else
            voxel_means = mean(dataOut(:, result.baseline_index), 2);
        end
        voxel_means_mat = repmat(voxel_means, 1, size(dataOut, 2));
        dataOut = 100*dataOut./voxel_means_mat;
        scaling_transform = voxel_means;
    case 'grand_mean_scaling'
        grand_mean = mean(dataOut, 'all', 'omitnan');
        dataOut = dataOut./grand_mean;
        scaling_transform = grand_mean;
    otherwise 
        error('Not an accepted method for scaling');
end

%% reshape back to original data dims
dataOut = reshape(dataOut, [n_width, n_depth, n_timepoints]);
Dop_scaled = permute(dataOut, [2, 1, 3]);
    
