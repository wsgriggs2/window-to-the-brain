function plotAngiogram(backgroundImage,X_img_mm, Z_img_mm, varargin)
% plotAngiogram(backgroundImage,X_img_mm, Z_img_mm, varargin)
% Plot the vascular anatomy, i.e., angiogram, from the current field of
% view
%
% INPUTS
%   backgroundImage         (zPix x xPix) array/image
%   X_img_mm                x axis indices in mm
%   Z_img_mm                z axis indices in mm
%   varargin
%       title               title for the plot
%       colormap            colormap to use for the visualization
%       show_colormap       bool; do you want a colorbar on the side?
%       colorbar_title      string; title for the colorbar (if used)
%       data_aspect_ratio   (1x3) vector to define voxel size relationship
%
% OUTPUTS
%   none

 
%% Input parser
p = inputParser;
p.addOptional('title', 'Average Doppler Image');
p.addOptional('colormap', 'gray');
p.addOptional('show_colormap', true);
p.addOptional('colorbar_title', '');
p.addOptional('data_aspect_ratio', [1 1 1]);
p.parse(varargin{:});
inputs = p.Results;


%% Display wallpaper
imagesc(X_img_mm, Z_img_mm, backgroundImage);
daspect(inputs.data_aspect_ratio);
axis image;
xlabel('mm');
ylabel('mm');
colormap(inputs.colormap);
title(inputs.title, 'interpreter', 'none');

if inputs.show_colormap
    cb = colorbar;
    if ~isempty(inputs.colorbar_title)
       cb.Label.String = inputs.colorbar_title;
    end
end