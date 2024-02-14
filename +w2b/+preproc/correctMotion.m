function [imageOut,template1] = correctMotion(imageIn, template)
% [imageOut,template1] = correctMotion(imageIn, template)
% Function that uses NormCorre to correct for motion across session
%
% INPUTS
%   imageIn         (Z x X x Time) vector; 3D doppler        
%   template        if you have a custom template you want to align all
%                   Doppler frames against
%
% OUTPUTS
%   imageOut        (Z x X x Time) vector; 3D doppler aligned
%   template1       template used to align against


%% Use NormCorre for motion correction

% Initialize parallel processing toolbox
gcp;

% If template defined, then use a template
if ~exist('template','var')
    useTemplate = false;
else
    useTemplate = true;
end

% Get size of the 3D doppler data
[yPixels, xPixels, nWindows] = size(imageIn);
Y = imageIn;

% Convert to single precision
Y = single(Y);
Y = Y - min(Y(:));


%% set parameters
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200, 'shifts_method', 'cubic');
options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200, 'shifts_method', 'cubic');


%% perform motion correction
if useTemplate
    tic; [M1,shifts1,template1,options_rigid] = normcorre(Y,options_rigid,template); toc %#ok<*ASGLU>
else
    tic; [M1,shifts1,template1,options_nonrigid] = normcorre_batch(Y,options_nonrigid); toc
end


%% pass back the result
disp('returning rigid-body motion corrected data')
imageOut = M1;
