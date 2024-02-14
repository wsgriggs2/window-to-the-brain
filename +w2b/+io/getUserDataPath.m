function data_path = getUserDataPath()
% path = getUserDataPath()
% Retrieve path to the user's data
%
% INPUTS:
%   none
% 
% Outputs:
%   path:               string; Absolute path to where fUS data is stored


%% Check if root has already been set
if isfile('data_path.mat')
    % If path has already been specified, then load that.
    load('data_path.mat', 'data_path');
else
    % If path has not been set, then ask user to specify a location
    data_path = w2b.io.specifyDataPath;
end