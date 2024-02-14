function data_path = specifyDataPath()
% specify_data_path Use GUI to specify path to the user's data and save
%                   path
%
% INPUTS:
%   none
%
% Outputs:
%   base_path:      string; Absolute path to where fUS data is stored


title_string = 'Specify path to the downloaded data folder:';
disp(title_string);
data_path = uigetdir(matlabroot, title_string);
save('data_path.mat', 'data_path');