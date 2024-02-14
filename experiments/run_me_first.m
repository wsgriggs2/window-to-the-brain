%% Setup script
% This script will let the user set the path to the downloaded data. The
% data can be found on CaltechDATA at https://doi.org/10.22002/f3y3k-em558

% Specify the data path
path = w2b.io.specifyDataPath;


% If no path set, let the user know
if ~ischar(path)
    error('You did not set a path to the data. Please re-run this script.');
end

% Tell the user what path they set
fprintf("You set '%s' as the path to the data. If this incorrect, rerun this script.\n", ...
    path);