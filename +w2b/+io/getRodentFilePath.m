function [path_name, session_name, run_name] = getRodentFilePath(varargin)
% [path_name, session_name, run_name] = getRodentFilePath(varargin)
% Retrieve path to the requested rodent datafile
%
% INPUTS:
%   varargin
%       'project_record_filename'       string; Name of the rodent metadata file
%       'suppress_screen_output'        bool; Do you want to output a list of the
%                                       metadata to the command line?
%       'project_record_index'          double; If you know the session/run #, you
%                                       can skip the session/run selector.
% 
% Outputs:
%   path_name                           string; Absolute path to where fUS 
%                                       data is stored
%   session_name                        string; Name of the session; YYYYMMDD
%   run_name                            string; Name of the specific run

%% Input parser
p = inputParser;
addOptional(p, 'project_record_filename','rodent_session_record.json');
addOptional(p, 'suppress_screen_output', false);
addOptional(p, 'project_record_index', NaN);
parse(p, varargin{:});
inputs = p.Results;

%% load session/run metadata
project_record_filename = inputs.project_record_filename;
ProjectRecord = w2b.io.loadJSONAsTable(project_record_filename);

%% get the sessions & runs desired
% If none specified, ask the user to select.
if isnan(inputs.project_record_index)
    indx = w2b.io.rodentSessionSelector('project_record_filename', project_record_filename, ...
        'suppress_screen_output', inputs.suppress_screen_output);
else
    indx = inputs.project_record_index;
end

% Get the saved information about where the data is stored.
data_path = w2b.io.getUserDataPath;

% Build the full path name
path_name = fullfile(data_path, 'rodent', ProjectRecord.session_folder{indx}, ProjectRecord.run_folder{indx});

% Extract the session and run names
session_name = ProjectRecord.session_folder{indx};
run_name = ProjectRecord.run_folder{indx};