function project_record_index = rodentSessionSelector(varargin)
% project_record_index = rodentSessionSelector(varargin)
% Select desired session/run combo based upon metadata file
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
%   base_path:      string; Absolute path to where fUS data is stored

%% Input parser
p = inputParser;
addOptional(p, 'project_record_filename','rodent_session_record.json');
addOptional(p, 'suppress_screen_output', false);
addOptional(p, 'selection_mode', 'single')
parse(p, varargin{:});
inputs = p.Results;

%% load session/run metadata
project_record_filename = inputs.project_record_filename;
ProjectRecord = w2b.io.loadJSONAsTable(project_record_filename);
if ~inputs.suppress_screen_output
    disp(ProjectRecord)
end


%% create the list of available sessions/runs
runStrings = cell(size(ProjectRecord,1),1);
for i = 1:size(ProjectRecord,1)
    run_folder = ProjectRecord.run_folder{i};
    Condition = ProjectRecord.condition(i);
    session_folder = ProjectRecord.session_folder{i};
    runStrings{i} = ['Session ' session_folder, ',Run ' run_folder ', Condition ' num2str(Condition)];
end


%% use a GUI to select which ones you want
[project_record_index,~] = listdlg('PromptString','Select Session/Run(s) to load:',...
    'SelectionMode',inputs.selection_mode,...
    'ListString',runStrings, ...
    'ListSize', [600 300]);