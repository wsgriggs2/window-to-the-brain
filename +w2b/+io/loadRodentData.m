function data = loadRodentData(varargin)
% data = loadRodentData(varargin)
%
% INPUTS:
%   See w2b.io.getRodentFilePath for list of possible varargin
%
% Outputs:
%   data:      struct; Data associated with a single session/run combo

%% Get the full pathname for a session/run's data
[path_name, session, run] = w2b.io.getRodentFilePath(varargin{:});
filename = 'Dop.mat';
full_filename = fullfile(path_name, filename);


%% Load the requested data
load(full_filename, 'Dop');
load(fullfile(path_name, 'UF.mat'), 'UF');


%% Infer timestamps.
% Collected data at 1 Hz
timestamps = [1:size(Dop, 3)]';


%% Create angiogram to show vascular anatomy and field of view
angiogram = w2b.util.makeAngiogram(Dop, 5);


%% Package as struct
data.dop = Dop;
data.angiogram = angiogram;
data.UF = UF;
data.timestamps = timestamps;
data.session = session;
data.run = run;
