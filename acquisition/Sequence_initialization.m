%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Functional Ultrasound Imaging of Freely Moving Adult Human Neural activity through an Acoustically Transparent Cranial Prosthetic                  %%
%   Claire Rabut, Sumner L. Norman, Whitney S. Griggs, Jonathan J. Russin, Kay Jann, Vasileios Christopoulos, Charles Liu, Richard A. Andersen, Mikhail G. Shapiro   %%

%% Parameter definition for functional Ultrasound Imaging Acquisition Script %

% These parameters will be pass into the Verasonics_sequence.p script to
% collect fUSI data %

% The Verasonics_sequence.p script needs to be executed on a high-frequency
% Verasonics 128 or 256 channels scanner
% UTA-260-S or UTA-260-D
% Vantage version: Vantage-4-4-0


% Define path to save data
path_to_vantage_folder = '\Users\Downloads\Vantage-4.4.0\' ; % Update with your path
path_save = 'D:\fUSI_sequence\'; % Update with your path
SaveName = 'Test';


% Acquisition Parameters

Number_Doppler_images = 300;            % Number of Doppler images to acquire for brain activity recording
ini_depth = 10;                         % Initial depth [mm]
final_depth = 50;                       % Final depth [mm]
Transmission_Frequency = 7.5;           % Transmit frequency in MHz
DutyCycle = 0.67;                       % Duty cycle 
NbHalfCycle = 2;                        % Number of half cycles in emission
Time_loop = 1;                          % Time loop in seconds
Angles_plane_waves = (-6:3:6)*pi/180;   % Angles of plane waves
Number_bmodes_per_Doppler = 300;        % Number of B-modes that constitute a Doppler image
BmodesFrameRate = 500;                  % Framerate of B-modes


% Script to send the variables to:
% "Verasonics_sequence.p" needs to be put in the root of the Vantage 4.4.0
% folder, i.e., "pathToVantageFolder\Verasonics_sequence.p".
Verasonics_sequence = fullfile(path_to_vantage_folder, 'Verasonics_sequence.p');
run(Verasonics_sequence);
