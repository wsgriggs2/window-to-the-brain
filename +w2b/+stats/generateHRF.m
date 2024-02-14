function hrf = generateHRF(duration, varargin)
% Generate hemodynamic response function
% 
% Boynton HRF - From Geoffrey M. Boynton, Stephen A. Engel, Gary H. Glover and David J. Heeger
%   (1996). Linear Systems Analysis of Functional Magnetic Resonance Imaging in
%   Human V1. J Neurosci 16: 4207-4221
%
% Code based upon nipy.nitime.fmri.hrf implementation

%% Parse optional inputs
% Default values based upon Kagan et al. 2010 fMRI paper for NHPs
p = inputParser;
addOptional(p, 'tau', 1); % Time constant of the gamma function
addOptional(p, 'n', 3); %Phase delay of the gamma function. 
addOptional(p, 'delta', 1); % A pure delay. Allowing for an additional delay from onst of time-series to the beginning of gamma HRF
addOptional(p, 'Fs', 1); % Sampling rate
parse(p, varargin{:});
result = p.Results;

%% Calculate how many timepoints
t_max = duration - result.delta;

%% Generate the timepoints
t = [zeros(1, round(result.delta * result.Fs)), linspace(0, t_max, round(t_max * result.Fs))];

%% Compute t_tau, a variable used multiple times.
t_tau = t / result.tau;

%% Generate the HRF according to the boynton formula
% This is already normalized.
hrf.shape = (t_tau .^ (result.n - 1) .* exp(-1 .* (t_tau)) ./ (result.tau * factorial(result.n - 1)));

%% Store the parameters along with the hrf
hrf.fs = result.Fs;
hrf.tau = result.tau;
hrf.delta = result.delta;
hrf.n = result.n;