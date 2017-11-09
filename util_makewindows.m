% pe_makewindows() - Epoch data.
% 
% Inputs: 
%   data      (vector) - time series to epoch
%   windowlen (int)    - length of each epoch in samples
%   noverlap  (int)    - number of samples overlapping between windows
%   srate     (int)    - sampling rate (Hz)
%
% Output:
%   datawin   (matrix) - epoched data
%   time      (vector) - time at center of epochs
% 
% Usage: [datawin, time] = pe_makewindows(data, windowlen, noverlap, srate);

function [datawin, time] = util_makewindows(data, windowlen, noverlap, srate)
if nargin < 4
    error('Need data, window size, sample overlap, and sampling rate.')
end

if ~isvector(data)
    error('Only one channel of data expected.')
end

% make data a column vector
data = data(:);

% total number of data points
nsample = length(data); 

% number of windows to make
nwindow = fix((nsample - noverlap)/(windowlen - noverlap));

% sample inds at center of windows
windowind = 1 + (0:(nwindow-1))*(windowlen-noverlap);

% time points at center of windows
time = ((windowind-1)+((windowlen)/2)')/srate; 
time = time(:);

% window the data
datawin = nan(windowlen, nwindow);
for w = 0:nwindow-1
    ind0 = (windowlen - noverlap)*w+1;
    datawin(:,w+1) = data(ind0:ind0+windowlen-1);
end
end