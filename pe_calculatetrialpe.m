% pe_calculatetrialpe() - Calculate permutation entropy (PE) across trials
%                         Calls pe_calculatechannelpe.m
%
% Inputs: 
%   data        (npnt x ntrial) - epochs for which to calculate PE
%   nsample     (int)           - number of samples to be converted to a symbol
%   noverlap    (int)           - number of samples overlapping between symbols; default is nsample - 1
%   srate       (int)           - sampling rate (Hz)
%   weighted    (boolean)       - if true, calculated PE weighted by amplitude; else, non-weighted PE
%
% Outputs:
%   peseries    (vector) - permutation entropy time series for data
%   time        (vector) - center of time periods used to calculate peseries
%
% Usage: 
%   [peseries, time] = pe_calculatetrialpe(data, nsamp, noverlap, srate, weighted);
%
function [peseries, time, symbolseries] = pe_calculatetrialpe(data, nsamp, noverlap, srate, weighted)
if nargin < 5
    error('All inputs are required.')
end

if isvector(data)
    error('Epoched data expected. Use pe_calculatechannelpe.m for single time series.')
end

% default noverlap
if isempty(noverlap)
    noverlap = nsample - 1;
end

% get number of samples and trials
[npnt, ntrial] = size(data);

% get indices of samples to use per symbol
[sampind, time] = util_makewindows(1:npnt, nsamp, noverlap, srate);
npe = length(time);

% calculate PE across trials per window
peseries = nan(npe, 1);
symbolseries = nan(npe, ntrial);

for ipe = 1:npe
    sampcurr = data(sampind(:, ipe), :);
    
    % get window across trials, unwrap to column, and pass to pe_calculatechannelpe 
    % use same window length and number of trials as the number of windows to calculate PE over
    [peseries(ipe), ~, symbolseries(ipe, :)] = pe_calculatechannelpe(sampcurr(:), nsamp, 0, ntrial, srate, weighted);
end
end