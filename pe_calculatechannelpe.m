% pe_calculatechannelpe() - Calculate permutation entropy (PE) of a single time series
%
% Inputs: 
%   data        (vector)  - time series for which to calculate PE
%   nsample     (int)     - number of samples to convert into a symbol
%   noverlap 	(int)     - number of samples overlapping between symbols; default is nsample - 1
%   nsymbol     (int)     - number of symbols over which to calculate permutation entropy
%   srate       (int)     - sampling rate (Hz)
%   weighted    (boolean) - if true, calculated PE weighted by amplitude; else, non-weighted PE
%
% Outputs:
%   peseries    (vector) - (weighted) PE time series for data
%   time        (vector) - center of time periods used to calculate peseries
%   symbolwin   (matrix) - windowed symbols used to calculate (weighted) PE
%
% Usage: 
%   [peseries, time] = pe_calculatechannelpe(data, nsample, noverlap, nsymbol, srate, weighted);
%
function [peseries, time, symbolwin, varwin] = pe_calculatechannelpe(data, nsample, noverlap, nsymbol, srate, weighted)
if nargin < 6
    error('All inputs expected.')
end

if ~isvector(data)
    error('One vector of data is expected.')
end

% default noverlap
if isempty(noverlap)
    noverlap = nsample - 1;
end

% epoch the data into windows of length nsample
% each column is a window
[datawin, time] = util_makewindows(data(:), nsample, noverlap, srate);
nwin = size(datawin, 2); % total number of windows/symbols to be constructed in data

% calculate variance in each window
varwin = var(datawin, 1);

% rank the data in each window
rankwin = nan(nsample, nwin);
for iwin = 1:nwin
    [~, ~, rankwin(:, iwin)] = unique(datawin(:, iwin));
end

% convert the ranks to symbols (literally, combine ranks into one number)
symbols = zeros(1, nwin);
for isamp = 1:nsample
    symbols = symbols + (rankwin(isamp,:) .* (10 ^ ((nsample - isamp + 1) - 1)));
end

% epoch symbols into windows of length nsymbol
symbolwin = util_makewindows(symbols(:), nsymbol, 0, srate);
varwin = util_makewindows(varwin(:), nsymbol, 0, srate);
time = mean(util_makewindows(time(:), nsymbol, 0, srate));
npe = size(symbolwin, 2); % number of PE samples to be calculated

% calculate permutation entropy per time period
peseries = zeros(1, npe);

% for each time period... 
for ipe = 1:npe
    % loop over unique symbols in each time period
    uniquesymbol = unique(symbolwin(:, ipe));
    
    for isym = 1:length(uniquesymbol)
        
        % get indices of particular symbol
        symbolscurr = symbolwin(:, ipe) == uniquesymbol(isym);
        
        if weighted
            % calculate partial variance explained by particular symbol in t
            psymbol = sum(varwin(symbolscurr, ipe))/sum(varwin(:, ipe));
        else
            % calculate probability of symbol occurence in t
            psymbol = sum(symbolscurr)/length(symbolscurr);
        end
        
        % calculate entropy for t
        if psymbol ~= 0 % log2(0) = -inf
            peseries(ipe) = peseries(ipe) - (psymbol * log2(psymbol));    
        end
    end
end
end