function [mdisten, varargout] = MDistEn(X,m,tau,varargin)

% Syntax:
%   [mdisten] = MDistEn(X,m,r,n,tau)
%
%   INPUTS:
%   X: matrix of the c-variate time series - a matrix of size c (number of 
%   channels) x N (number of sample for each channel)
%   m: vector of embedding dimensions - size c x 1, has to follow the same 
%   order of the channels
%   tau: vector of time delays - size c x 1, has to follow the same order
%   of the channels in X
%   B (optional): number of bins in the histogram to draw the distribution
%
%   OUTPUTS:
%   mdisten: values of multivariate distribution entropy 
%   (optional): number of points in the multivariate phase space 
%
% ________________________________________________________________________
%
% File:                         MDistEn.m
% Last revised:                 12 Apr 2024
% ________________________________________________________________________
%
% Copyright (C) 2024 Andrea Gargano, Mimma Nardelli
%
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% You can contact the author by e-mail: 
% andrea.gargano@phd.unipi.it, andrea.grgn96@gmail.com
% ________________________________________________________________________
%
% This method was first proposed in:
% A. Gargano, M. Nardelli, E.P. Scilingo, "Exploring Multivariate Dynamics 
% of Emotions Through Time-Varying Self-Assessed Arousal and Valence Rating"
% 
%
% If you use this program in support of published research, please include
% a citation of the reference above. If you use this code in a software 
% package, please explicitly inform the end users of this copyright notice 
% and ask them to cite the reference above in their published research.
% ________________________________________________________________________
%
% ________________________________________________________________________

%% sanity checks

% check number of inputs
narginchk(3,4)

% check assignment of optional inputs
if ~isempty(varargin)
    [B] = deal(varargin{:});
end

% check number of outputs
nargoutchk(1,2)

% check size of inputs
[C, ~] = size(X);
[c_m, r_m] = size(m);
[c_t, r_t] = size(tau);
% check consistency of m and tau
if ~(r_m == r_t)
    error('Incorrect call to MDistEn: check dimensions of inputs!')
end
% check consistency of m, tau, and X
if ~ ((c_m == c_t) && (c_m == C))
    error('Incorrect call to MDistEn: check dimensions of inputs!')
end

%% Initialization

% choose maximum embedding dimension
m_max = max(m);
% ensure having at least 2 dimensions 
if m_max == 1
    m_max = m_max+1;
    warning('Maximum embedding dimension is now 2!')
end

%% Multivariate phase space reconstruction 

% Multichannel signals may have different amplitude ranges: 
% normalization via zscore
X = zscore(X,0,2);

% reconstruct the multivariate phase space
[X_mv, N_de] = multivariate_phase_space(X, m_max, tau);
% number of reconstructed points 
varargout{1} = N_de;

%% DistEn of multivariate phase space

% compute Chebychev distance 
dv = pdist(X_mv, 'chebychev');

% check external setting on number of bins (B)
if ~exist('B','var')
    % compute number of bins via Freedman-Diaconis Rule 
    
    % number of point-distances
    N = size(dv,2);
    
    % robust measure of spread via iqr
    width = iqr(dv);
    % estimate bin width
    width = 2*width*N^(-1/3);

    % compute number of bins in base 2
    B = ceil(range(dv)/width);
    B = 2^ceil(log2(B));
end

% esimate probability density by histogram
num = hist(dv, linspace(0, 1, B));
freq = num./sum(num);
% % use the following lines if the function hist is no more supported 
% % (on average both methods give similar results, single values may be 
% % subject to variations)
% edges = linspace(0, 1, B);
% halfedge = mean(diff(edges));
% edges = [-halfedge, halfedge+edges];
% freq = histcounts(dv, edges,"Normalization","probability");

% mdisten calculation
mdisten = -sum(freq.*log2(freq+eps)) ./ log2(B);

end
