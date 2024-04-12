function [X_med, l_med] = multivariate_phase_space(X,m,tau)
% Syntax:
% [Xsc_med, l_med] = multivariate_phase_space(X,m,tau)
%
% Reconstruct the (multivariate) median trajectory of the multichannel 
% input data X given the maximum embedding dimension (m) and the array of 
% channel-specific time delays (tau). This results into an m-dimensional 
% median trajectory encompassing all the number of channels in X.
%
% INPUTS: 
% - X, normalized c-variate time series - a matrix of size c (number of 
%   channels) x N (number of sample for each channel)
% - m, maximum embedding dimension - size 1
% - tau, vector of time delays - size c x 1, has to follow the same order
%   of the channels in X
%
% OUTPUTS: 
% - X_med, median reconstucted trajectory 
% - l_med, length of the median trajectory
%
% ________________________________________________________________________
%
% File:                         multivariate_phase_space.m
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
% If you use this program in support of published research, please include
% a citation of the reference above. If you use this code in a software 
% package, please explicitly inform the end users of this copyright notice 
% and ask them to cite the reference above in their published research.
% ________________________________________________________________________
%
% ________________________________________________________________________

%% sanity checks

% check number of inputs
narginchk(3,3)

% check number of outputs
nargoutchk(2,2)

% check size of inputs
[C, N] = size(X);
[r_t, ~] = size(tau);

% check consistency of tau and X
if ~ (r_t == C)
    error('Incorrect call to multivariate_phase_space: check dimensions of inputs!')
end
%% initialization

% cell to store reconstructed normalized trajectories for each channel
rnt = cell(1,C); 

%% reconstruct phase space trajectories

% for each channel
for c = 1:C
    % select signal
    sig = X(c,:);
    % compute number of vectors in the reconstructed signal
    N_ps = N - (m-1)*tau(c);
    % compute reconstructed trajectories
    ind = hankel(1:N_ps, N_ps:N);
    rnt{c} = sig(ind(:, 1:tau(c):end));
    % select minimum length between channels
    if c == 1
        l_med = N_ps;
    else
        if l_med > N_ps
            l_med = N_ps;
        end
    end
end

%% compute median trajectory

% inizialization of median trajectory for each dimension
X_med = zeros(l_med,m);
% initialization of accumulation array for each channel
acc = nan(l_med,C);

% for each dimension
for v = 1:m    
    % for each channel
    for c = 1:C
        % add each channel in a single matrix 
        acc(:,c) = rnt{c}(1:l_med,v);
    end    
    % compute the median trajectory
    X_med(:,v) = median(acc,2);
end

end