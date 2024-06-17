function g = getPCF(X, dr)

% Calculate PCF from agent locations X in a unit square, with bin width dr for pair distances
% X is a n x 2 matrix where n the number of agents and the two columns are
% the x and y coordinates in [0,1] for each agent 

% Vector of bin edges for distance between agents in a pair
r = 0:dr:1;

nAgents = size(X, 1);

deltaX = mod(X(:, 1)-X(:, 1)'+0.5, 1) - 0.5;        % matrix of horizontal distances dx between agent i and agent j with periodic boundaries
deltaY = mod(X(:, 2)-X(:, 2)'+0.5, 1) - 0.5;        % matrix of vertical distances dy between agent i and agent j with periodic boundaries
deltaR = sqrt(deltaX.^2 + deltaY.^2);               % matrix of Euclidean distances dr between agent i and agent j with periodic boundaries 
deltaR = deltaR(deltaR > 0);                        % exclude self pairs

% Number of pairs in each distance bin
nPairs = histcounts(deltaR, r);

rMid = 0.5*(r(1:end-1)+r(2:end));   % bin mid points

% Calculate PCF by dividing nPairs by the expected number of pairs in a Poisson process with the same average density
% This is the pair denasity (= agent density squared) times the area of the
% annulus with midpoint radius rMid and width dr 
g = nPairs./(nAgents^2 * 2*pi*dr*rMid); 



