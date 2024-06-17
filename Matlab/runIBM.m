function [nAgents, X] = runIBM(t, par)

kFunc = @(r2, sd)(1/(2*pi*sd^2) * exp(-r2/(2*sd^2)) );      % Gaussian kernel function

% X is a n x 2 matrix where n the number of agents and the two columns are
% the x and y coordinates in [0,1] for each agent:
X = rand(par.n0, 2);        % initial agent locations are uniformly random

nSteps = length(t)-1;
dt = t(2)-t(1);

nAgents = zeros(nSteps+1, 1);
nAgents(1) = par.n0;

for iStep = 1:nSteps
    prolFlag = rand(nAgents(iStep), 1) < dt * par.b;    % simulate which agents profilerate
    nDaughters = sum(prolFlag);
    xDaughter = mod(X(prolFlag, :) + par.sigmaDisp*randn(nDaughters, 2), 1);        % simulate daughter locations by adding a draw from the dispersal kernel to the parent location and applicying periodic boundaries

    deltaX = mod(X(:, 1)-X(:, 1)'+0.5, 1) - 0.5;     % dx with periodic boundaries
    deltaY = mod(X(:, 2)-X(:, 2)'+0.5, 1) - 0.5;
    Z = kFunc(deltaX.^2+deltaY.^2, par.sigmaComp);   % evaluate competition kernel at pair distances
    pDeath = dt * (par.Mu1 + par.Mu2*sum(triu(Z, 1)+tril(Z, -1), 2));       % sum to calculate contribution to overall death rate (excluding self pairs by removing diagonla of Z)
    dieFlag = rand(nAgents(iStep), 1) < pDeath;     % simulate which agents die

    X = [X(~dieFlag, :); xDaughter];        % new list of agent locations (excluding dead agents and including new daughters)
    nAgents(iStep+1) = size(X, 1);          % record number of agents for the next time step
end

