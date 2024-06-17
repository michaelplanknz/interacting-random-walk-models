function par = getPar()

% Get model parameters and numerical settings

par.b = 0.4;        % birth rate
par.Mu1 = 0.2;      % density independent death rate
par.Mu2 = 0.001;    % neighbour dependent death rate

par.sigmaDisp = 0.15;   % SD of dispersal kernel (default)
par.sigmaComp = 0.15;   % SD of competition kernel (default)
par.sigmaShortRange = 0.015;  % SD of kernel (short-range case)


par.tMax = 100;     % run time
par.n0 = 20;        % initial number of agents
par.nMax = 1e4;     % max number of agents
par.dt = 0.01;      % time step


par.xiMax = 0.5;    % max value of xi variable in SMD
par.dxi = 0.005;    % step size of xi variable in SMD

