clear
close all

rng(28449);     % For reproducibility

% Get model parameters and numerical settings
par = getPar();

% Bin width for empirical PCF calculation
dr = 0.02;

% Time vector
t = 0:par.dt:par.tMax;

% Mean-field solution
myODE = @(t, u)( (par.b-par.Mu1)*u - par.Mu2*u^2 );
[~, u] = ode45(myODE, t, par.n0);

% Set up separate parameter structures for short-range competition and
% short-range dispersal cases
parSRC = par;
parSRD = par;
parSRC.sigmaComp = par.sigmaShortRange;
parSRD.sigmaDisp = par.sigmaShortRange;


% Spatial moment dynamics solutions
n = 2*par.xiMax/par.dxi+1;  % 2nd moment u2 will be described by a nxn matrix
IC = par.n0^2*ones(n);  % initial condition for u2  = u1^2 uniformly
[~, Z] = ode45(@(t, y)getSMD(t ,y, par), t, IC(:));    % solve SMD equations for 2nd moment
Zm = reshape(Z(end, :), n, n);              % extract matrix for u2 at end of simulation
gf = Zm((n+1)/2, (n+1)/2:end)/Zm(end, end);    % calculate PCF at end of simulation
z1 = sqrt(Z(:, end));                       % extract time-dependent solution for u1 from the far-field of u2
[~, Z] = ode45(@(t, y)getSMD(t ,y, parSRC), t, IC(:)); % repeated for SRC and SRD cases
Zm = reshape(Z(end, :), n, n);
gfSRC = Zm((n+1)/2, (n+1)/2:end)/Zm(end, end);
z1SRC = sqrt(Z(:, end));
[~, Z] = ode45(@(t, y)getSMD(t ,y, parSRD), t, IC(:));
Zm = reshape(Z(end, :), n, n);
gfSRD = Zm((n+1)/2, (n+1)/2:end)/Zm(end, end);
z1SRD = sqrt(Z(:, end));

% Run IBM
[nAgents, X] = runIBM(t, par);
[nAgentsSRC, XSRC] = runIBM(t, parSRC);
[nAgentsSRD, XSRD] = runIBM(t, parSRD);

% Calculate empirical PCF from IBM results
pcf = getPCF(X, dr);
pcfSRC = getPCF(XSRC, dr);
pcfSRD = getPCF(XSRD, dr);

save('results.mat')


r = dr/2:dr:1-dr/2;
xi = 0:par.dxi:(n-1)/2*par.dxi;
%%
h = figure(1);
h.Position = [  711    45   836   902];
tiledlayout(4, 3, "TileSpacing", "compact")
nexttile([2 3]);
plot(t, nAgents)
hold on
plot(t, nAgentsSRC)
plot(t, nAgentsSRD)
h = gca;
h.ColorOrderIndex = 1;
plot(t, z1, 'HandleVisibility', 'off')
plot(t, z1SRC, 'HandleVisibility', 'off')
plot(t, z1SRD, 'HandleVisibility', 'off')
plot(t, u, 'k--')
ylim([0 inf])
xlabel('time')
ylabel('number of agents')
legend('\sigma_d=0.15, \sigma_b=0.15', '\sigma_d=0.015, \sigma_b=0.15', '\sigma_d=0.15, \sigma_b=0.015', 'mean field', 'Location', 'southeast')
title('(a)')
nexttile;
plot(X(:, 1), X(:, 2), '.')
xlabel('x')
ylabel('y')
title('(b) \sigma_d=0.15, \sigma_b=0.15')
nexttile;
h = gca;
hold on
h.ColorOrderIndex = 2;
plot(XSRC(:, 1), XSRC(:, 2), '.')
xlabel('x')
ylabel('y')
title('(c) \sigma_d=0.015, \sigma_b=0.15')
h.Box='on';
nexttile;
h = gca;
hold on
h.ColorOrderIndex = 3;
plot(XSRD(:, 1), XSRD(:, 2), '.')
xlabel('x')
ylabel('y')
title('(d) \sigma_d=0.15, \sigma_b=0.015')
h.Box='on';
nexttile;
plot(r, pcf, xi , gf)
yline(1, 'k:')
xlim([0 0.4])
ylim([0 8])
xlabel('|\xi|')
ylabel('g(\xi)')
title('(e) \sigma_d=0.15, \sigma_b=0.15')
nexttile;
plot(r, pcfSRC, xi , gfSRC)
yline(1, 'k:')
xlim([0 0.4])
ylim([0 8])
xlabel('|\xi|')
ylabel('g(\xi)')
title('(f) \sigma_d=0.015, \sigma_b=0.15')
nexttile;
plot(r, pcfSRD, xi , gfSRD)
yline(1, 'k:')
xlim([0 0.4])
ylim([0 8])
xlabel('|\xi|')
ylabel('g(\xi)')
title('(g) \sigma_d=0.15, \sigma_b=0.015')
saveas(gcf, 'slm.png')