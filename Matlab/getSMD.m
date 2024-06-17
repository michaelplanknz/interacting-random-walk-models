function du2dt = getSMD(t, u2, par)

l = length(u2);

% convert u2 into a square matrix
U2 = reshape(u2, sqrt(l), sqrt(l));
   
% extract u1 as the far-field of u2 by averaging the 4 corner elements
u1 = mean(mean( sqrt(u2([1 end], [1 end])) ));

% make a padded version of U2 (padded with u1^2) for convolutions
n = size(U2);
nPad = 3*(n-1)+1;
U2Pad = u1^2*ones(nPad);
U2Pad(n:2*n-1, n:2*n-1 ) = U2;

% Set up matrices of x and y coordinates
x = par.dxi*(-(n-1)/2:(n-1)/2 );
y = x;
[X, Y] = meshgrid(x, y);

% Evaluate kernels across the domain
kFunc = @(r2, sd)(1/(2*pi*sd^2) * exp(-r2/(2*sd^2)) );
kComp = kFunc(X.^2+Y.^2, par.sigmaComp);
kDisp = kFunc(X.^2+Y.^2, par.sigmaDisp);

% Calcute the 5 terms in the eqn for du2/dt
termA = par.b*u1*kDisp;
termB = par.b * par.dxi^2*conv2(kDisp, U2Pad, 'same');
termC = -par.Mu1*U2;
termD = -par.Mu2*kComp.*U2;

% For the 5th term use the 4/1/1 power-2 closure and evaluate the integral ariding
% from each term of the closure
I1 = 4/5 * (1/u1) * par.dxi^2 * U2 * sum(sum(kComp.*U2));
I2 = 1/5 * (1/u1) * par.dxi^2 * U2 .* conv2(kComp, U2Pad, 'same'); 
I3 = 1/5 * (1/u1) * par.dxi^2 * conv2(kComp.*U2, U2Pad, 'same');
I4 = -1/5 * u1^3;
termE = -par.Mu2 * (I1+I2+I3+I4);

% Sum terms to get du2/dt
M = 2*(termA+termB+termC+termD+termE);

% Reshape as a column vector
du2dt = M(:);




