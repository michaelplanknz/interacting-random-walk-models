# Random walk models in the life sciences: including local interactions

Code for example results shown in review article "Random walk models in the life sciences: including local interactions" by Plank, Simpson & Baker.

The /Julia/ folder contains Julia code to reproduce the results on lattice-based models shown in Figures 1 and 2.

The /Matlab/ folder contains Matlab code to reproduce the results on the spatial logistic model shown in Figure 4. 

# How to use this repository

To reproduce the results on the spatial logistic model, navigate to the /Matlab/ folder and run the Matlab script *main.m*. This runs simulations of the ABM and solves the spatial moment dynamics equation for $u_2(\xi,t)$ described in Box 4 of the paper, as well as the mean-field model (which is the logisitc growth equation). 

The spatial moments solution is calculated numerically by discretising the integro-differential equation for $u_2(\xi,t)$ onto a two-dimensional grid of values for $\xi$. The third moment is approximated using the 4/1/1 power-2 closure. The solution for the second moment at each grid point is then calculated using the Matlab ODE solver *ode45*. Note that the equation for the first moment $u_1(t)$ is not solved explicitly; instead $u_1(t)$ is calculated from the solution for $u_2(\xi,t)$ in the limit of large $|\xi|$ using the relationship $\lim_{\xi\to\infty} u_2(\xi,t) = u_1(t)^2$.  

The pair correlation function is computed empirically from the ABM and from the spatial moments solution via $g(\xi,t) = u_2(\xi,t)/u_1(t)^2$.

Results are plotted and saved as a .png graphic with filename *slm.png*.



