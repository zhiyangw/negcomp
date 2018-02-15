# negcomp one block case 
Key Functions:
plotNewtonPaths.m shows how to use the runSimulation function to plot different Newton Paths with different initial guess
runSimulation.m takes the previous time step solution, first guess and prop object to run a FIM simulation
	also takes in chop fraction and the maximum number of newton steps
initialize.m initializes the different properties for the reservoir and simulation, outputs the prop object that is used throughout
R.m - Residual for mass and energy equation
J.m - Jacobian for residual
