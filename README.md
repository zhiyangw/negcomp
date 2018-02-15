# negcomp one block case 
Key Functions: <br />
plotNewtonPaths.m shows how to use the runSimulation function to plot different Newton Paths with different initial guess <br />
runSimulation.m takes the previous time step solution, first guess and prop object to run a FIM simulation <br />
	also takes in chop fraction and the maximum number of newton steps <br />
initialize.m initializes the different properties for the reservoir and simulation, outputs the prop object that is used throughout <br />
R.m - Residual for mass and energy equation <br />
J.m - Jacobian for residual <br />
