# Non-Linear-Finite-element-program-
A finite element code to solve the problem of creep of a thick-walled pipe under internal pressure with visco-elastic behaviour
### Problem Statement 
The problem of creep of a thick-walled pipe under internal pressure p is considered
as sketched . The pressure rises linearly up to its final value pmax and
is then hold until final time. Plain strain (along z direction is 0) conditions are
assumed.
<p align="left">
<img src="https://github.com/sai-karthikeya-vemuri/Non-Linear-Finite-element-program-/blob/master/problem_statement1.png"  />
<img src="https://github.com/sai-karthikeya-vemuri/Non-Linear-Finite-element-program-/blob/master/problem_statement2.png"  />
</p>


### User Manual
All the files(.m) have to be saved in a single folder and opened in GNU Octave.Required parameters are given in the array params in the order of
[E[MPa], poisson’s ratio ,Q [MPa], T[s], a[mm], b[mm], pmax[MPa], tL[s],
tf[s]].\
The program after running prints the number of iterations it took for each
time step to converge and radial and cylindrical components of stress and
over stress.

### Testing
Verification is done for Q=0(linear elastic case) and the number of iterations
it took for Newton-Raphson to converge for each time step is 1.The displacements also matched with analytical solution.
