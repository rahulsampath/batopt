
%Number of time steps
N = 10; 

%Time step 
dt = 0.1;

%Number of elements in each of the 3 segments.
Ne = [10; 10; 10];

%Total number of grid points
M = 1 + sum(Ne);

%grid spacing
h = 0.1;

cp = [1.0; 1.0; 1.0];

Ka = [10; 10; 10];
Kb = [0; 0; 0];

maxOptIters = 100;

elemMmat = elemMassMat();
elemKmat = elemStiffMat();

Bmat = formBmat(elemMmat, M, Ne, h, dt, cp);

for iter = 1:maxOptIters
end

%2. Loop (Max Optimization Iterations or Gradient is zero)
   %1. Solve Forward Problem
	%1. Time Stepping (Forward in time)
	    %1. Compute A matrix for current time-step
	    %2. Compute C vector for current time-step.
            %3. Solve for Temperatures at current time-step.
   %2. Solve Adjoint Problem
	%1. Time Stepping (Backward in time)	 
	    %1. Compute dO/dT 
            %2. Compute dF/dT
            %3. Form Adjoint matrix
            %4. Solve for Lagrange multipliers (Lambdas)
   %3. Compute DO/Dg (Gradient)
   %4. BFGS update
	%1. New search direction p = - H (DO/Dg)
        %2. gNew = gOld + alpha p (line search) 
        %3. dg = gNew - gOld
        %4. y = newGradient - oldGradient 
        %5. rho = 1/(y dot dg)
        %6. update H using BFGS formula


 
 
 

