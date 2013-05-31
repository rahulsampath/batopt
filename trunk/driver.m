
%Assume uniform grid in space and time.

%0. Compute B.
%1. Initial Guess for g (K1b, K2b, K3b)
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


 
 
 

