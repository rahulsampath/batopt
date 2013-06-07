
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

%Convective Coefficient
H = 1.0;

%Bulk Temperature
Tinf = 300;

%Initial Temperature
I = 300*ones(M, 1);

%Observations: Taken at x = 0 and all time > 0.
Tstar = zeros(N, 1);
for p = 1:N
    Tstar(p) = 300 + (p*2.0);
end

%Specific Heat (includes density factor)
cp = [1.0; 1.0; 1.0];

%Diffusion Coefficient (Linear Model: K = Ka + Kb*t)
Ka = [10; 10; 10];
Kb = [0; 0; 0];

maxOptIters = 100;

elMmat = elemMmat();
elKmat = elemKmat();

Bmat = formBmat(elMmat, M, Ne, h, dt, cp);

dfdG1vec = zeros(M, N);
dfdG2vec = zeros(M, N);
dfdG3vec = zeros(M, N);
for iter = 1:maxOptIters
   T = solveForward(Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, Tinf, I, dt);
   Lvec = solveAdjoint(T, Tstar, Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, dt);
   for p = 1:N
      dfdG1vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 1);
      dfdG2vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 2);
      dfdG3vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 3);
   end
   dOdGvec = formDODGvec(Lvec, dFdG1vec, dFdG2vec, dFdG3vec);
   %4. BFGS update
end

%2. Loop (Max Optimization Iterations or Gradient is zero)
	%1. New search direction p = - H (DO/Dg)
        %2. gNew = gOld + alpha p (line search) 
        %3. dg = gNew - gOld
        %4. y = newGradient - oldGradient 
        %5. rho = 1/(y dot dg)
        %6. update H using BFGS formula


 
 
 

