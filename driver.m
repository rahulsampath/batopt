
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
Tinit = 300*ones(M, 1);

%Observations: Taken at x = 0 and all time > 0.
Tstar = zeros(1, N);
for p = 1:N
    Tstar(p) = 300 + (p*2.0);
end

%Specific Heat (includes density factor)
cp = [1.0; 1.0; 1.0];

%Diffusion Coefficient (Linear Model: K = Ka + Kb*t)
Ka = [10; 10; 10];
Kb = [0; 0; 0];

elMmat = elemMmat();
elKmat = elemKmat();

%Mass matrix
Bmat = formBmat(elMmat, M, Ne, h, dt, cp);

%Approximate Inverse of Hessian Matrix
invHess = eye(3, 3);

%Maximum optimization iterations
maxOptIters = 100;

dfdG1vec = zeros(M, N);
dfdG2vec = zeros(M, N);
dfdG3vec = zeros(M, N);
T = solveForward(Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, Tinf, Tinit, dt);
deltaT = T(1, :) - Tstar;
obj = sum(deltaT.^2);
if obj <= 1.0e-12
   disp('Objective function value is below the tolerance!');
else
   Lvec = solveAdjoint(T, Tstar, Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, dt);
   for p = 1:N
      dfdG1vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 1);
      dfdG2vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 2);
      dfdG3vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 3);
   end
   dOdGvec = formDODGvec(Lvec, dFdG1vec, dFdG2vec, dFdG3vec);
   gradOnorm = sqrt((dOdGvec')*dOdGvec);
   if gradOnorm <= 1.0e-12
      disp('Norm of the gradient is below the tolerance!') 
   else   
      for iter = 1:maxOptIters
          dG = -invHess*dOdGvec;
          alpha = 1.0;
          KbNew = Kb + dG;
          T = solveForward(Bmat, elKmat, M, N, Ne, h, H, Ka, KbNew, Tinf, Tinit, dt);
          deltaT = T(1, :) - Tstar;
          objNew = sum(deltaT.^2);
          if objNew >= obj
             while alpha > 1.0e-12
                   alpha = 0.5*alpha;
                   dG = 0.5*dG;
                   KbNew = Kb + dG;
                   T = solveForward(Bmat, elKmat, M, N, Ne, h, H, Ka, KbNew, Tinf, Tinit, dt);
                   deltaT = T(1, :) - Tstar;
                   objNew = sum(deltaT.^2);
                   if objNew < obj
                      break;
                   end
             end
             if objNew >= obj
                disp('Line Search Failed!')
                break;
             end
          end
          Kb = KbNew;
          obj = objNew;
          if obj <= 1.0e-12
             disp('Objective function value is below the tolerance!');
             break;
          end
          Lvec = solveAdjoint(T, Tstar, Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, dt);
          for p = 1:N
              dfdG1vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 1);
              dfdG2vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 2);
              dfdG3vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 3);
          end
          dOdGvecNew = formDODGvec(Lvec, dFdG1vec, dFdG2vec, dFdG3vec);
          yVec = dOdGvecNew - dOdGvec;
          dOdGvec = dOdGvecNew;
          gradOnorm = sqrt((dOdGvec')*dOdGvec);
          if gradOnorm <= 1.0e-12
             disp('Norm of the gradient is below the tolerance!') 
             break;
          end
          rho = 1.0/((yVec') * dG);
          %BFGS update
          invHess = ((eye(3,3) - (rho*dG*(yVec')))*invHess*(eye(3,3) - (rho*yVec*(dG')))) + (rho*dG*(dG'));
      end
   end
end
 
 
 

