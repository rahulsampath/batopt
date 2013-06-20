
clc;
clear all;

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

%Specific Heat (includes density factor)
cp = [1.0; 1.0; 1.0];

%Diffusion Coefficient (Linear Model: K = Ka + Kb*t)
Ka = [10; 10; 10];
Kb = [2; 3; 4]; %MMS solution

elMmat = elemMmat();
elKmat = elemKmat();

%Mass matrix
Bmat = formBmat(elMmat, M, Ne, h, dt, cp);

T = solveForward(Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, Tinf, Tinit, dt);

% %Observations: Taken at x = 0 and all time > 0.
% Tstar = zeros(1, N);
% for p = 1:N
%     Tstar(p) = 300 + (p*2.0);
% end

Tstar = T(1, :);


%Initial guess
Kb = [0; 0; 0];

%Approximate Inverse of Hessian Matrix
invHess = eye(3, 3);

%Maximum optimization iterations
maxOptIters = 100;

dfdG1vec = zeros(M, N);
dfdG2vec = zeros(M, N);
dfdG3vec = zeros(M, N);
% STEP 2 for computing T at all time steps
T = solveForward(Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, Tinf, Tinit, dt);
deltaT = T(1, :) - Tstar;
obj = sum(deltaT.^2);
if obj <= 1.0e-12
   disp('Objective function value is below the tolerance!');
else
    % STEP 3 solve Adjoint for Lambda, Lvec
   Lvec = solveAdjoint(T, Tstar, Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, dt);
   for p = 1:N
      dfdG1vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 1);
      dfdG2vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 2);
      dfdG3vec(:, p) = formDFDGvec(elKmat, M, Ne, h, T, dt, p, 3);
   end
   % STEP 4 compute gradient DODg (dOdg is zero in this case)
   dOdGvec = formDODGvec(Lvec, dfdG1vec, dfdG2vec, dfdG3vec);
   gradOnorm = sqrt((dOdGvec')*dOdGvec);
   if gradOnorm <= 1.0e-12
      disp('Norm of the gradient is below the tolerance!') 
   else   
     % BFGS Iteration 
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
          dOdGvecNew = formDODGvec(Lvec, dfdG1vec, dfdG2vec, dfdG3vec);
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
 
 
 

