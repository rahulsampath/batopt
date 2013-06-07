function T = solveForward(Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, Tinf, Tinit, dt)

T = zeros(M, N);

Cvec = formCvec(M, h, H, Tinf, dt, 1);
Amat = formAmat(Bmat, elKmat, M, Ne, h, H, Ka, Kb, dt, 1);
rhs = Cvec + (Bmat*Tinit);
%Linear Solve
T(:, 1) = Amat\rhs; 

for p = 2:N
    Cvec = formCvec(M, h, H, Tinf, dt, p);
    Amat = formAmat(Bmat, elKmat, M, Ne, h, H, Ka, Kb, dt, p);
    rhs = Cvec + (Bmat*T(:, (p - 1))); 
    T(:, p) = Amat\rhs; 
end


