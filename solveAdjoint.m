function Lvec = solveAdjoint(T, Tstar, Bmat, elKmat, M, N, Ne, h, H, Ka, Kb, dt)

Lvec = zeros(M, N);

Evec = formEvec(T, Tstar, M, N);
Amat = formAmat(Bmat, elKmat, M, Ne, h, H, Ka, Kb, dt, N);
Lvec(:, N) = Amat\Evec;

for p = (N-1):-1:1
    Evec = formEvec(T, Tstar, M, p);
    Amat = formAmat(Bmat, elKmat, M, Ne, h, H, Ka, Kb, dt, p);
    rhs = Evec + (Bmat*Lvec(:, (p + 1)));
    Lvec(:, p) = Amat\rhs;
end

