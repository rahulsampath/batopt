function Evec = formEvec(T, Tstar, M, p)

Evec = zeros(M, 1);

Evec(1) = 2.0*(T(1, p) - Tstar(p));


