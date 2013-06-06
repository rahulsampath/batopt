function Amat = formAmat(Bmat, elKmat, M, Ne, h, Ka, Kb, dt, p)

Kmat = formKmat(elKmat, M, Ne, h, Ka, Kb, dt, p);

Amat = Kmat + Bmat;

Amat(1, 1) = Amat(1, 1) + h;
Amat(M, M) = Amat(M, M) - h;

