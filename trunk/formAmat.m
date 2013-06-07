%System matrix for forward and adjoint equations
function Amat = formAmat(Bmat, elKmat, M, Ne, h, H, Ka, Kb, dt, p)

Kmat = formKmat(elKmat, M, Ne, h, Ka, Kb, dt, p);

Amat = Kmat + Bmat;

%Robin Matrix Correction
Amat(1, 1) = Amat(1, 1) + H;
Amat(M, M) = Amat(M, M) - H;

