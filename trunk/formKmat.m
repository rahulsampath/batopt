%Stiffness Matrix
function Kmat = formKmat(elKmat, M, Ne, h, Ka, Kb, dt, p)

K = Ka + (Kb*p*dt);

Kmat = zeros(M, M);

nd0 = 1;
for i = 1:3
    scaling = 2.0*K(i)/h;
    for e = 1:(Ne(i))
        nd1 = nd0 + 1;
        Kmat(nd0, nd0) = Kmat(nd0, nd0) + (scaling*elKmat(1, 1));
        Kmat(nd0, nd1) = Kmat(nd0, nd1) + (scaling*elKmat(1, 2));
        Kmat(nd1, nd0) = Kmat(nd1, nd0) + (scaling*elKmat(2, 1));
        Kmat(nd1, nd1) = Kmat(nd1, nd1) + (scaling*elKmat(2, 2));
        nd0 = nd1;
    end
end



