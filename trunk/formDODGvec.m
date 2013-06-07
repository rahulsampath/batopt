function dOdGvec = formDODGvec(Lvec, dFdG1vec, dFdG2vec, dFdG3vec)

dOdGvec = zeros(3, 1);

dOdGvec(1) = -sum(sum(Lvec.*dFdG1vec));
dOdGvec(2) = -sum(sum(Lvec.*dFdG2vec));
dOdGvec(3) = -sum(sum(Lvec.*dFdG3vec));


