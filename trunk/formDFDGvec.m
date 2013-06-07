%Derivative of Constraints with respect to control variables
function dFdGvec = formDFDGvec(elKmat, M, Ne, h, T, dt, p, compId)

dFdGvec = zeros(M, 1);

scaling = 2.0*p*dt/h;
tmpMat = zeros((1 + Ne(compId)), (1 + Ne(compId)));
if compId == 1
   stId = 1;
   inVec = T(stId:(stId + Ne(1)), p);
   nd0 = 1;
   for e = 1:(Ne(compId))
       nd1 = nd0 + 1;
       tmpMat(nd0, nd0) = tmpMat(nd0, nd0) + (scaling*elKmat(1, 1));
       tmpMat(nd0, nd1) = tmpMat(nd0, nd1) + (scaling*elKmat(1, 2));
       tmpMat(nd1, nd0) = tmpMat(nd1, nd0) + (scaling*elKmat(2, 1));
       tmpMat(nd1, nd1) = tmpMat(nd1, nd1) + (scaling*elKmat(2, 2));
       nd0 = nd1;
   end
   dFdGvec(stId:(stId + Ne(1))) = tmpMat*inVec;
elseif compId == 2
   stId = 1 + Ne(1);
   inVec = T(stId:(stId + Ne(2)), p);
   nd0 = 1;
   for e = 1:(Ne(compId))
       nd1 = nd0 + 1;
       tmpMat(nd0, nd0) = tmpMat(nd0, nd0) + (scaling*elKmat(1, 1));
       tmpMat(nd0, nd1) = tmpMat(nd0, nd1) + (scaling*elKmat(1, 2));
       tmpMat(nd1, nd0) = tmpMat(nd1, nd0) + (scaling*elKmat(2, 1));
       tmpMat(nd1, nd1) = tmpMat(nd1, nd1) + (scaling*elKmat(2, 2));
       nd0 = nd1;
   end
   dFdGvec(stId:(stId + Ne(2))) = tmpMat*inVec;
else
   stId = 1 + Ne(1) + Ne(2);
   inVec = T(stId:(stId + Ne(3)), p);
   nd0 = 1;
   for e = 1:(Ne(compId))
       nd1 = nd0 + 1;
       tmpMat(nd0, nd0) = tmpMat(nd0, nd0) + (scaling*elKmat(1, 1));
       tmpMat(nd0, nd1) = tmpMat(nd0, nd1) + (scaling*elKmat(1, 2));
       tmpMat(nd1, nd0) = tmpMat(nd1, nd0) + (scaling*elKmat(2, 1));
       tmpMat(nd1, nd1) = tmpMat(nd1, nd1) + (scaling*elKmat(2, 2));
       nd0 = nd1;
   end
   dFdGvec(stId:(stId + Ne(3))) = tmpMat*inVec;
end



