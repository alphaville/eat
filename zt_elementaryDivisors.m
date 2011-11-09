tol=1E-7;
A =[
     1     2     0     0     0     0     0     0     0     0
    -2     4     0     0     0     0     0     0     0     0
     0     0     1     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0
     0     0     0     0     1     0     0     0     0     0
     0     0     0     0     0     1     2     0     0     0
     0     0     0     0     0    -2     4     0     0     0
     0     0     0     0     0     0     0     1     2     0
     0     0     0     0     0     0     0    -2     4     0
     0     0     0     0     0     0     0     0     0     1];
 [V nMulti]= zz_elementaryDivisors(A,tol);
 assert(nMulti==12);
 assert(length(V)==2);
 for i=1:2
     if (V(i).polynomial)
         % ***finish test
 end
 