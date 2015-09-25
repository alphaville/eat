tol = 1e-4;
Z =[
     0     1     1     1     1
     0     0     1     1     1
     0     0     0     1     1
     0     0     0     0     1
     0     0     0     0     0 ];
 [poly evals expo]=zz_minPolAlg(Z);
 assert(norm(zz_minPolExhaust(Z)-poly)<=tol);
 assert(norm(poly-[1 0 0 0 0 0])==0);
 assert(evals==0);
 assert(expo==5);
 assert(norm(polymatrixval(poly,Z),'fro')<1E-10);
 disp('Test 1 :: passed [OK]');
 
 
 Z =[
     0     0     1     1     1
     0     0     0     1     1
     0     0     0     0     1
     0     0     0     0     0
     0     0     0     0     0 ];
 poly=zz_minPolAlg(Z);
 assert(norm(zz_minPolExhaust(Z)-poly)<=tol);
 disp('Test 2 :: passed [OK]');
 
 Z =[
     0     0     1     1     1
     0     0     0     0     1
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0 ];
 poly=zz_minPolAlg(Z);
 assert(norm(zz_minPolExhaust(Z)-poly)<=tol);
 disp('Test 3 :: passed [OK]');
 Z=zeros(10,10);
 [poly evals expo]=zz_minPolAlg(Z);
 assert(norm(poly-[1 0])==0)
 assert(evals==0);
 assert(expo==1);
 assert(norm(polymatrixval(poly,Z),'fro')<1E-10);
 disp('Test 4 :: passed [OK]');
 
 Z=eye(10);
 [poly evals expo]=zz_minPolAlg(Z);
 assert(norm(poly-[1 -1])==0)
 assert(evals==1);
 assert(expo==1);
 assert(norm(polymatrixval(poly,Z),'fro')<1E-10);
 disp('Test 5 :: passed [OK]');
 
 N=5;
 Z=zeros(N);
 for i=1:N-1
     Z(i,i+1)=1;
 end
 Z(N,:)=-(1:N);
 poly=zz_minPolAlg(Z);
 assert(norm(polymatrixval(poly,Z),'fro')<1E-10);
 disp('Test 6 :: passed [OK]');
 
 
 A =[
     4     2     0     0     0     0     0     0     0     0
    -2     4     0     0     0     0     0     0     0     0
     0     0     1     1     0     0     0     0     0     0
     0     0     0     1     1     0     1     0     0     0
     0     0     0     0     1     0     0     1     0     0
     0     0     0     0     0     4     2     0     0     0
     0     0     0     0     0    -2     4     0     0     0
     0     0     0     0     0     0     0     4     2     0
     0     0     0     0     0     0     0    -2     4     0
     0     0     0     0     0     0     0     0     0     1];
 poly=zz_minPolAlg(A);
 assert(norm(polymatrixval(poly,A),'fro')<1E-10);
 assert(norm(zz_minPolExhaust(A)-poly)<=tol);
 disp('Test 7 :: passed [OK]');
 
 %Testing minPolExhaust
 N=100;
 for i=1:N
     F=rand(2,2);
     F=blkdiag(F,F,ones(2,2));
     poly=zz_minPolExhaust(F,1E-4);
 end
 disp('Test 8 :: passed [OK]');