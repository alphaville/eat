Z =[
     0     1     1     1     1
     0     0     1     1     1
     0     0     0     1     1
     0     0     0     0     1
     0     0     0     0     0 ];
 [poly evals expo]=zz_minPolAlg(Z);
 assert(norm(poly-[1 0 0 0 0 0])==0);
 assert(evals==0);
 assert(expo==5);
 assert(norm(polymatrixval(poly,Z),'fro')<1E-10);
 disp('Test 1 :: passed [OK]');
 
 Z=zeros(10,10);
 [poly evals expo]=zz_minPolAlg(Z);
 assert(norm(poly-[1 0])==0)
 assert(evals==0);
 assert(expo==1);
 assert(norm(polymatrixval(poly,Z),'fro')<1E-10);
 disp('Test 2 :: passed [OK]');
 
 Z=eye(10);
 [poly evals expo]=zz_minPolAlg(Z);
 assert(norm(poly-[1 -1])==0)
 assert(evals==1);
 assert(expo==1);
 assert(norm(polymatrixval(poly,Z),'fro')<1E-10);
 disp('Test 3 :: passed [OK]');
 
 N=5;
 Z=zeros(N);
 for i=1:N-1
     Z(i,i+1)=1;
 end
 Z(N,:)=-(1:N);
 poly=zz_minPolAlg(Z);
 assert(norm(polymatrixval(poly,Z),'fro')<1E-10);
 disp('Test 4 :: passed [OK]');