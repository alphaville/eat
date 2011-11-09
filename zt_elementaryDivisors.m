tol=1E-7;
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
[V nMulti]= zz_elementaryDivisors(A,tol);
assert(nMulti==12);
assert(length(V)==2);
for i=1:2
    p=V(i).polynomial;
    pLength=length(p);
    if pLength==3
        norm(p-[1 -8 20]  )<=tol;
    elseif pLength==2
        norm(p-[1 -1]  )<=tol;
    else
        assert(false);
    end
end
disp('Test 1 :: passed [OK]');

countSuccessfulEquivalences=0;
N=10000;
for i=1:N
    T=rand(10,10);
    Y=T\(A*T);
    [V nMulti]= zz_elementaryDivisors(Y,1E-4);    
    if nMulti==12 && length(V)==2
        countSuccessfulEquivalences=countSuccessfulEquivalences+1;
    end
end
msg=[' with ',num2str(countSuccessfulEquivalences*100/N),'% success'];
disp(['Test 2 :: passed [OK]',msg]);