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
[DEval sizeDistinct]=zz_distinctEvals(A,tol);
assert(sizeDistinct==3);
flag1=0;flag2=0;flag3=0;
for i=1:3
    if abs(DEval(i,1)-(2.5000+1.3229i))^2<=tol
        flag1=1;
        assert(DEval(i,2)==3);
    end
    if abs(DEval(i,1)-(2.5000-1.3229i))^2<=tol
        flag2=1;
        assert(DEval(i,2)==3);
    end
    if abs(DEval(i,1)-1)^2<=tol
        flag3=1;
        assert(DEval(i,2)==4);
    end
end
assert(flag1==1);assert(flag2==1);assert(flag3==1);
disp('Test 1 :: passed [OK]');

A=diag(ones(20,1));
[DEval sizeDistinct]=zz_distinctEvals(A,tol);
assert(sizeDistinct==1);
assert(DEval(1,1)==1);
disp('Test 2 :: passed [OK]');