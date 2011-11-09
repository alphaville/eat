function ind = zz_eigenIndex(A,lambda)
%zz_eigenIndex
%
% Returns the index of an eigenvalue.

ind=0;
n=size(A,1);
M=A-lambda*eye(n);
k=rank(M);
S=M;
for i=1:n    
    flag=rank(S*M)-rank(S); % rank[M^(i+1)]-rank[M^i]
    S = S*M; % S=M^(i+1)
    if flag==0
        ind=i;
        break;
    end
end