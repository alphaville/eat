%CARTPROD is a combination generator using multiple sets. Providing a
%pair of sets S={s_i} and T={t_i} produced the set S x T = {(s_i ,
%t_j)}_i,j, that is their cartesian product.
%
%   Syntax PROD = cartprod(X)
%   where X is a cell each entry of which is a cell.
%
%   Example of use. Let A,B,C,D be matrices and S,T are the sets S={A,B}
%   and T = {C,D}. Let X be the cell: X = {S,T}. The following function
%   call:
%                           PROD = cartprod(X)
%   will return the cell PROD with elements PROD{1}={A,C}, PROD{2}={A,D},
%   PROD{3}={B,C} and PROD{4}={B,D}, that is the cartesian product of S and
%   T as a cell. Of course the elements in S and T can be vectors of
%   numbers and the wrapper cell X may contain more than two subcells.
%   Here is a running example:
%   A = [1 2;3 4];
%   B = [1;1]
%   C = 100;
%   D = 'word';
%   X = {{A,B},{C,D}};
%   R = cartprod(X);
%
%   Then R{i}{j} represents the j-th coordinate of the i-th element of the
%   cartesian product. The pairs contained in R will therefore be:
%   (A,C),(A,D),(B,C),(B,D)
%
function R = cartprod(X)
if (~iscell(X))
    error('Input argument is not a cell. Type "help cartprod" for details');
end
nX = length(X);
if nX==1
    lX=length(X{1});
  R=cell(1,lX);
  for i=1:lX
      R{i}=X{1}{i};
  end
else
    R = cartesian(X{1},X{2});
    if (nX > 2)
        for i=3:nX
            R = cartesian(R,X{i});
        end
    end
end

function LR=cartesian(L,R)
if (nargin~=2)
    error('Illegal number of input arguments');
end
if (~iscell(L) || ~iscell(R))
    error('Illegal argument: Both arguments must be cells');
end
nL = length(L);
nR = length(R);

LR = cell(1,nL*nR);
k = 1;
for i=1:nL
    for j=1:nR
        if (iscell(L{i}))
            for p=1:length(L{i})
                LR{k}{p} = L{i}{p};
            end
            LR{k}{p+1}=R{j};
        else
            LR{k} = {L{i}, R{j}};
        end
        k= k +1;
    end
end


