function charpoly = zz_charPolyAlg(A,tol)

if nargin==1
    tol=1E-6;
end
V=zz_elementaryDivisors(A,tol);
charpoly=1;
for i=1:length(V)
    charpoly=conv(charpoly,prv_expandDivisor(V(i)));
end

function poly = prv_expandDivisor(divisor)
    m=divisor.multiplicity;
    p=divisor.polynomial;
    poly=1;
    for i=1:m
        poly=conv(poly,p);
    end
