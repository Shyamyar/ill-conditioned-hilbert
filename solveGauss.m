function x = solveGauss(A,b)
s = length(A);
U = A;
c = b;
for j = 1:(s-1)
    for i = s:-1:j+1
        m = U(i,j)/U(j,j);
        U(i,:) = U(i,:) - m*U(j,:);
        c(i) = c(i) - m*c(j);
    end
end
x = backsub(U,c);