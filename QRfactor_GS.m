function [Q,R]=QRfactor_GS(A)
% MODIFIEDGRAMSCHMIDT modified Gram-Schmidt orthogonalization
%   [Q,R]=ModifiedGramSchmidt(A); computes the modified Gram-Schmidt
%   orthogonalization of the vectors in the columns of the matrix A
[m,n] = size(A);
Q=zeros(m,n);
R=zeros(m,n);
V=A;
for i=1:n
    R(i,i)=norm(V(:,i),2);
    Q(:,i)=V(:,i)/R(i,i);
    for j=i+1:n
        R(i,j)=Q(:,i)'*V(:,j);
        V(:,j)=V(:,j)-R(i,j)*Q(:,i);
    end
end