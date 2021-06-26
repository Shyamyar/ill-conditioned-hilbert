function [Q,R]=QRfactor_HH(A)
[m,n]=size(A);
R=A; %Start with R=A
Q=eye(m); %Set Q as the identity matrix
for k=1:m-1
    x=zeros(m,1);
    x(k:m,1)=R(k:m,k);
    g=norm(x);
    v=x;
    v(k)=x(k)+g;
    %Orthogonal transformation matrix that eliminates one element
    %below the diagonal of the matrix it is post-multiplying:
    s=norm(v);
    if s~=0, w=v/s;
        u=2*R'*w;
        R=R-w*u'; %Product HR
        Q=Q-2*Q*w*w'; %Product QR
    end
end