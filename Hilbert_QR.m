clear;
clc;
n=input('Order of Hilbert Matrices (Provide Array of Values): ');
s=length(n);
Hilbert(s)=struct();
f = 1;
for k=n
    Hilbert(f).H_n=zeros(k);
    Hilbert(f).Q_GS=zeros(k);
    Hilbert(f).R_GS=zeros(k);
    Hilbert(f).Q_HH=zeros(k);
    Hilbert(f).R_HH=zeros(k);
    for i=1:k
        for j=1:k
            Hilbert(f).H_n(i,j)=1/(i+j-1);
        end
    end
    [Hilbert(f).Q_GS,Hilbert(f).R_GS] = QRfactor_GS(Hilbert(f).H_n);
    [Hilbert(f).Q_HH,Hilbert(f).R_HH] = QRfactor_HH(Hilbert(f).H_n);
    f = f + 1;
end
%{
for k=n
    fprintf ('Hilbert Matrix (%d)\n',k)
    disp(Hilbert(k).H_n)
    [Q(k).Q_n,R(k).R_n]=qr(Hilbert(k).H_n);
    fprintf ('Q(%d)\n',k)
    disp(Q(k).Q_n)
    fprintf ('R(%d)\n',k)
    disp(R(k).R_n)
end
%}