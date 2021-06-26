function H_n = Hilbert_n(n)
H_n=zeros(n);
for i=1:n
    for j=1:n
        H_n(i,j)=1/(i+j-1);
    end
end