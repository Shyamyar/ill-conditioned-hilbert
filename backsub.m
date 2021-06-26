function x = backsub(U,c)
s=length(U);
x = zeros(s,1);
x(s) = c(s)/U(s,s);               
for i = s-1:-1:1                    
    sum = 0;clc
    
    for j = s:-1:i+1                
        sum = sum + U(i,j)*x(j);    
    end 
    x(i) = (c(i)- sum)/U(i,i);
end 