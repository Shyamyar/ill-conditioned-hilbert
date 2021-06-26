clear
clc
prompt = input('Do you want to provide Order of Hilbert Matrix (Y or N): ', 's');
if prompt == 'Y' || prompt == 'y'
    n = input('Order of Hilbert Matrix: ');
    k = n;
else
    k = [2,3,4,10,20];
end

f = 2;
table = {'n','GE_Err','GS_Err','HH_Err','M_Err','Inv_Err','Div_Err','Condition Number'};
table2 = {'n','GS_norm','HH_norm','M_norm'};

for n = k
    %Creating Hilbert Matrix
    H_n = Hilbert_n(n);
    
    %Creating 'x'vector
    x_cor = zeros(n,1);
    for i=1:n
        x_cor(i)= (-1)^i * (i / (i + 1));
    end
    
    %Determining 'y' vector
    y = H_n * x_cor;
    
    % Solving using Matrix Functions and Expressions
    x_inv = inv(H_n) * y;
    x_div = H_n \ y;

    % QR by Matlab Funtion
    [Q_M,R_M] = qr(H_n);
    b_M = Q_M' * y;
    x_M = backsub(R_M,b_M);
    orthoM = Q_M'*Q_M;
    normM = norm(orthoM-eye(n),'fro');

    % Gaussian Elimination
    x_GE = solveGauss(H_n,y);
    
    % QR by Gram Schmidt
    [Q_GS,R_GS]=QRfactor_GS(H_n);
    b_GS = Q_GS' * y;
    x_GS = backsub(R_GS,b_GS);
    orthoGS = Q_GS'*Q_GS;
    normGS = norm(orthoGS-eye(n),'fro');

    % QR by Householder Method
    [Q_HH,R_HH]=QRfactor_HH(H_n);
    b_HH = Q_HH' * y;
    x_HH = backsub(R_HH,b_HH);
    orthoHH = Q_HH'*Q_HH;
    normHH = norm(orthoHH-eye(n),'fro');
    
    % Error Calculation
    Err_M = norm(x_cor - x_M, inf);
    Err_inv = norm(x_cor - x_inv, inf);
    Err_div = norm(x_cor - x_div, inf);
    Err_GE = norm(x_cor - x_GE, inf);
    Err_GS = norm(x_cor - x_GS, inf);
    Err_HH = norm(x_cor - x_HH, inf);
    
    % Print results if only one input
    if length(k) == 1
        fprintf ('Hilbert Matrix (%d)\n',n)
        disp(H_n)
        x_cor
        x_GE
        x_GS
        x_HH
        x_M;
        x_inv;
        x_div;
    end
    
    % Creating Error Table
    table{f,1} = n;
    table{f,2} = Err_GE;
    table{f,3} = Err_GS;
    table{f,4} = Err_HH;
    table{f,5} = Err_M;
    table{f,6} = Err_inv;
    table{f,7} = Err_div;
    table{f,8} = cond(H_n);
    
    table2{f,1} = n;
    table2{f,2} = normGS;
    table2{f,3} = normHH;
    table2{f,4} = normM;
   
    f = f + 1;
end

% Display Error Tables
disp(table)
disp(table2)

% Plot Error Graphs
if length(k)>1
    figure(1);
    hold on;
    plot(k,log10(cell2mat(table(2:end,2))))
    plot(k,log10(cell2mat(table(2:end,3))))
    plot(k,log10(cell2mat(table(2:end,4))))
    plot(k,log10(cell2mat(table(2:end,5))))
    plot(k,log10(cell2mat(table(2:end,6))))
    plot(k,log10(cell2mat(table(2:end,7))))
    legend("GE","GS","HH","M","Inv","Div", 'FontSize', 10)
    xlabel("n", 'FontSize', 12);
    ylabel("Log10 Error: log(||x*-x||)", 'FontSize', 12);
    title("Error Vs Order of Hilbert Matrix (n) on Log10 Scale", 'FontSize', 14);
    hold off;

    figure(2);
    plot(k,log10(cell2mat(table(2:end,8))))
    xlabel("n", 'FontSize', 12);
    ylabel("Log 10 Condition Number", 'FontSize', 12);
    title("Condtion Number Vs. Order of Hilbert Matrix (n) on Log10 Scale", 'FontSize', 14);
    
    figure(3);
    hold on;
    plot(k,log10(cell2mat(table2(2:end,2))))
    plot(k,log10(cell2mat(table2(2:end,3))))
    plot(k,log10(cell2mat(table2(2:end,4))))
    legend("GS","HH","M", 'FontSize', 12)
    xlabel("n", 'FontSize', 12);
    ylabel("Log10 Frobenius Norm", 'FontSize', 10);
    title("Frobenius Norm for Q'Q=I on Log10 Sacle", 'FontSize', 14);
    hold off;

end