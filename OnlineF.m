%% Calculating The F in Online Mode.
function[f]=OnlineF(xc,X1,X2,Sigma,A,B)

M = max(size(A));
Num = zeros(1,1);
Den = zeros(1,1);
whilecondition = 1;

while whilecondition
    
    for i=1:M
        Num = Num+A(i)*exp(-((X1-xc(1,i))^2+(X2-xc(2,i)) ^2)/Sigma);
        Den = Den+B(i)*exp(-((X1-xc(1,i))^2+(X2-xc(2,i))^2)/Sigma);
    end

    if (Num) > 1e-20 && (Den) > 1e-20
        whilecondition = 0;
    else
        if Sigma<0.01
            Sigma = 1.1*Sigma;
            Num = 0;
            Den = 0;
        else
            whilecondition = 0;
        end
    end
end

f = Num/Den;
end