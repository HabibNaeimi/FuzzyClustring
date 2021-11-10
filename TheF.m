%% Calculating The F in offline Clustring.

function F=TheF(x_c,A,B,Samples,Sigma)
nea = numel(A);
Num = zeros(nea,1);
Den = zeros(nea,1);

for i=1:numel(A)
    Num(i) = A(i).*prod(exp(-(Samples(:,1:end-1)-x_c(i,:)).^2/Sigma^2));
    Den(i) = B(i).*prod(exp(-(Samples(:,1:end-1)-x_c(i,:)).^2/Sigma^2));
end
    
a = sum(Num);
b = sum(Den);
F = a/(b+eps);
end