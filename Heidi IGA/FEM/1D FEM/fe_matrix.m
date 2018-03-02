function [fe]=fe_matrix(a,b,L,type,E,A,iM)

X_p=linspace(a,b,type);

for j=1:type
%Gauss Quadrature
n=6;        %Gauss point quadrature order
%Location matrix/Integration points (Xi)
Xi=[0.9324695142 0.6612093865 0.2386191861 -0.2386191861 -0.6612093865 -0.9324695142]';
%Weight factors (W)
W=[0.1713244924 0.3607615730 0.4679139346 0.4679139346 0.3607615730 0.1713244924]';

J=0.5*(b-a);
I_hat=0;
for i=1:n
    x=0.5*(b+a)+0.5*Xi(i)*(b-a);
    I_hat=I_hat+W(i)*(x^(j-1)*sin(x*(pi()/L)));
end
I(j)=J*I_hat;
end

fe=iM'*I';
