function [ke,iM]=ke_matrix(a,b,type,E,A)

X_p=linspace(a,b,type);    %Nodal positions

M=zeros(type);
if type==2
    p2p=[0 0; 0 b-a];
else if type==3
    p2p=[0 0 0; 0 b-a b^2-a^2; 0 b^2-a^2 (4/3)*(b^3-a^3)];
    end
end
for i=1:type
    M(:,i)=X_p.^(i-1);
end
iM=M^-1;
%Integral
BtB=iM'*p2p*iM;
ke=E*A*BtB;