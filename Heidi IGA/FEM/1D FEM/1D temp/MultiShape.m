%%Multi-order Shape Functions

function [x,N,B]=MultiShape(x0,L,nen)
% x0=0;     %Starting point of element
% L=10;     %Length of element
% nen=4;    %Number of elemental nodes

xp=linspace(x0,L,nen)';   %Position of elemental nodes
x=linspace(x0,L)';
p2=zeros(100,nen);

for i=1:nen
    M(:,i)=xp.^(i-1);
    p1(:,i)=x.^(i-1);               %For N-matrix
    if i==1
        p2(:,i)=0;                  %For B-matrix
    else
        p2(:,i)=(i-1)*x.^(i-2);     
    end
end

N=p1*M^-1;
B=p2*M^-1;

end