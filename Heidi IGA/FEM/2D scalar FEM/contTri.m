%%MainTri Contour Plot
function []=contTri (x_e,y_e,T_e)
% x_e = [0 1 0];
% y_e = [0 0 1];
% T_e = [1 2 3];

res = 5;
nl = 20;            %Number of lines in contour plot

xi = linspace(0,1,res);
eta = linspace(0,1,res);
[Xi,Eta] =meshgrid(xi,eta);

for i=1:size(Eta)
    Eta(:,i)=Eta(:,i)*(1-Xi(1,i));
end 

N1 = Xi;            %Isoparametric shape Functions
N2 = Eta;
N3 = 1-Xi-Eta;

X=N1*x_e(1)+N2*x_e(2)+N3*x_e(3);
Y=N1*y_e(1)+N2*y_e(2)+N3*y_e(3);
T=N1*T_e(1)+N2*T_e(2)+N3*T_e(3);

contour(X,Y,T,nl)           %Contour plot

for i=1:length(x_e)
    nodeTemp = text(x_e(i)+0.1,y_e(i)-0.1,num2str(round(T_e(i),2)));        %Label temperature values at nodes
    nodeTemp.Color = 'r';
end 
end