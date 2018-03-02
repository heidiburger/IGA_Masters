%%FEM Project 2
%%May 2017
%%Heidi Burger
%%
clear
%Inputs
k=5;             %Wcm^-1C^-1
offset=0.1;      %OFfset of nodal labels
% BC=[3 4 11];     %Dirichlet nodes
% 
% x=2;             %Arbitrary point for analysis
% y=1;
% %%
% %Pre-processing
% ICA=[1 2 3 4];        %Element node numbers
% X=[1 2 3 2];                    %Nodal coordinates
% Y=[1 1 2 2];

BC=[3 4 11];     %Dirichlet nodes

ICA=[2 3 4 5;1 2 5 6;4 11 10 9;4 9 8 5;5 8 7 6];
X=[0 0 0 4 4 4 6 (8-sqrt(3)) 7 8 8];
Y=[4 1 0 0 3 4 4 3 (4-sqrt(3)) 2 0];
x=2;
y=1;

ne=length(ICA);     %5 elements
nn=length(X);       %11 nodes

[T]=tempQuad(x,y,ICA,X,Y,BC);       %Function calculating the nodal temperature values

T_e=zeros(ne,4);
for i=1:ne
    x_e=X(ICA(i,:));
    y_e=Y(ICA(i,:));
    T_e=T(ICA(i,:));
    figure(1)
    hold on
    contQuad(x_e,y_e,T_e)          %Contour plot of temperature values accross domain.
end
for i=1:11
    nNr=text(X(i)+0.1,Y(i)+0.02,num2str(i));   %Label nodes
end
hold off

