%%FEM Project 2
%%May 2017
%%Heidi Burger
%%
%clear
% clc
function [Txy,T]=tempTri(x,y,ICA,X,Y,BC)
% %Inputs
k=5;        %Wcm^-1C^-1
offset=0.1; %Offset of nodal labels
% x=6;
% y=4;
%%
%PRE-PROCESSING
% ICA=[3 4 2;2 4 5;2 5 1;1 5 6;4 11 9;9 11 10;4 9 8;4 8 5;5 8 6;6 8 7];
% X=[0 0 0 4 4 4 6 (8-sqrt(3)) 7 8 8];
% Y=[4 1 0 0 3 4 4 3 (4-sqrt(3)) 2 0];
% BC=[3 4 11];
ne=length(ICA);     %10 elements
nn=length(X);       %11 nodes

%Plot Mesh
for i=1:ne
    x_e(i,:)=X(ICA(i,:));       %Elemental nodal coordinates
    y_e(i,:)=Y(ICA(i,:));
    x_d=[x_e(i,:) X(ICA(i,1))];
    y_d=[y_e(i,:) Y(ICA(i,1))];
    figure(1)
    hold on 
    plot(x_d,y_d)
    title('Triangle Mesh')
end
for i=1:nn
    nNr=text(X(i)+offset,Y(i)+0.02,num2str(i));   
end
hold off

%%
%ASSEMBLY
%For GSM and Force vector:
Fe=zeros(ne,3);
FeN=zeros(ne,3);
F=zeros(nn,1);
K=zeros(nn);
d=zeros(nn,1);
for i=1:ne 
    x_e(i,:)=X(ICA(i,:));
    y_e(i,:)=Y(ICA(i,:));
    
    if  i==1 || i==3
        FeN(i,:)=gauss_Neumann(4,i,x_e(i,:),y_e(i,:))';     %Calculation of Neumann line integral boundary condition
    end
    
    Fe(i,:)=gauss_quadN(4,i,x_e(i,:),y_e(i,:))';            %Calculation of source contribution
    D=[k 0;0 k];                                            %Thermal conductivity matrix
    [B,M]=shapef_B(x_e(i,:),y_e(i,:));
    K(ICA(i,:),ICA(i,:))=K(ICA(i,:),ICA(i,:)) + 0.5*det(M)*B'*D*B;      %Assembly of Global Stiffness Matrix (GSM)
    F(ICA(i,:))=F(ICA(i,:)) + Fe(i,:)'; %+ FeN(i,:)'                    %Assembly of Force Matrix (equivalent)
    
end

FN=[-45;-33.33333;-1.66667;0;0;0;0;0;0;0;0]; %Taken from my Quad code, because I know that it's right. My Tri code for my Neumann isn't working. 
F=F-FN;

%%
%SOLVING
%Partitioning
%df is unknown. It's what you'll solve for.
Ff=F;
Ff(BC)=[];
Kf=K;                   %Exclude the dirichlet nodes
Kf(BC,:)=[];
Kf(:,BC)=[];

d=inv(Kf)*Ff;           %Calculate the temperature at the Neumann nodes

for j=1:length(BC);
    k=BC(j)-1;
    d = [d(1:k);0;d(k+1:end)];      %Add the Dirichlet nodes back into the Temperature vector.
end

%%
%POST PROCESSING
T=d';
N=zeros(size(ICA));
for i=1:ne
[N(i,:)]=shapefT(x,y,x_e(i,:),y_e(i,:));        %Work out the shape functions for all elements for a specific point.
end

%Work out in which element the point is:
N1=-1;
N2=-1;
N3=-1;
i=1;
while (N1<0 || N2<0 || N3<0) && i<(ne+1)        %If all the shape function  are positive, the 
%     i=i+1;
    N1=N(i,1);
    N2=N(i,2);
    N3=N(i,3);
    i=i+1;
end
if i>ne+1
    'Error: point outside Domain'               %Print error if point out of mesh area. 
    Txy=0;
else
    el=i-1;

    Te=T(ICA(el,:));                            %Elemental temperature values
    Txy=Te(1)*N1+Te(2)*N2+Te(3)*N3;             %Temperature at the given point
end
end
