%%FEM Project 2
%%May 2017
%%Heidi Burger
%%
%clear
% clc
function [T]=tempTri(x,y,ICA,X,Y,BC)
%Inputs
k=5;            %Wcm^-1C^-1
offset=0.1;     %Offset of nodal labels
BC=[3 4 11];     %Dirichlet nodes

ICA=[2 3 4 5;1 2 5 6;4 11 10 9;4 9 8 5;5 8 7 6];
X=[0 0 0 4 4 4 6 (8-sqrt(3)) 7 8 8];
Y=[4 1 0 0 3 4 4 3 (4-sqrt(3)) 2 0];
x=2;
y=1;

ne=length(ICA);     %Number of elements
nn=length(X);       %Number of nodes

%Plot Mesh
for i=1:ne
    x_e(i,:)=X(ICA(i,:));
    y_e(i,:)=Y(ICA(i,:));
    x_d=[x_e(i,:) X(ICA(i,1))];
    y_d=[y_e(i,:) Y(ICA(i,1))];
    figure(1)
    hold on 
    plot(x_d,y_d)
    title('Quadrilatrial Mesh')
end
for i=1:nn
    nNr=text(X(i)+offset,Y(i)+0.02,num2str(i));   
end
hold off

%%
%ASSEMBLY
%For GSM and Force vector:
Fe=zeros(ne,4);
FeN=zeros(ne,4);
F=zeros(nn,1);
K=zeros(nn);
d=zeros(nn,1);
for i=1:ne 
    x_e(i,:)=X(ICA(i,:));       %Elemental nodal coordinates
    y_e(i,:)=Y(ICA(i,:));
    
    if  i==1 || i==2
        FeN(i,:)=gauss_NeumannQ(4,i,x_e(i,:),y_e(i,:))';        %Calculation of Neumann line integral boundary condition
    end
    
    [Fe(i,:),Kint]= gauss_quadNQ(4,i,x_e(i,:),y_e(i,:),k);      %Calculation of source contribution
    D=[k 0;0 k];                                                %Thermal conductivity matrix
    K(ICA(i,:),ICA(i,:))=K(ICA(i,:),ICA(i,:)) + Kint;           %Assembly of Global Stiffness Matrix (GSM)
    F(ICA(i,:))=F(ICA(i,:)) - FeN(i,:)'+ Fe(i,:)';              %Assembly of Force Matrix (equivalent)
end

%%
%%SOLVING
%%Penalty Method 
% for j=1:length(BC)
%     beta=K(BC(j),BC(j))*10^7
%     K(BC(j),BC(j))=beta;
%     d(BC(j))=0;                 
%     F(BC(j))=beta*d(BC(j));
%     d=inv(K)*F;
% end

%Partition Method
Ff=F;
Ff(BC)=[];
Kf=K;           %Exclude the dirichlet nodes
Kf(BC,:)=[];
Kf(:,BC)=[];

d=inv(Kf)*Ff;   %Calculate the temperature at the Neumann nodes

for j=1:length(BC);
    k=BC(j)-1;
    d = [d(1:k);0;d(k+1:end)];      %Add the Dirichlet nodes back into the Temperature vector.
end
%%
%POST PROCESSING Not needed for this part of the project. 
T=d';
% % N=zeros(size(ICA));
% % for i=1:ne
% %     N1 = (1/4)*(1-Xi).*(1-Eta);
% %     N2 = (1/4)*(1+Xi).*(1-Eta);
% %     N3 = (1/4)*(1+Xi).*(1+Eta);
% %     N4 = (1/4)*(1-Xi).*(1+Eta);
% %     
% %     [N(i,:)]=[N1 N2 N3 N4];
% % end
% % 
% % %Work out in which element the point is:
% % N1=-1;
% % N2=-1;
% % N3=-1;
% % N4=-1;
% % i=1;
% % while (N1<0 || N2<0 || N3<0 || N4<0) && i<(ne+1)
% % %     i=i+1;
% %     N1=N(i,1);
% %     N2=N(i,2);
% %     N3=N(i,3);
% %     N3=N(i,4);
% %     i=i+1;
% % end
% % if i>ne
% %     'Error: point outside Domain'
% %     Txy=0;
% % else
% %     el=i-1;
% % 
% %     T_e=T(ICA(el,:));
% %     Txy=N1*T_e(1)+N2*T_e(2)+N3*T_e(3)+N4*T_e(4);  
% % end
end

