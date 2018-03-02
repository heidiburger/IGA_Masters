%%FEM Project 2
%%May 2017
%%Heidi Burger
%%
clear
% clc
%%Inputs
x=2;
y=1;
%%
%PRE-PROCESSING
ICA=[3 4 2;2 4 5;2 5 1;1 5 6;4 11 9;9 11 10;4 9 8;4 8 5;5 8 6;6 8 7];
X=[0 0 0 4 4 4 6 (8-sqrt(3)) 7 8 8];
Y=[4 1 0 0 3 4 4 3 (4-sqrt(3)) 2 0];

ne=length(ICA);     %10 elements
nn=length(X);       %11 nodes
BC=[3 4 11];

[Txy,T]=tempTri(x,y,ICA,X,Y,BC);

T_e=zeros(ne,3);
for i=1:ne
    x_e(i,:)=X(ICA(i,:));
    y_e(i,:)=Y(ICA(i,:));
    T_e(i,:)=T(ICA(i,:));
    figure(1)
    hold on
    contTri (x_e(i,:),y_e(i,:),T_e(i,:))
end
for i=1:11
    nNr=text(X(i)+0.1,Y(i)+0.02,num2str(i));   
end
hold off

%%
%Temperature and heat flux along the line AB
n=12; %Number of lines of linspace.
x=linspace(0,6,n);
y=4;
%Temperature
for i=1:length(x)
    [Txy(i),T]=tempTri(x(i),y,ICA,X,Y,BC);
end

figure(2)
plot(x,Txy)
title('Temperature Plot for line AB')
xlabel('Position along line AB (cm)')
ylabel('Temperature (Degrees Celcius)')

flux=zeros(2,length(x));

for i=1:length(x)
    k=x(i);
    if k<4
        [B,M]=shapef_B(x_e(4,:),y_e(4,:));
        flux(:,i)=B*T_e(4,:)';
    else
        [B,M]=shapef_B(x_e(10,:),y_e(10,:));
        flux(:,i)=-k*B*T_e(10,:)';    
    end    
end
    
figure(3)
q=quiver(x,y,flux(1,:),flux(2,:))
q.MaxHeadSize = 0.032;
title('Heat Flux along line AB')
xlabel('Position along line AB (cm)')
ylabel('Heat Flux (W/cm^2)')