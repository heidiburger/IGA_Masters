%%Shape Functions Tutorial
%%Question 4: Fourier's Law
order=5;                    %Sekect 2,3 or 5.
k=0.2;                      %W*m^-1*K^-1
X=[0 2.5 5 7.5 10];         %X position along bar in m

if order==2
    T=[273 286]';
end
if order==3
    T=[273 288 286]';
end    
if order==5;
    T=[273 281 288 288 286]';    %Temperature in K
end 

[x,N,B]=MultiShape(X(1),X(end),length(T));
Ta=N*T(1:length(T));
Flux=k*B*T(1:length(T));

%Shape Functions
figure(1)
plot(x,N)
title('Shape Functions')
xlabel('Bar Position (m)')
ylabel('Weight')
grid on

%Temperature Interpolation
figure(2)
subplot(1,2,1)
plot(x,Ta)
title('Temperature Interpolation')
xlabel('Bar Position (m)')
ylabel('Temperature (K)')
legend('Quadratic')

%Heat Flux Interpolation
subplot(1,2,2)
plot(x,Flux)
title('Heat Flux Interpolation')
xlabel('Bar Position (m)')
ylabel('Flux (W*m^-1)')
legend('Quadratic')
