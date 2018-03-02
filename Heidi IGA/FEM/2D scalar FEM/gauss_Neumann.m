function [FeN]=gauss_Neumann(n,nne,x_e,y_e)

J=(y_e(3)-y_e(1))/1; %Ratio between real element lenght and the isoparametric element length.

%Location matrix/Integration points (Xi,Eta)
%Weight factors (W)
switch (n)
    case 1
        Xi=     [1/3];
        Eta=    [1/3];
        W=      [1];
    case 3
        Xi=     [1/6, 2/3, 1/6];
        Eta=    [1/6, 1/6, 2/3];
        W=      [1/3, 1/3, 1/3];
    case 4
        Xi=     [1/3, 1/5, 3/5, 1/5];
        Eta=    [1/3, 1/5, 1/5, 3/5];
        W=      [-9/16, 25/48, 25/48, 25/48];
end

%Do integral:
        FeN=0;
        N=[0 0 0];
    for i=1:n
        Eta(i)=0;                               %Let Eta=0 for a line integral problem. Make sure that your nodes are numbered correctly to correspond to the Eta=0 line. 
        N=[Xi(i), Eta(i), 1-Xi(i)-Eta(i)];      %Shape Functions
        x=N(1)*x_e(1) + N(2)*x_e(2) + N(3)*x_e(3);      %Coordinates in the isoparametric space
        y=N(1)*y_e(1) + N(2)*y_e(2) + N(3)*y_e(3);
        
        FeN=   FeN + J*W(i)*(N)'*(-10*y);       %Integration for Neumann condition
    end
end