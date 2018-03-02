function [I]=gauss_quadN(n,nne,x_e,y_e)

Ae=     0.5*( x_e(1)*(y_e(2)-y_e(3)) + x_e(2)*(y_e(3)-y_e(1)) + x_e(3)*(y_e(1)-y_e(2)) );   %Area of a triangle
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
        I=0;
        N=[0 0 0];
    for i=1:n
        N=[Xi(i), Eta(i), 1-Xi(i)-Eta(i)];              %Shape Functions
        x=N(1)*x_e(1) + N(2)*x_e(2) + N(3)*x_e(3);      %Scaled xand y coordinates in the parametric space. 
        y=N(1)*y_e(1) + N(2)*y_e(2) + N(3)*y_e(3);
        
        I=   I + Ae*W(i)*(N)'*(0.5*(8-x)^2*(4-y)^2);    %Integration
    end

end