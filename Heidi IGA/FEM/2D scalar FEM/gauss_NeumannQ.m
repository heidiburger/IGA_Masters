function [FeN]=gauss_NeumannQ(n,nne,x_e,y_e)

J=(y_e(1)-y_e(2))/2; %%Ratio between real space and isoparametric space. Scales the integral. 
                    
%Location matrix/Integration points (Xi,Eta)
%Weight factors (W)
switch (n)
    case 1
        Xi=     [0];
        Eta=    [0];
        W=      [2];
    case 2
        Xi=     [-1/sqrt(3), 1/sqrt(3)];
        Eta=    [-1/sqrt(3), 1/sqrt(3)];
        W=      [1, 1];
    case 3
        Xi=     [-sqrt(5)/3, 0, sqrt(5)/3];
        Eta=    [-sqrt(5)/3, 0, sqrt(5)/3];
        W=      [5/9, 8/9, 5/9];        
    case 4
        Xi=     [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116];
        Eta=    Xi;
        W=      [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451];
end 

%Do integral:
        FeN=0;
        N=[0 0 0 0];

    for i=1:n
        Eta(i)=-1;  %%Make sure node numbering of the edge along which you integrate corresponds to the Eta=-1 line .
        N=0.25*[(1-Xi(i))*(1-Eta(i)), (1+Xi(i))*(1-Eta(i)), (1+Xi(i))*(1+Eta(i)), (1-Xi(i))*(1+Eta(i)) ];       %Shape Functions
        
        x=N(1)*x_e(1) + N(2)*x_e(2) + N(3)*x_e(3) + N(4)*x_e(4);            %Coordinates in the isoparametric space
        y=N(1)*y_e(1) + N(2)*y_e(2) + N(3)*y_e(3) + N(4)*y_e(4);
        
        q= [ (x+y) 2*(x+y) ];       %Enter Flux Function manually.
        normal= (1/sqrt(2))*[-1 1];      %Work out normal of edge.
        
        FeN=   FeN + J*W(i)*(N)'*(-10*y);           %Integration for Neumann contribution

    end
end