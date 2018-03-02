function [Fe,Kint]=gauss_quadNQ(n,nne,x_e,y_e,k)

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
        Fe=0;
        Kint=0;
        N=[0 0 0 0];
        D=[k 0;0 k];            %Thermal conductivity matrix. k=5 Wcm^-1C^-1
    for i=1:n
        for j=1:n
            N1=0.25*(1-Xi(j))*(1-Eta(i));       %Shape Functions that change wrt to Xi and Eta for each element. 
            N2=0.25*(1+Xi(j))*(1-Eta(i));
            N3=0.25*(1+Xi(j))*(1+Eta(i));
            N4=0.25*(1-Xi(j))*(1+Eta(i));
            
            N=[N1 N2 N3 N4];

            x=N(1)*x_e(1) + N(2)*x_e(2) + N(3)*x_e(3) + N(4)*x_e(4);        %Coordinates in the isoparametric space
            y=N(1)*y_e(1) + N(2)*y_e(2) + N(3)*y_e(3) + N(4)*y_e(4);

            XY=[x_e',y_e'];                                                 %Nodal position matrix
            GN=0.25*[-(1-Eta(i)) (1-Eta(i)) (1+Eta(i)) -(1+Eta(i));         %Derivative of Shape Function matrix
                     -(1-Xi(j)) -(1+Xi(j)) (1+Xi(j)) (1-Xi(j))];
            J=GN*XY;                                                        %Jacobian
            B=inv(J)*GN;                                                    %B matrix

            Fe = Fe + det(J)*W(i)*W(j)*(N)'*(0.5*(8-x)^2*(4-y)^2);          %Integration of source contribution

            Kint = Kint + det(J)*W(i)*W(j)*B'*D*B;                          %Assembly of elemental K matrix
            
        end
    end
end