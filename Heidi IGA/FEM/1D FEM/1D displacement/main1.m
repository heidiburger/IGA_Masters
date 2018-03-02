function [err_tot]=main1(nel,type)
%%FEM Project 1
%%
%%Pre-processing
%Geometry
X0=0;       %Starting position of rod in m
L=1;        %Length of rod in m
A=1;        %Area in m^2
%Properties
E=100;      %Young's Modulus in N/m^2
P=10;       %Traction stress
%Code
%nel=20;     %Number of elements
%type=2;     %Type of element: 2-linear 3-quadratic
if type==2
    nn=nel+1;   %Number of nodes for linear element
else
    nn=2*nel+1; %Number of nodes for quadratic element
end
BC=1;       %Number of Dirichet boundary conditions
le=L/nel;   %Length of element
%%
%%Assembly
%Initialise matrices
K=zeros(nn);
F=zeros(nn,1);
d=zeros(nn,1);

%Build GSM and F
for i=1:nel
    if type==2
            asm=[i i+1];
    else if type==3
            asm=[2*i-1 2*i 2*i+1];
        end
    end
    
    a=le*(i-1);
    b=le*(i);
%Construct K from ke
    [ke,iM]=ke_matrix(a,b,type,E,A);
    K(asm,asm)=K(asm,asm)+ke;
%Construct F from fe
    fe=fe_matrix(a,b,L,type,E,A,iM);
 	F(asm)=F(asm)+fe;       %Body function integral set in ke_matrix function
    
end

    F(end)=F(end)+A*P;
%%
%%Solving
%BC's
d(BC+1:end,1)=nan;          %To ensure that code breaks if I do something wrong.
F(1:BC,1)=nan;
%Partition K
K_e=K(1:BC,1:BC);
K_ff=K(BC+1:end,BC+1:end);
K_ef=K(1,BC+1:end);
%Partition F
F_e=F(1:BC);
F_f=F(BC+1:end);
%Partition d
d_e=d(1:BC);
d_f=d(BC+1:end);

%Solving for free displacements and essential forces
d_f=inv(K_ff)*(F_f-K_ef'*d_e);
F_e=K_e*d_e+K_ef*d_f;

%Construct Force and Displacement vector
K(1:BC,1:BC)=K_e;
K(BC+1:end,BC+1:end)=K_ff;
K(1,BC+1:end)=K_ef;

F(1:BC)=F_e;
F(BC+1:end)=F_f;

d(1:BC)=d_e;
d(BC+1:end)=d_f;

%%Plot displacements
x=linspace(X0,L,nn);
y=(1/(A*E))*((L/pi())^2*sin(pi()*x/L)+x*(P*A+L/pi()));
% figure(1)
% plot(x,d,'b',x,y,'r')
% title('Plot of displacement against bar length: FEM solution (blue) and analytical (red)')
% xlabel('Bar length in m')
% ylabel('Displacement in m')

%Stress
for i=1:nel
    if type==2
            asm=[i i+1];
    else if type==3
            asm=[2*i-1 2*i 2*i+1];
        end
    end
    a=le*(i-1);
    b=le*(i);
    
    [N,B]=stress_matrix(a,b,type,E,A,d,1);
    stress_nel=E*B*d(asm);

    x=linspace(a,b,100);
%     figure(2)
%     clear figure(2)
%     hold on
%     plot(x,stress_nel)
%     title('Plot of stress against bar length')
%     xlabel('Bar length in m')
%     ylabel('Stress in Nm^2')
end
% hold off
%%
%%Error Analysis

% Exact Solution to u(x)
% uex=(1/(A*E))*((L/pi())^2*sin(pi()*x/L)+x*(P*A+L/pi()));

% FEM Solution to Analytical
% uh=Nd

% % err=zeros(1,nel);
% % err_top=0;
% % err_bottom=0;
% % 
% % for j=1:nel
% %     if type==2
% %             asm=[j j+1];
% %     else if type==3
% %             asm=[2*j-1 2*j 2*j+1];
% %         end
% %     end
% %     a=le*(j-1);
% %     b=le*(j);
% %     n=3;        %Gauss point quadrature order
% % %%
% % %Gauss Quadrature to solve for ||uex||L2, bottom of err function.
% % 
% % %Location matrix/Integration points (Xi)
% %         Xi=[0 0 0;1/sqrt(3) -1/sqrt(3) 0;sqrt(0.6) 0 -sqrt(0.6)]';
% % %Weight factors (W)
% %         W=[2 0 0;1 1 0;5/9 8/9 5/9]';
% % 
% %         for i=1:n
% %             x(i)=0.5*(b-a)+0.5*Xi(i,n)*(b+a);
% %         end
% %         J=0.5*(b+a);
% %         I_hat=0;
% %         for i=1:n
% %             I_hat=I_hat+W(i,n)*(1/(A*E))*((L/pi())^2*sin(pi()*x(i)/L)+x(i)*(P*A+L/pi()));
% %         end
% %         err_bottom=err_bottom+J*I_hat;
% % %%
% % %Gauss Quadrature to solve for ||uex-uh||L2, top of err function.
% %         
% % %Location matrix/Integration points (Xi)
% %         Xi=[0 0 0;1/sqrt(3) -1/sqrt(3) 0;sqrt(0.6) 0 -sqrt(0.6)]';
% % %Weight factors (W)
% %         W=[2 0 0;1 1 0;5/9 8/9 5/9]';
% % 
% %         for i=1:n
% %             x(i)=0.5*(b-a)+0.5*Xi(i,n)*(b+a);
% %         end
% %         J=0.5*(b+a);
% %         I_hat=0;
% %         for i=1:n
% %             [ke,iM]=ke_matrix(a,b,type,E,A);
% %             if type==2
% %             I_hat=I_hat+W(i,n)*((1/(A*E))*((L/pi())^2*sin(pi()*x(i)/L)+x(i)*(P*A+L/pi()))-([1 x(i)]*iM*d(asm)))^2;
% %             else
% %             I_hat=I_hat+W(i,n)*((1/(A*E))*((L/pi())^2*sin(pi()*x(i)/L)+x(i)*(P*A+L/pi()))-([1 x(i) x(i)^2]*iM*d(asm)))^2;                
% %             end
% %         end
% %         err_top=err_top+J*I_hat;
% %     
% % 
% % end
% % err_tot=sqrt(err_top)/sqrt(err_bottom);


%Total Error is the sum of all errors for each element.

err_top=0;
err_bottom=0;
for j=1:nel
        if type==2
            asm=[j j+1];
        else if type==3
            asm=[2*j-1 2*j 2*j+1];
        end
        end
    a=le*(j-1);
    b=le*(j);
%For Linear
%uh=(1/(b-a))*[-1 1]*d(asm)
        n=3;        %Gauss point quadrature order
        
%Location matrix/Integration points (Xi)
        Xi=[0 0 0;1/sqrt(3) -1/sqrt(3) 0;sqrt(0.6) 0 -sqrt(0.6)]';
%Weight factors (W)
        W=[2 0 0;1 1 0;5/9 8/9 5/9]';
        
%Gauss Quadrature to solve for ||uex||en, bottom of err function.
        for i=1:n
            x(i)=0.5*(b-a)+0.5*Xi(i,n)*(b+a);
        end
        J=0.5*(b+a);
        I_hat=0;
        for i=1:n
            I_hat=I_hat+W(i,n)*((1/(A*E))*((L/pi)*cos(x(i)*pi/L)+(P/A)+(L/pi)))^2;
            elnum=j
            %check=I_hat*J
        end
        err_bottom=err_bottom+(E/2)*J*I_hat;
%Gauss Quadrature to solve for ||uex-uh||en, top of err function.        
        for i=1:n
            x(i)=0.5*(b-a)+0.5*Xi(i,n)*(b+a);
        end
        J=0.5*(b+a);
        I_hat=0;
        
        for i=1:n
            if type==2
                uex = (1/(E*A))  *  ((L/pi)*cos(x(i)*pi/L) + P*A + (L/pi));
                uh = (1/(b-a))*[-1 1] * d(asm);
                I_hat = I_hat + W(i,n)*(uex-uh)^2;
            else
                [ke,iM]=ke_matrix(a,b,type,E,A);
                
                uex = (1/(E*A))  *  ((L/pi)*cos(x(i)*pi/L) + P*A + (L/pi));
                uh = [0 1 2*x(i)] * iM * d(asm);
                I_hat = I_hat + W(i,n)*(uex-uh)^2;
                
            end
        end
        err_top=err_top+(E/2)*J*I_hat;
        errr=E*0.5*J*I_hat
        
end

err_tot=sqrt(err_top)/sqrt(err_bottom);

%err_tot=sum(err_en)

end