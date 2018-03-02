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
    nel=20;     %Number of elements
    type=3;     %Type of element: 2-linear 3-quadratic
    if type==2
        nn=nel+1;   %Number of nodes for linear element
    else
        nn=2*nel+1; %Number of nodes for quadratic element
    end
    BC=[1 nn];       %Number of Dirichet boundary conditions
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

%%Plot displacements
    x=linspace(X0,L,nn);
    y=-1/6*x.^3 + 1/6*x;%(1/(A*E))*((L/pi())^2*sin(pi()*x/L)+x*(P*A+L/pi()));
    figure(1)
    subplot(2,1,1)

    plot(x,d,'b',x,y,'r')
    title('Plot of displacement against bar length')
    xlabel('Bar length in m')
    ylabel('Displacement in m')

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
        figure(1)
        subplot(2,1,2)
        hold on
        plot(x,stress_nel)
        title('Plot of stress against bar length')
        xlabel('Bar length in m')
        ylabel('Stress in N/m^2')
    end
    hold off
%%
%%Energy Error Method
    for j=1:nel
            if type==2
                asm=[j j+1];
        else if type==3
                asm=[2*j-1 2*i 2*j+1];
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
                I_hat=I_hat+W(i,n)*((L/pi())^2*sin(pi()*x(i)/L)+x(i)*(P*A+L/pi()))^2;
            end
            err_bottom=(E/2)*J*I_hat;
    %Gauss Quadrature to solve for ||uex-uh||en, top of err function.        
            for i=1:n
                x(i)=0.5*(b-a)+0.5*Xi(i,n)*(b+a);
            end
            J=0.5*(b+a);
            I_hat=0;
            for i=1:n
                if type==2
                    I_hat=I_hat+W(i,n)*((L/pi())^2*sin(pi()*x(i)/L)+x(i)*(P*A+L/pi())-((1/(b-a))*[-1 1]*d(asm)))^2;
                else
                    [ke,iM]=ke_matrix(a,b,type,E,A);
                    I_hat=I_hat+W(i,n)*((L/pi())^2*sin(pi()*x(i)/L)+x(i)*(P*A+L/pi())-([0 1 2*x(i)]*iM*d(asm)))^2;
                end
            end
            err_top=(E/2)*J*I_hat;

        err_en(j)=sqrt(err_top)/sqrt(err_bottom);

    end

    err_tot=sum(err_en);
