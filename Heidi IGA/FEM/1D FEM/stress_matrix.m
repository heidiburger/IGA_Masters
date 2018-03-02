function [N,B]=stress_matrix(a,b,type,E,A,d,nel)

X_p=linspace(a,b,type);    %Nodal positions
nn=size(d);

    x=linspace(a,b)';
    p1=zeros(100,type);
    p2=zeros(100,type);
    
    for i=1:type
        M(:,i)=X_p.^(i-1);
        p1(:,i)=x.^(i-1);               %For N-matrix
        if i==1
            p2(:,i)=0;                  %For B-matrix
        else
            p2(:,i)=(i-1)*x.^(i-2);     
        end
    end
    iM=M^-1;
    N=p1*iM;
    B=p2*iM;
    

 

