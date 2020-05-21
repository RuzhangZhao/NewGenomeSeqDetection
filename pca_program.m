%Input h: 12-dim natural vectors of all vertexes; x: is the linear combination of the integer
%output t w the vertex of polyhedron
function W=pca_program(h,x)
c=h'*x;
c(1:4)=round(c(1:4));
T=[];
n0=size(h,1);
p=pca(h);  %PCA analysis
%dichotomy to find the vertex of polyhedron
for j=1:12
    flag=1;
    t1=0;
    t2=200;
    while(flag)
    a1=c-(t1+t2)/2*p(:,j);
    lb=zeros(n0,1);
    ub=ones(n0,1);
    Aeq=[ub';h'];
    beq=[1;a1];
    [~,fval,exitflag,~]=linprog(ub,[],[],Aeq,beq,lb,ub); 
    if abs(fval-1)<0.000001&&exitflag==1
        t1=(t1+t2)/2;
    else t2=(t1+t2)/2;
    end
    if abs(t1-t2)<0.000001
        flag=0;
        T=[T,t1];
    end
    end
end
for j=1:12
    t(:,j)=c-T(j)*p(:,j);
end
%the another direction
T=[];
for j=1:12
    flag=1;
    t1=0;
    t2=200;
    while(flag)
    a1=c+(t1+t2)/2*p(:,j);
    lb=zeros(n0,1);
    ub=ones(n0,1);
    Aeq=[ub';h'];
    beq=[1;a1];
    [~,fval,exitflag,~]=linprog(ub,[],[],Aeq,beq,lb,ub); 
    if abs(fval-1)<0.000001&&exitflag==1
        t1=(t1+t2)/2;
    else t2=(t1+t2)/2;
    end
    if abs(t1-t2)<0.000001
        flag=0;
        T=[T,t1];
    end
    end
end
for j=1:12
    w(:,j)=c+T(j)*p(:,j);
end
W=[t,w];

