function [Mean,CV,Corre]=RunODESimulation4(tt, M,ParaBestFit,OtherPara)

cutoff=0.000001;
%disp(ParaBestFit);
options = odeset('RelTol',1e-3,'AbsTol',1);%,'Stats','on');%,'OutputFcn',@odeplot);
[T,solution] = ode45(@ODESimulation, tt, OtherPara.InitialValue(1:M));%,options);
solution=solution';
for ii=1:M
Mean(ii,:)=solution(ii,:);
Var(ii,:)=Mean(1,:);%solution(ii+M,:)-Mean(ii,:).^2;
end
CV=Mean;%./max(Mean,cutoff);

kkk=0;
for ii=1:M
    for jj=ii+1:M
    kkk=kkk+1;      
    Corre(kkk,:)=Mean(1,:);%(solution(kkk+2*M,:)-Mean(ii,:).*Mean(jj,:))./max(sqrt(Var(ii,:)),cutoff)./max(sqrt(Var(jj,:)),cutoff);
    end
end


%%
function ODE = ODESimulation(TimePoint, x)
%Variable vector
M=OtherPara.M;

l=1;kTransition=zeros(M,M);lambda=zeros(1,M);
for i=1:M
    if OtherPara.lambdaMatrix(i)~=0
        lambda(i)=ParaBestFit(l);%
        l=l+1;
    end
end
for i=1:M
    if OtherPara.muMatrix(i)~=0
        mu(i)=ParaBestFit(l);%
        l=l+1;
    end
end
for i=1:M
    for j=1:M
        if OtherPara.Connectivity(i,j)~=0
        kTransition(i,j)=ParaBestFit(l);%
        l=l+1;
        end
    end
    
end


NFit=ParaBestFit(end-5);



for i=1:size(kTransition,2)   
K(i)=lambda(i)-mu(i)-(sum(kTransition(i,:),2)-kTransition(i,i));
KPlus(i)=lambda(i)+mu(i)+(sum(kTransition(i,:),2)-kTransition(i,i));
end


%ODEs generation capacity new 2:
for i=1
    G=log2(x(i)./OtherPara.InitialValue(1))*lambda(i)/(lambda(i)-kTransition(i,2));
ODE(i)=(lambda(i))*x(i)*(1-G/NFit)-mu(i)*x(i)+(-(sum(kTransition(i,:),2)-kTransition(i,i)))*x(i)+(kTransition(:,i)'*x(1:M)-kTransition(i,i)*x(i));
end
for i=2:M
ODE(i)=(lambda(i))*x(i)-mu(i)*x(i)+(-(sum(kTransition(i,:),2)-kTransition(i,i)))*x(i)+(kTransition(:,i)'*x(1:M)-kTransition(i,i)*x(i));
end

disp('ODEs carryign capacity new 2');

ODE=ODE';

end
end
