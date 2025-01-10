function  FitValue= ODE_Fitting5(ParaToFit,TimePoint,OtherPara,MaxMean,label)%,OtherPara,DataTimePoint)%%
%Need ODE45:
M=OtherPara.M;%

%display(ParaToFit);
%OtherPara.InitialValue

options = odeset('RelTol',1e-3,'AbsTol',1);%,'NonNegative',[1:5]);%,'Stats','on');%,'OutputFcn',@odeplot);
[T,solution] = ode45(@ODENumericalFormula, TimePoint, OtherPara.InitialValue(1:5));%,options);

%Use Y(i) to represent variables...
%matlab generate an numerical ODE formula...
%%
function ODE = ODENumericalFormula(TimePoint, x)
%Variable vector
%Corre2 = sym('Corre',[M,M]);
% kk=1;
% Corre2=zeros(M,M);
% for j=1:M
%     for i=j+1:M
%         Corre2(i,j)=x(2*M+kk);
%         kk=kk+1;
%     end
% end
% kk=1;
% for i=1:M
%     for j=i+1:M
%         Corre2(i,j)=x(2*M+kk);
%         kk=kk+1;
%     end
% end

l=1;kTransition=zeros(M,M);lambda=zeros(1,M);
for i=1:M
    if OtherPara.lambdaMatrix(i)~=0
        lambda(i)=ParaToFit(l);%
        l=l+1;
    end
end
for i=1:M
    if OtherPara.muMatrix(i)~=0
        mu(i)=ParaToFit(l);%
        l=l+1;
    end
end
for i=1:M
    for j=1:M
        if OtherPara.Connectivity(i,j)~=0
        kTransition(i,j)=ParaToFit(l);%
        l=l+1;
        end
    end
    
end

NFit=ParaToFit(end);


for i=1:size(kTransition,2)   
K(i)=lambda(i)-mu(i)-(sum(kTransition(i,:),2)-kTransition(i,i));
KPlus(i)=lambda(i)+mu(i)+(sum(kTransition(i,:),2)-kTransition(i,i));
end

%disp(kTransition);disp(K(2));disp(lambda(2)-mu(2));disp(K(1));

% %ODEs Carrying capcity 1:
% for i=1:M
% ODE(i)=K(i)*x(i)+(kTransition(:,i)'*x(1:M)-kTransition(i,i)*x(i));
% end
% i=1;
% ODE(1)=(lambda(i)-mu(i))*x(i)*(1-x(i)/NFit)+(-(sum(kTransition(i,:),2)-kTransition(i,i)))*x(i)+(kTransition(:,i)'*x(1:M)-kTransition(i,i)*x(i));


%ODEs generation capacity new 2:
for i=1
    G=log2(x(i)./OtherPara.InitialValue(1))*lambda(i)/(lambda(i)-kTransition(i,2));
ODE(i)=(lambda(i))*x(i)*(1-G/NFit)-mu(i)*x(i)+(-(sum(kTransition(i,:),2)-kTransition(i,i)))*x(i)+(kTransition(:,i)'*x(1:M)-kTransition(i,i)*x(i));
end
for i=2:M
%ODE(i)=(lambda(i))*x(i)*(1-sum(x)/NFit)-mu(i)*x(i)+(-(sum(kTransition(i,:),2)-kTransition(i,i)))*x(i)+(kTransition(:,i)'*x(1:M)-kTransition(i,i)*x(i));
ODE(i)=(lambda(i))*x(i)-mu(i)*x(i)+(-(sum(kTransition(i,:),2)-kTransition(i,i)))*x(i)+(kTransition(:,i)'*x(1:M)-kTransition(i,i)*x(i));
end



% if length(ParaToFit)>=l
%     %disp('HSPC');
% ODE(1)=ODE(1)+ParaToFit(end);%HSPC with constant production;
% end

% for i=1:M
%     Matrix1=kTransition.*Corre2;%\sum k_mi*<n_m*n_i>
%     ODE(i+M)=KPlus(i)*x(i)+(kTransition(:,i)'*x(1:M)-kTransition(i,i)*x(i))...
%             +2*K(i)*x(i+M)+2*(sum(Matrix1(:,i),1)-Matrix1(i,i));
% end
% 
% kk=0;
% for i=1:M
%     for j=i+1:M
%     kk=kk+1;      
%     ODE(kk+2*M)=(K(i)+K(j))*x(kk+2*M)...
%         +(kTransition(:,i)'*Corre2(:,j)-kTransition(j,i)*Corre2(j,j))+kTransition(j,i)*x(j+M)-kTransition(j,i)*x(j)...
%         +(kTransition(:,j)'*Corre2(:,i)-kTransition(i,j)*Corre2(i,i))+kTransition(i,j)*x(i+M)-kTransition(i,j)*x(i);
%     end
% end
ODE=ODE';
%disp(ODE);
end

%%

solution=solution';
%solution=solution(:,2:end);%Delete t=0 time point if staring from 0
%display(size(solution));
cutoff=0.000001;


for ii=1:M
Mean(ii,:)=solution(ii,:);
%Var(ii,:)=solution(ii+M,:);
Var(ii,:)=Mean(ii,:);%solution(ii+M,:)-Mean(ii,:).^2;
end
CV=Mean;%sqrt(Var)./max(Mean,cutoff);
%display(CV);
kkk=0;
for ii=1:M
    for jj=ii+1:M
    kkk=kkk+1;      
    Corre(kkk,:)=Mean(ii,:);%solution(kkk+2*M,:)-Mean(ii,:).*Mean(jj,:))./sqrt(max(Var(ii,:),cutoff))./sqrt(max(Var(jj,:),cutoff));
    end
end
%display(Corre);
 sumation=max(sum(Mean,1),cutoff);%%

 for ii=1:M
 Ratio(ii,:)=Mean(ii,:)./sumation;%%
 end
 
%  if label==7 || label==8
%  for ii=1:M
%         Mean(ii,:)=Mean(ii,:)/MaxMean(ii);
%  end
%  end
%disp(Mean);
%disp(OtherPara.NumDataToUse);
%display(reshape(Mean(:,OtherPara.NumDataToUse)',1,[]));
FitValue=...[reshape(Ratio(:,OtherPara.NumDataToUse)',1,[]),...
    [reshape(Mean(:,OtherPara.NumDataToUse)',1,[])/OtherPara.NormalFactor];%,...
   % reshape(CV(:,OtherPara.NumDataToUse)',1,[])];%,reshape(Corre(:,OtherPara.NumDataToUse)',1,[])];%%
%FitValue=[Mean(2:4,end)',CV(2:4)',Corre(4:6)'];
%display(FitValue);%%
%display(ParaToFit);

FitValue=real(FitValue);
%FitValue=[FitValue,FitValue];
%disp(size(FitValue));
%display(ParaToFit);

end