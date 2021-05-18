%(c) 2021 Signaling Systems Lab UCLA
%All rights reserved. 
%This MATLAB code implements the mathematical modeling on the differentiation 
%process of dendritic cell population. It fitted the birth-death ODE to 
%the cell number measurements and inferred the differential rates.
%Contact: jamestang23@gmail.com

%A detailed application of this package is given in the main text. 
%%
function ParaFitSummary=Main5(filefolder,nCT,TimePoints,Initial,DataToFit0,DataToFit,Scan,label,SampleSize,CellTypeLabel,WTorKO,PlotSum,ProliferationMode)
tic;
%%
%Input
global ParaBestFit OtherPara ParaToFit tt M 
syms t
tt=linspace(TimePoints(1),TimePoints(end),10000);%0.01:0.01:8;%[[0.001:0.001:0.1],[0.1:0.01:8]];
filename=[filefolder,'\Topology',num2str(WTorKO),'_',num2str(label),'.mat'];
OtherPara.NumDataToUse=1:1:size(TimePoints,2);%1:1:size(TimePoints,2);%[];%1;%size(TimePoints,2);%1:1:size(TimePoints,2);% for absolute statistics etc
OtherPara.NormalFactor=DataToFit0.NormalFactor;
MaxMean=OtherPara.NormalFactor;

%Number of cell states;
M=nCT;
ParaRange=2;
%Structure of the system 
%Topological connnectivity matrix 
Connectivity=[0,1,0,0,0;
    0,0,1,1,1;
    0,0,0,0,0;
    0,0,0,0,0;
    0,0,0,0,0];
lambdaMatrix=ones(1,M);
muMatrix=ones(1,M);
kTemp=reshape(Connectivity',[1,M*M]);
kTemp(kTemp==0)=[];
lambdaTemp=lambdaMatrix;lambdaTemp(lambdaTemp==0)=[];
muTemp=muMatrix;muTemp(muTemp==0)=[];
NumberPara=size([lambdaTemp,muTemp,kTemp],2);

%Initial value: may serve as a cell source for the differentiated states

InitialValue=...[Initial.RelativeMean,
    [Initial.AbsoluteMean,Initial.SecondMoment,Initial.Corre];
OtherPara.label=label;
OtherPara.folder=filefolder;
OtherPara.lambdaMatrix=lambdaMatrix;
OtherPara.muMatrix=muMatrix;
OtherPara.M=M;
OtherPara.InitialValue=InitialValue;
OtherPara.Connectivity=Connectivity;

%% The data to be fitted:
%Mean: time seris, abosolut number, CV, Correlation

%bootstrapping to get uncertainty
errorbar0.RelativeMean=[0.02,0.02,0.02,0.02,0.02,0.02,...
    0.05,0.05,0.05,0.05,0.05,0.05,...
    0.02,0.02,0.02,0.02,0.02,0.02];%...
errorbar0.AbsoluteMean=[10,10,10];%,...
errorbar0.CV=[0.1,0.4,0.5];%,...
errorbar0.Corre=[0.1,0.2,0.05];%Mannel input for now, should use bootstrapping, or directly from data.
errorbar=[errorbar0.RelativeMean,errorbar0.AbsoluteMean,errorbar0.CV,errorbar0.Corre];


OtherPara.TimePoint=TimePoints;%[1,2,3,4,6,8];%d
DataTimePoint0.RelativeMean=ones(M,1)*TimePoints;
DataTimePoint0.AbsoluteMean=ones(M,1)*TimePoints;
DataTimePoint0.CV=ones(M,1)*TimePoints;
DataTimePoint0.Corre=ones(M*(M-1)/2,1)*TimePoints;
if exist(filename)
    load(filename);  
else
%% fit the parameters: by numerical ODE
RsdnrmNow=Inf;
ParaFitSummary=[];
SuccessFit=0;SuccessScan=0;
for iter=1:Scan
    display(iter);

kTransition=(1+(rand(M,M)-0.5))*1;%rand(M,M)*ParaRange;%1+(rand(M,M)-0.5);
kTransition=kTransition.*Connectivity;

if label==3 || 4
   lambda=(DataToFit0.lambdaFitRatio+ones(1,M)).*(rand(1,M)-0.5); 
   mu=ones(1,M).*(rand(1,M)-0.5); 
else
    lambda=(1+(rand(1,M)-0.5))*1;%rand(1,M)*ParaRange;%1+(rand(1,M)-0.5);
    mu=(1+(rand(1,M)-0.5))*1;%rand(1,M)*ParaRange;%1+(rand(1,M)-0.5);
end

lambda=lambda.*lambdaMatrix;
mu=mu.*lambdaMatrix;


kFit=reshape(kTransition',[1,M*M]);kFit(kFit==0)=[];
lambdaFit=lambda;lambdaFit(lambdaFit==0)=[];
muFit=mu;muFit(muFit==0)=[];
%Carrying capcity
NFit=1e6*(rand);%based on WT


ProductHSPC=(1+(rand(1,1)-0.5))*1;
ParaToFit=[lambdaFit,muFit,kFit,NFit];%,ProductHSPC];%Add influx need modify ODE-Fitting2, RunODESimulation.m
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');

%%%%%%%%%%%%%%%%%%%%%%%%%Generation capacity begin%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%comment the current optimization constrain and uncomment the other types of constrain.
% %Final: strictly fix proliferation by Ki67; allow death [0.1,10] times change
disp('Constrain net-proliferation rates, add constrain on generation capacity to avoid numerical error');
lb = [DataToFit0.lambdaFitRatio*0.5,ones(1,length(muFit))*0,ones(1,length(kFit))*0.00001,0];
ub = [DataToFit0.lambdaFitRatio/0.5,ones(1,length(muFit))*0,ones(1,length(kFit))*100, 1e10];

%Fix KO parameters to WT optimal: [0.980424257789821,0.947991883028425,0.416834380504602,0.499999169241122,0.125006056905165,0,0,0,0,0,1.00000098391786e-05,0.00781985600711960,0.0597208403987025,0.469729998051095,1.83975303107402]
%Fix KO parameters to WT mean: [0.979154802632275,0.898823690578527,0.404904995231959,0.499876659606393,0.185845413456722,1.77138628042335e-05,0.0215128899748691,0.0595736537478091,0.407419780958218]

% disp('Fix net-proliferation rates to WT');
% if WTorKO==2
% lb = [0.979154802632275,0.898823690578527,0.404904995231959,0.499876659606393,0.185845413456722,ones(1,length(muFit))*0,ones(1,length(kFit))*0.00001,0];
% ub = [0.979154802632275,0.898823690578527,0.404904995231959,0.499876659606393,0.185845413456722,ones(1,length(muFit))*0,ones(1,length(kFit))*100, 1e10];
% end

% disp('Fix differentiation rates to WT');
% if WTorKO==2
% lb = [DataToFit0.lambdaFitRatio*0.5,ones(1,length(muFit))*0,1.77138628042335e-05,0.0215128899748691,0.0595736537478091,0.407419780958218,0];
% ub = [DataToFit0.lambdaFitRatio/0.5,ones(1,length(muFit))*0,1.77138628042335e-05,0.0215128899748691,0.0595736537478091,0.407419780958218, 1e10];
% end

% disp('Constrain Generation capacity to WT');
% if WTorKO==2
% lb = [DataToFit0.lambdaFitRatio*0.5,ones(1,length(muFit))*0,ones(1,length(kFit))*0.00001,1.83975303107402];
% ub = [DataToFit0.lambdaFitRatio/0.5,ones(1,length(muFit))*0,ones(1,length(kFit))*100, 1.83975303107402];
% end

% disp('Fix  net-proliferation and Constrain Generation capacity to WT');
% if WTorKO==2
% lb = [0.980424257789821,0.947991883028425,0.416834380504602,0.499999169241122,0.125006056905165,ones(1,length(muFit))*0,ones(1,length(kFit))*0.00001,1.83975303107402];
% ub = [0.980424257789821,0.947991883028425,0.416834380504602,0.499999169241122,0.125006056905165,ones(1,length(muFit))*0,ones(1,length(kFit))*100, 1.83975303107402];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%Generation capacity end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


objfcn = @(ParaToFit,TimePoint)ODE_Fitting5(real(ParaToFit),OtherPara.TimePoint,OtherPara,MaxMean,label);
[ParaFit,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat] = lsqcurvefit(objfcn,real(ParaToFit),OtherPara.TimePoint,DataToFit,lb,ub,options);

[chi2, chi2_old] = chi_squared(DataToFit,Rsd(1:length(DataToFit))+DataToFit,size(ParaToFit,2));%,errorbar);
[chi2HSPC, chi2HSPC_old] = chi_squared(DataToFit(1:7),Rsd(1:7)+DataToFit(1:7),size(ParaToFit,2));%HSPC
DataFitted=Rsd(1:length(DataToFit))+DataToFit;
AICc=2*length(ParaToFit)+2*chi2_old+2*length(ParaToFit)*(length(ParaToFit)+1)...
    /(SampleSize-length(ParaToFit)-1);
AICc_2=2*length(ParaToFit)+2*chi2+2*length(ParaToFit)*(length(ParaToFit)+1)...
    /(SampleSize-length(ParaToFit)-1);
ParaFitSummary=[ParaFitSummary',[ParaFit,chi2HSPC,chi2,AICc,AICc_2]']';

%evaluate likelihood by counting the number with \chi^{2}<1
if ExFlg~=0
    SuccessScan=SuccessScan+1;
end
if chi2_old<1
    SuccessFit=SuccessFit+1;
end
    
if Rsdnrm<RsdnrmNow
    ParaBestFit=[ParaFit,chi2HSPC,chi2,AICc,AICc_2,iter];
    RsdnrmNow=Rsdnrm;
end   
assignin('base', 'ParaFitSummaryTemp', ParaFitSummary);
end
toc;
save(filename,'OtherPara','DataTimePoint0','DataToFit0','DataFitted','DataToFit','ParaFitSummary','ParaBestFit','SuccessScan','SuccessFit','Scan');
end
%% Run ODE with best fit and plot
[Mean,CV,Corre]=RunODESimulation4(tt, M,ParaBestFit,OtherPara);

OtherPara.label=label;
OtherPara.folder=filefolder;

%Load the WT and corresponding KO for the ratio
if WTorKO==1 
    filename2=[filefolder,'\Topology',num2str(2),'_',num2str(label),'.mat'];
     ParaBestFit1=ParaBestFit;
    if exist(filename2)       
       load(filename2);   
    end
    ParaBestFitTemp=[ParaBestFit1;ParaBestFit];
elseif WTorKO==2 
    filename2=[filefolder,'\Topology',num2str(1),'_',num2str(label),'.mat'];
    ParaBestFit2=ParaBestFit;
   if exist(filename2)     
     load(filename2);
    end
    ParaBestFitTemp=[ParaBestFit;ParaBestFit2];
end
load(filename);

ParaBestFit=ParaBestFitTemp;

%% plot figures
PlotSumRates3(tt,DataTimePoint0,DataToFit0,Mean,CV,Corre,OtherPara,ParaBestFit,ParaFitSummary,CellTypeLabel,WTorKO);
PlotFigure2(tt,DataTimePoint0,DataToFit0,Mean,CV,Corre,OtherPara,ParaBestFit,CellTypeLabel,WTorKO);
CheckRates3(tt,DataTimePoint0,DataToFit0,Mean,CV,Corre,OtherPara,ParaBestFit,ParaFitSummary,CellTypeLabel,WTorKO);
return;


end