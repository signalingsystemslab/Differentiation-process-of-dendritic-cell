%(c) 2021 Signaling Systems Lab UCLA
%All rights reserved. 
%This MATLAB code implements the mathematical modeling on the differentiation 
%process of dendritic cell population. It fitted the birth-death ODE to 
%the cell number measurements and inferred the differential rates.
%Contact: jamestang23@gmail.com

%A detailed application of this package is given in the main text. 
%%
clear;
filefolder0=pwd;
Scan=10;%100;% The number of numerical runs to search the optimal parameters
CarryingCapcity=1;% Use generation capacity  for HSPC

FixProliferationbyNew=1;%Fix net-proliferation by estimation from experiment
FixProliferationWithFixRatio=0;%Fix ratio for all the net-proliferation
PlotSum=0;%1 is to plot summariezed rates.
ProliferationMode=1;

if FixProliferationbyNew==1
     %Now here
    filefolder0=[filefolder0,'\Final_WTKO'];
end
%Fix ratio for all the net-proliferation: effective 1 paramter
if FixProliferationWithFixRatio==1
 filefolder0=[filefolder0,'FixRatio'];
end
 
%Carrying capacity: manually specify the name of the output file, need to
%change the way of parameter constrain for optimiation in Main5.m: search 
%"Generation capacity" in Main5.m and comment the current optimization
%constrain and uncomment the other types of constrain.
if CarryingCapcity==1
  %filefolder0=[filefolder0,'\NewGenerationCapacity_',num2str(CarryingCapcity)];%'\New_1000_1'];%'\20190806_Data2_WTHET_1_1000'];%,'\20190806_Data2_WTKO_1_1000'];%'\20190806_Data2_WTHET_1_1000'];%'\20190806_Data2_WTKO_1_const_1000'];%'\20190801_WTHET_1_1000'];%'\20190801_WTKO_1_1000']
  %filefolder0=[filefolder0,'\NewGenerationCapacityFixProliferation_',num2str(CarryingCapcity)];%'\New_1000_1'];%'\20190806_Data2_WTHET_1_1000'];%,'\20190806_Data2_WTKO_1_1000'];%'\20190806_Data2_WTHET_1_1000'];%'\20190806_Data2_WTKO_1_const_1000'];%'\20190801_WTHET_1_1000'];%'\20190801_WTKO_1_1000']
  filefolder0=[filefolder0,'\NewGenerationCapacityFixDifferential_',num2str(CarryingCapcity)];%'\New_1000_1'];%'\20190806_Data2_WTHET_1_1000'];%,'\20190806_Data2_WTKO_1_1000'];%'\20190806_Data2_WTHET_1_1000'];%'\20190806_Data2_WTKO_1_const_1000'];%'\20190801_WTHET_1_1000'];%'\20190801_WTKO_1_1000']
  %filefolder0=[filefolder0,'\NewGenerationCapacityFixGenerationCapacity_',num2str(CarryingCapcity)];%'\New_1000_1'];%'\20190806_Data2_WTHET_1_1000'];%,'\20190806_Data2_WTKO_1_1000'];%'\20190806_Data2_WTHET_1_1000'];%'\20190806_Data2_WTKO_1_const_1000'];%'\20190801_WTHET_1_1000'];%'\20190801_WTKO_1_1000']
  % filefolder0=[filefolder0,'\NewGenerationCapacityFixProliferAndGenerationCapacity_',num2str(CarryingCapcity)];%'\New_1000_1'];%'\20190806_Data2_WTHET_1_1000'];%,'\20190806_Data2_WTKO_1_1000'];%'\20190806_Data2_WTHET_1_1000'];%'\20190806_Data2_WTKO_1_const_1000'];%'\20190801_WTHET_1_1000'];%'\20190801_WTKO_1_1000']
end

if ~exist(filefolder0)
mkdir(filefolder0);
end
if CarryingCapcity~=0
LabelChoice=[2];
end
kk=1;
for label=LabelChoice%1:5%5:8%4:4%1:3
%WTorKO=2;
CellCount_WT=[];CellCount_KO=[];
%CellTypeLabel={'HSPC','preDC','prepDC','pDC','cDC1','cDC2'};% need merge=1
CellTypeLabel={'HSPC','preDC','pDC','cDC1','cDC2'};% need merge=0
CellType=CellTypeLabel;

% % WT&KO
  [CellCount_WT_raw, txt]= xlsread('YL524-data.xlsx',1,'B12:P18');
 [TimePoints_WT, txt]= xlsread('YL524-data.xlsx',1,'A12:A18');
 [CellCount_KO_raw, txt]= xlsread('YL524-data.xlsx',1,'B32:P38');
  [TimePoints_KO, txt]= xlsread('YL524-data.xlsx',1,'A32:A38');
  

SampleSize=3;  
%CellCount_WT=zeros(size(CellType,2),size(TimePoints_WT,1),size(CellCount_WT_raw,1));
for k=1:SampleSize %replicate: didn't use the last one, because there is NaN
    for i=1:size(CellType,2)
        CellCount_WT(i,:,k)=CellCount_WT_raw(:,SampleSize*(i-1)+k:SampleSize*(i-1)+k);
        CellCount_KO(i,:,k)=CellCount_KO_raw(:,SampleSize*(i-1)+k:SampleSize*(i-1)+k);
    end
end
%% input proliferation rates 
 if FixProliferationbyNew==1
 %Fix net-proliferation by estimation from experiment
  Ki_67WT(1,1:5)= [4, 2, 1, 1,1];%scRNAseqAssump;
  Ki_67WT(2,1:5)=[4, 2, 1, 1,1];
  Ki_67KO(1,1:5)= [5,2,1.3,1,1];%scRNAseqAssump;
  Ki_67KO(2,1:5)=[5,2,1.3,1,1];
  Normalization=max(max(Ki_67WT));
  Ki_67WT=Ki_67WT/Normalization;%Normalize by WT maximum 
  Ki_67KO=Ki_67KO/Normalization;%Normalize by WT
 end
 Ki_67WT(1,:)=[];Ki_67KO(1,:)=[]; 
 Ki_67WTUse=Ki_67WT;Ki_67KOUse=Ki_67KO;

%% Specify the time-wise regime
consistency=0;

TimeLeft=1;
if label==1
    TimeRight=7;  
elseif label==2
     TimeLeft=1;
     TimeRight=7;
%region=1;
end
region=label;
 if CarryingCapcity~=0
      TimeLeft=1;
    TimeRight=7;
 end

DataRegime=TimeLeft:TimeRight;

if consistency==1
    DataRegime=[1,4,6];
end

CellCount_WT=CellCount_WT(:,DataRegime,:);
CellCount_KO=CellCount_KO(:,DataRegime,:);
TimePoints_WT=TimePoints_WT(DataRegime);
TimePoints_KO=TimePoints_KO(DataRegime);

filefolder=[filefolder0,'\Region_',num2str(region)];
if ~exist(filefolder)
    mkdir(filefolder);
end


nCT=size(CellCount_WT,1);
CellType=fliplr(CellType);
TimePoints_WT=TimePoints_WT';TimePoints_KO=TimePoints_KO';

for j=1:size(TimePoints_WT,2)
    Correlation_WT(:,:,j)=corrcoef(reshape(CellCount_WT(:,j,:),size(CellType,2),SampleSize)'); %correlation for each time point
    Correlation_KO(:,:,j)=corrcoef(reshape(CellCount_KO(:,j,:),size(CellType,2),SampleSize)');
end

%% Collect data for simulation
for WTorKO=1:2
if WTorKO==1
    TimePoints=TimePoints_WT;MeanCountUse=mean(CellCount_WT,3);
    VarToUse=var(CellCount_WT,1,3);StdToUse=std(CellCount_WT,1,3);
    CorreToUse=Correlation_WT;
   DataToFit0.lambdaFitRatio=Ki_67WTUse;     
else
    %nCT=size(CellCount_KO,1);
    TimePoints=TimePoints_KO;MeanCountUse=mean(CellCount_KO,3);
    VarToUse=var(CellCount_KO,1,3); StdToUse=std(CellCount_KO,1,3);
    CorreToUse=Correlation_KO;
    DataToFit0.lambdaFitRatio=Ki_67KOUse;
end
DataToFit0.RelativeMean=MeanCountUse./sum(MeanCountUse,1);%[];%,...   
DataToFit0.AbsoluteMean=MeanCountUse;
DataToFit0.CV=StdToUse./MeanCountUse;%[];
DataToFit0.Corre=[];
for j=1:size(TimePoints_WT,2)
    k=1;
for i=1:size(CellType,2)
    for ii=i+1:size(CellType,2)
    DataToFit0.Corre(k,j)=CorreToUse(i,ii,j);
    if j==1
        Initial.Corre(k)=DataToFit0.Corre(k,1)*StdToUse(i,1)*StdToUse(ii,1)+MeanCountUse(i,1).*MeanCountUse(ii,1);%zeros(nCT*(nCT-1)/2,1)';
    end
    k=k+1;
    end
end
end

Initial.AbsoluteMean=DataToFit0.AbsoluteMean(:,1)';
Initial.SecondMoment=(StdToUse(:,1).*StdToUse(:,1)+MeanCountUse(:,1).*MeanCountUse(:,1))';%zeros(nCT,1)';
DataToFit0.NormalFactor=max(max(DataToFit0.AbsoluteMean));
DataToFit=...[reshape(DataToFit0.RelativeMean',1,[]),...
    [reshape(DataToFit0.AbsoluteMean',1,[])/DataToFit0.NormalFactor];%,reshape(DataToFit0.CV',1,[])];%,reshape(DataToFit0.Corre',1,[])];

%% Fit the data to the model to infer parameters
ParaFitSummary=Main5(filefolder,nCT,TimePoints,Initial,DataToFit0,DataToFit,Scan,label,SampleSize,CellTypeLabel,WTorKO,PlotSum,ProliferationMode);
close all;

ParaFitSummaryMaster{kk}=ParaFitSummary;
kk=kk+1;
end
end


%% Plot the violin distribution of fitting accuracy 
figure('position', [00, 00, 800, 600])
title('Fitting accuracy');
colors = jet(length(ParaFitSummaryMaster));
for ii=1:length(ParaFitSummaryMaster)
    scaleViolin=mean(ParaFitSummaryMaster{ii}(:,end-2));
    ViolinChi=[ParaFitSummaryMaster{ii}(:,end-2)];
violinPlot2(ViolinChi*100,'x',[1*ii .7 1 1.8],'facecolor',colors(ii, : )...[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1]
,'edgecolor','none','bw',1,'mc',[],'medc',[])

ChiSquareSummary(:,ii)=ViolinChi*100;%ParaFitSummaryMaster{ii}(:,end-2)';
end
alpha(1);
 if CarryingCapcity~=0
xlabelstring={'WT','KO'};
 else
 xlabelstring={'WT 1','KO 1','WT 2','KO 2'};
 end

set(gca, 'xlim', [0 length(ParaFitSummaryMaster)+1],'XTick',1:length(ParaFitSummaryMaster),'XTickLabel',string(xlabelstring));
set(gca,'TickLabelInterpreter','none');
ylabel('Mean relative error (%)');
xtickangle(45)
ylim([0 200]);
hold off
set(gca,'FontSize',24,'linewidth',2);
figurename=[filefolder0,'\FittingGoodnessAll.jpg'];
print(gcf, '-djpeg', '-r300',figurename); 
figurename=[filefolder0,'\FittingGoodnessAll.svg'];
%print(gcf, '-dsvg', '-r300',figurename);
close all;
ExcelName=[filefolder0,'\ChiSquareSummaryAll.csv'];
csvwrite(ExcelName,ChiSquareSummary);

