 function PlotSumRates3(tt,DataTimePoint0,DataToFit0,Mean,CV,Corre,OtherPara,ParaBestFit,ParaFitSummary,CellTypeLabel,WTorKO)

%% Plot Rates 
Species=size(DataToFit0.RelativeMean,1);

ProliferationRate=ParaFitSummary(:,1:Species);
DeathRate=ParaFitSummary(:,Species+1:Species+Species);
Net_ProliferationRate=ProliferationRate-DeathRate;
DifferentiationRate=ParaFitSummary(:,2*Species+1:end-5);
GenerationCapacity=ParaFitSummary(:,end-4);
chi2=ParaFitSummary(:,end-2);
[chiValues,Index] = sort(chi2);
chiValues=abs(log(chiValues));

DifferentiationRate=DifferentiationRate(Index,:);
Net_ProliferationRate=Net_ProliferationRate(Index,:);
GenerationCapacity=GenerationCapacity(Index,:);

%%
figure('position', [00, 00, 800, 600])
Cutoff=min(50,round(length(chi2)*0.5)); %top-100 ranked fitting
for kk=1:Species
    
Net_ProliferationRateMean(kk)=mean(Net_ProliferationRate(1:Cutoff,kk),1);
Net_ProliferationRateStd(kk)=std(Net_ProliferationRate(1:Cutoff,kk),1);
end
color=hsv(length(Net_ProliferationRateMean));
c = categorical(CellTypeLabel,CellTypeLabel); 
%b=bar(categorical(CellTypeLabel),Net_ProliferationRateMean);
b =bar(c,Net_ProliferationRateMean);%,'FaceColor','flat');
for k = 1:length(Net_ProliferationRateMean)
    b.FaceColor = 'flat';
    b.CData(k,:) = color(k,:);%[.5 0 .5];
    %b(k).CData = color(k,:);
end

hold on

er = errorbar(1:Species,Net_ProliferationRateMean,Net_ProliferationRateStd,'linewidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
%er.linewidth = 2;
title('Net-preliferation rates');
ylim([0 1.5]);
%%ylim([0 1.2]);
hold off
set(gca,'FontSize',22,'linewidth',2);

figurename=[OtherPara.folder,'\Net_PreliRatesSummary_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
print(gcf, '-djpeg', '-r300',figurename);%%print(gcf, '-dsvg',figurename);
figurename=[OtherPara.folder,'\Net_PreliRatesSummary_',num2str(WTorKO),'_',num2str(OtherPara.label),'.svg'];
%print(gcf, '-dsvg', '-r300',figurename);
Net_ProliferationSummary=[Net_ProliferationRateMean',Net_ProliferationRateStd'];
ExcelName=[OtherPara.folder,'\Net_PreliRatesSummary_',num2str(WTorKO),'_',num2str(OtherPara.label),'.csv'];
csvwrite(ExcelName,Net_ProliferationSummary);

GenerationCapacityCutoff=GenerationCapacity(1:Cutoff);
ExcelName=[OtherPara.folder,'\GenerationCapacity_',num2str(WTorKO),'_',num2str(OtherPara.label),'.csv'];
csvwrite(ExcelName,GenerationCapacityCutoff);
%%
figure('position', [00, 00, 800, 600])
for kk=1:Species
Temp11(kk)=mean(ProliferationRate(1:Cutoff,kk),1);
Net_ProliferationRateStd(kk)=std(ProliferationRate(1:Cutoff,kk),1);
end
color=hsv(length(Temp11));
c = categorical(CellTypeLabel,CellTypeLabel); 
b =bar(c,Temp11);%,'FaceColor','flat');
for k = 1:length(Temp11)
    b.FaceColor = 'flat';
    b.CData(k,:) = color(k,:);%[.5 0 .5];
    %b(k).CData = color(k,:);
end
hold on
er = errorbar(1:Species,Temp11,Net_ProliferationRateStd,'linewidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
%er.linewidth = 2;
title('Preliferation rates');
hold off
set(gca,'FontSize',22,'linewidth',2);
figurename=[OtherPara.folder,'\PreliRatesSummary_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
%%%saveas(gcf,figurename); 

%%
figure('position', [00, 00, 800, 600])
%c=jet(size(Net_ProliferationRate),1);
%Cutoff=10; %top-100 ranked fitting
for kk=1:Species
Temp11(kk)=mean(DeathRate(1:Cutoff,kk),1);
Net_ProliferationRateStd(kk)=std(DeathRate(1:Cutoff,kk),1);
end
color=hsv(length(Temp11));
c = categorical(CellTypeLabel,CellTypeLabel); 
b =bar(c,Temp11);%,'FaceColor','flat');
for k = 1:length(Temp11)
    b.FaceColor = 'flat';
    b.CData(k,:) = color(k,:);%[.5 0 .5];
    %b(k).CData = color(k,:);
end
hold on
er = errorbar(1:Species,Temp11,Net_ProliferationRateStd,'linewidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
%er.linewidth = 2;
title('Death rates');
hold off
set(gca,'FontSize',22,'linewidth',2);
figurename=[OtherPara.folder,'\DeathRatesSummary_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
%%%saveas(gcf,figurename); 


%%
%subplot(1,2,2)
figure('position', [00, 00, 800, 700])
%c=jet(size(Net_ProliferationRate),1);
for kk=1:size(DifferentiationRate,2)
Net_DifferentiationRateMean(kk)=mean(DifferentiationRate(1:Cutoff,kk),1);
Net_DifferentiationRateStd(kk)=std(DifferentiationRate(1:Cutoff,kk),1);
end
CellTypeLabel2={'HSPC->preDC','preDC->pDC','preDC->cDC1','preDC->cDC2'};
c = categorical(CellTypeLabel2,CellTypeLabel2); 
%b=bar(categorical(CellTypeLabel),Net_ProliferationRateMean);
b=bar(c,Net_DifferentiationRateMean);
color2=parula(length(Net_DifferentiationRateMean));
for k = 1:length(Net_DifferentiationRateMean)
    b.FaceColor = 'flat';
    b.CData(k,:) = color2(k,:);%[.5 0 .5];
    %b(k).CData = color(k,:);
end
hold on

er = errorbar(1:size(DifferentiationRate,2),Net_DifferentiationRateMean,Net_DifferentiationRateStd,'linewidth',2);    
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
%er.linewidth = 2;
title('Differentiation rates');
%ylim([0 1.2]);
ylim([0 0.6]);
hold off
set(gca,'FontSize',22,'linewidth',2);


figurename=[OtherPara.folder,'\DiffRatesSummary_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
print(gcf, '-djpeg', '-r300',figurename);
figurename=[OtherPara.folder,'\DiffRatesSummary_',num2str(WTorKO),'_',num2str(OtherPara.label),'.svg'];
%print(gcf, '-dsvg', '-r300',figurename);
DifferentialRateSummary=[Net_DifferentiationRateMean',Net_DifferentiationRateStd'];
ExcelName=[OtherPara.folder,'\DifferentialRates_',num2str(WTorKO),'_',num2str(OtherPara.label),'.csv'];
csvwrite(ExcelName,DifferentialRateSummary);

if WTorKO==1
ParaFitMean=[Net_ProliferationRateMean,Net_DifferentiationRateMean];
assignin('base', 'ParaFitMean',ParaFitMean);
end

%%
%plot input proliferation parameters:
figure('position', [00, 00, 800, 600])
c = categorical(CellTypeLabel,CellTypeLabel); 
b=bar(c,DataToFit0.lambdaFitRatio)
for k = 1:length(Net_ProliferationRateMean)
    b.FaceColor = 'flat';
    b.CData(k,:) = color(k,:);%[.5 0 .5];
    %b(k).CData = color(k,:);
end
hold on
title('Input net-preliferation rates');
%%%ylim([-1 1]);

%%ylim([0 2]);
hold off
set(gca,'FontSize',22,'linewidth',2);

figurename=[OtherPara.folder,'\Input_Net_PreliRatesSummary_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
%%%saveas(gcf,figurename); 


 end
