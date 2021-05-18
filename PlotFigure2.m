function PlotFigure2(tt,DataTimePoint0,DataToFit0,Mean,CV,Corre,OtherPara,ParaBestFit,CellTypeLabel,WTorKO)

c=jet(size(Corre,1));
c1=hsv(size(Mean,1));
CorreLegend=cell(size(Corre,1),1);
kk=1;
for jj=1:size(Mean,1)
    for jjj=jj+1:size(Mean,1)
        CorreLegend{kk,1}=[num2str(jj),'-',num2str(jjj)];
        kk=kk+1;
    end
end

%% Plot effective proliferation rates from carrying capacity

figure('position', [00, 00, 900, 600])
%for ii=1:size(Mean,1)
G=log2(sum(Mean(1,:),1)./OtherPara.InitialValue(1))*ParaBestFit(WTorKO,1)/(ParaBestFit(WTorKO,1)-ParaBestFit(WTorKO,11));
ProliferationPropensity=(1-G/ParaBestFit(WTorKO,15));
assignin('base', 'ParaBestFit', ParaBestFit); 
plot(tt,ProliferationPropensity,'linewidth',2,'color','k');hold on;
%text(max(tt)*0.2,max(ylim)*0.95,['chi^2=',num2str(round(ParaBestFit(WTorKO,end-3)*100,1)),'%'],'FontSize',24)
ylim([0 1]);%set(gca,'yscale','log');
legend(['Capacity ',num2str(round(ParaBestFit(WTorKO,15),2))],'Location','best');
ylabel('Proliferation propensity');xlabel('Time (d)');xlim([min(tt) max(tt)]);
set(gca,'FontSize',24,'linewidth',2);
figurename=[OtherPara.folder,'\EffectiveProli_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
print(gcf, '-djpeg', '-r300',figurename);%%print(gcf, '-dsvg',figurename);
figurename=[OtherPara.folder,'\EffectiveProli_',num2str(WTorKO),'_',num2str(OtherPara.label),'.svg'];
%print(gcf, '-dsvg', '-r300',figurename);
filename=[OtherPara.folder,'\EffectiveProli_',num2str(WTorKO),'.mat'];
save(filename,'ProliferationPropensity');



%%
figure('position', [00, 00, 900, 600])
for ii=1:size(Mean,1)
plot(tt,Mean(ii,:),'linewidth',2,'color',c1(ii,:));hold on;
end
for ii=1:size(Mean,1)
%plot(DataTimePoint0.AbsoluteMean(ii,:),DataToFit0.AbsoluteMean(ii,:),'.', 'markersize',40,'color',c1(ii,:));
errorbar(DataTimePoint0.AbsoluteMean(ii,:),DataToFit0.AbsoluteMean(ii,:),DataToFit0.CV(ii,:).*DataToFit0.AbsoluteMean(ii,:),'.', 'markersize',40,'color',c1(ii,:));
end
legend(CellTypeLabel,'Location','bestoutside');
%text(max(tt)*0.2,max(ylim)*0.95,['chi^2=',num2str(round(ParaBestFit(WTorKO,end-3)*100,1)),'%'],'FontSize',24)
ylim([1e2 2e6]);set(gca,'yscale','log');
ylabel('Absolute cell number');xlabel('Time (d)');xlim([min(tt) max(tt)]);
set(gca,'FontSize',24,'linewidth',2);
figurename=[OtherPara.folder,'\CellNumberLog_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
print(gcf, '-djpeg', '-r300',figurename);%%print(gcf, '-dsvg',figurename);
figurename=[OtherPara.folder,'\CellNumberLog_',num2str(WTorKO),'_',num2str(OtherPara.label),'.svg'];
%print(gcf, '-dsvg', '-r300',figurename);

%return;

figure('position', [00, 00, 900, 600])
for ii=1:size(Mean,1)
plot(tt,Mean(ii,:),'linewidth',2,'color',c1(ii,:));hold on;
end
for ii=1:size(Mean,1)
%plot(DataTimePoint0.AbsoluteMean(ii,:),DataToFit0.AbsoluteMean(ii,:),'.', 'markersize',40,'color',c1(ii,:));
errorbar(DataTimePoint0.AbsoluteMean(ii,:),DataToFit0.AbsoluteMean(ii,:),DataToFit0.CV(ii,:).*DataToFit0.AbsoluteMean(ii,:),'.', 'markersize',40,'color',c1(ii,:));
end
legend(CellTypeLabel,'Location','bestoutside');
ylabel('Absolute cell number');xlabel('Time (d)');xlim([min(tt) max(tt)]);
set(gca,'FontSize',24,'linewidth',2);
figurename=[OtherPara.folder,'\CellNumber_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
print(gcf, '-djpeg', '-r300',figurename);%%print(gcf, '-dsvg',figurename);
figurename=[OtherPara.folder,'\CellNumber_',num2str(WTorKO),'_',num2str(OtherPara.label),'.svg'];
%print(gcf, '-dsvg', '-r300',figurename);

CellNumber=[tt',Mean'];
ExcelName=[OtherPara.folder,'\CellNumber_',num2str(WTorKO),'_',num2str(OtherPara.label),'.csv'];
csvwrite(ExcelName,CellNumber);




figure('position', [00, 00, 900, 600])
for ii=1:size(Mean,1)
plot(tt,CV(ii,:),'linewidth',2,'color',c1(ii,:));hold on;%plot(tt,CV(3,:),'linewidth',2);plot(tt,CV(4,:),'linewidth',2);
end
for ii=1:size(Mean,1)
if DataToFit0.CV(1,1)~=0
plot(DataTimePoint0.CV(ii,:),DataToFit0.CV(ii,:),'.', 'markersize',40,'color',c1(ii,:));
end
end
legend(CellTypeLabel,'Location','bestoutside');
ylabel('CV');xlabel('Time (d)');xlim([min(tt) max(tt)]);%ylim([0 6]);
set(gca,'FontSize',24,'linewidth',2);
% figurename=[OtherPara.folder,'\CV',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
% print(gcf, '-djpeg', '-r300',figurename);%%print(gcf, '-dsvg',figurename);

figure('position', [00, 00, 900, 600])
for ii=1:size(Corre,1)
plot(tt,Corre(ii,:),'linewidth',2,'color',c(ii,:));hold on;%plot(tt,Corre(5,:),'linewidth',2);plot(tt,Corre(6,:),'linewidth',2);
end
for ii=1:size(Corre,1)
if size(DataToFit0.Corre,1)>1
plot(DataTimePoint0.Corre(ii,:),DataToFit0.Corre(ii,:),'.', 'markersize',40,'color',c(ii,:));
end
end
legend(CorreLegend,'Location','bestoutside');
ylabel('Correlation coefficient');xlabel('Time (d)');xlim([min(tt) max(tt)]);%ylim([-0.2 10]);
set(gca,'FontSize',24,'linewidth',2);
% figurename=[OtherPara.folder,'\Correlation',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
% print(gcf, '-djpeg', '-r300',figurename);%%print(gcf, '-dsvg',figurename);


%% ratio plot 
sumation=sum(Mean,1);%%
 
figure('position', [00, 00, 900, 600])
 clf;
for ii=1:size(Mean,1) 
plot(tt,Mean(ii,:)./sumation,'linewidth',2,'color',c1(ii,:));hold on;
end
for ii=1:size(Mean,1)
plot(DataTimePoint0.RelativeMean(ii,:),DataToFit0.RelativeMean(ii,:),'.', 'markersize',40,'color',c1(ii,:));
end
legend(CellTypeLabel,'Location','bestoutside');
ylabel('Ratio');xlabel('Time (d)');xlim([min(tt) max(tt)]);
set(gca,'FontSize',24,'linewidth',2);
%% Plot Rates 
Species=size(Mean,1);

ProliferationRate=ParaBestFit(WTorKO,1:Species);
DeathRate=ParaBestFit(WTorKO,Species+1:Species+Species);
Net_ProliferationRate=ProliferationRate-DeathRate;
DifferentiationRate=ParaBestFit(WTorKO,2*Species+1:end-4);

%% Plot Rates WT/KO
ProliferationRateRatio=ParaBestFit(2,1:Species)./ParaBestFit(1,1:Species);
DeathRateRatio=ParaBestFit(2,Species+1:Species+Species)./ParaBestFit(1,Species+1:Species+Species);
Net_ProliferationRateRatio=(ParaBestFit(2,1:Species)-ParaBestFit(2,Species+1:Species+Species))./(ParaBestFit(1,1:Species)-ParaBestFit(1,Species+1:Species+Species));
DifferentiationRateRatio=ParaBestFit(2,2*Species+1:end-4)./ParaBestFit(1,2*Species+1:end-4);

return;
end