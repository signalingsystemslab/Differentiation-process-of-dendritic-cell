 function CheckRates3(tt,DataTimePoint0,DataToFit0,Mean,CV,Corre,OtherPara,ParaBestFit,ParaFitSummary,CellTypeLabel,WTorKO)

%% Plot Rates 
Species=size(DataToFit0.RelativeMean,1);

ProliferationRate=ParaFitSummary(:,1:Species);
DeathRate=ParaFitSummary(:,Species+1:Species+Species);
Net_ProliferationRate=ProliferationRate-DeathRate;
DifferentiationRate=ParaFitSummary(:,2*Species+1:end-4);
chi2=ParaFitSummary(:,end-2);
[chiValues,Index] = sort(chi2);
chiValues=abs(log(chiValues));
DifferentiationRate=DifferentiationRate(Index,:);
Net_ProliferationRate=Net_ProliferationRate(Index,:);

figure('position', [00, 00, 1800, 800])
%c=jet(size(Net_ProliferationRate),1);
FigSize=round(sqrt(Species))+1;
for kk=1:Species
    subplot(FigSize,FigSize,kk)
b=bar(Net_ProliferationRate(:,kk)','FaceColor','flat');
%text([1:length(Net_ProliferationRate(:,kk))], Net_ProliferationRate(:,kk)', num2str(round(Net_ProliferationRate(:,kk),1)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16)
for ii=1:size(chiValues)
    b.CData(ii,:) = [0.5 0.5 chiValues(ii)/max(chiValues)];
end

legend(num2str(mean(Net_ProliferationRate(:,kk))),'Location','best');
set(gca,'FontSize',16,'linewidth',2);
if kk==1
    %ylim([0 2]);
elseif kk==2
    %ylim([0 5]);
elseif kk==3
    %ylim([-3 0]);
elseif kk==4
    %ylim([-1 1]);
elseif kk==5
    %ylim([-2 1]);
end
end
figurename=[OtherPara.folder,'\Net_PreliRatesDistri_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
print(gcf, '-djpeg', '-r300',figurename);%%print(gcf, '-dsvg',figurename);
figurename=[OtherPara.folder,'\Net_PreliRatesDistri_',num2str(WTorKO),'_',num2str(OtherPara.label),'.svg'];
%print(gcf, '-dsvg', '-r300',figurename);


%subplot(1,2,2)
figure('position', [00, 00, 1800, 800])
FigSize=round(sqrt(size(DifferentiationRate,2)))+1;
for kk=1:size(DifferentiationRate,2)
    
    subplot(FigSize,FigSize,kk)
b=bar(DifferentiationRate(:,kk)','FaceColor','flat');
%text([1:length(Net_ProliferationRate(:,kk))], Net_ProliferationRate(:,kk)', num2str(round(Net_ProliferationRate(:,kk),1)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16)
for ii=1:size(chiValues)
    b.CData(ii,:) = [0.5 chiValues(ii)/max(chiValues) 0.5];
end
%b.FaceColor = 'r';
legend(num2str(mean(DifferentiationRate(:,kk))),'Location','best');
set(gca,'FontSize',16,'linewidth',2);
if kk==1
    %ylim([0 0.03]);
elseif kk==2
    %ylim([0 4]);
elseif kk==3
    %ylim([0 3]);
elseif kk==4
    %ylim([0 1]);
end
%text([1:length(DifferentiationRate(:,kk))], DifferentiationRate(:,kk)', num2str(round(DifferentiationRate(:,kk),1)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16)
end
figurename=[OtherPara.folder,'\DiffRatesDistri_',num2str(WTorKO),'_',num2str(OtherPara.label),'.jpg'];
saveas(gcf,figurename);
figurename=[OtherPara.folder,'\DiffRatesDistri_',num2str(WTorKO),'_',num2str(OtherPara.label),'.svg'];
%print(gcf, '-dsvg', '-r300',figurename);
 end
