%% Power vs accuracy with some area lines
areas = load('simPointAreas.mat').simPointAreas;
powers = load('simPointPower.mat').simPointPower;
load("Detection results/SNR_Results.mat");
FOM_SAR = 8.8e-6/(98*512);
% Split correct into separate categories

interictalCorrect = zeros(100,1);
ind = 1:3:length(Correct);
ictalCorrect = zeros(100,1);
for l=1:length(ind)
    i=ind(l);
    interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
    ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
    correctPercent = (2*interictalCorrect+ictalCorrect)/3;
end

correctSortedByArea = nan(5, 100);
powersSortedbyArea = nan(5,100);
for i = 1:100
    if (areas(i)<5000)
        correctSortedByArea(ceil((areas(i)-1000)/1000), i) = correctPercent(i);
        powersSortedbyArea(ceil((areas(i)-1000)/1000),i) = powers(i);
    else
        correctSortedByArea(5, i) = correctPercent(i);
        powersSortedbyArea(5,i) = powers(i);
    end
    correctSortedByPower(ceil((powers(i)/(1e-6))),i) = correctPercent(i);
    areasSortedbyPower(ceil((powers(i)/(1e-6))),i) = areas(i);
    
end
% figure
% hold on

scatter(100*correctSortedByArea', 1e6*powersSortedbyArea', 50, "filled")


plot_paper

load("SAR_Results.mat");
SARPowers = load("SARPowers.mat").powers;
SARAreas = load("SARAreas.mat").areas;
SARcorrectPercent = zeros(1,12);
ind = 1:3:length(Correct);
for l=1:length(ind)
    i=ind(l);
    SARcorrectPercent(l+2) = (2*Correct(i)+2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i)*UncertainButCorrect(i)+2*Uncertain(i+1)*UncertainButCorrect(i+1)+Uncertain(i+2)*UncertainButCorrect(i+2))/5; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l+2) = (2*Uncertain(i)*UncertainButCorrect(i)+2*Uncertain(i+1)*UncertainButCorrect(i+1)+Uncertain(i+2)*UncertainButCorrect(i+2))/5;
end
hold on
plot(100*SARcorrectPercent(3:12), 1e6*SARPowers(3:12), "k--")
ylim([.5, 5])
xlim([65, 100])
legend(["1000-2000 C_{unit}", "2000-3000 C_{unit}","3000-4000 C_{unit}", "4000-5000 C_{unit}", ">5000 C_{unit}", "Baseline system"], 'location', 'northwest')
plot_paper
xlabel("Correct detection (%)")
ylabel("Power (\muW)")

% figure
% hold on
% scatter(100*correctSortedByPower', areasSortedbyPower', 50, "filled")
% legend(["0-1","1-2", "2-3"], 'location', 'northwest')
% plot_paper

