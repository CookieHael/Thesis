clear
load("SAR_Results.mat");
correctPercent = zeros(1,12);
uncertainCorrectPercent = zeros(1,12);
ind = 1:3:length(Correct);
for l=1:length(ind)
    i=ind(l);
    correctPercent(l+2) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainButCorrect(i+1)+Uncertain(i+2)*UncertainButCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
    uncertainCorrectPercent(l+2) = (2*Uncertain(i)*UncertainButCorrect(i)+2*Uncertain(i+1)*UncertainButCorrect(i+1)+Uncertain(i+2)*UncertainButCorrect(i+2))/5;
end

figure
subplot(2,2,1)
hold on
plot(100*correctPercent)
plot(100*uncertainCorrectPercent)
hold off
legend(["Correct", "Uncertain but correct"], 'location', 'southeast')
plot_paper
xlabel("ADC resolution (bit)")
ylabel("Correct detection (%)")
xlim([3, 12])



correctPercent2 = zeros(1,12);
uncertainCorrectPercent2 = zeros(1,12);
ind = 1:3:length(Correct);
for l=1:length(ind)
    i=ind(l);
    correctPercent2(l+2) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainButCorrect(i+1)+Uncertain(i+2)*UncertainButCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
    uncertainCorrectPercent2(l+2) = (2*Uncertain(i+1)*UncertainButCorrect(i+1)+Uncertain(i+2)*UncertainButCorrect(i+2))/3;
end

% subplot(2,1,2)
% hold on
% plot(correctPercent2)
% plot(uncertainCorrectPercent2)
% hold off
% plot_paper
% xlim([3, 12])

RSNR_avg = zeros(1,12);
RSNR_normal = zeros(1,12);
RSNR_interictal = zeros(1,12);
RSNR_ictal = zeros(1,12);


SSIM_avg = zeros(1,12);
SSIM_normal = zeros(1,12);
SSIM_ictal = zeros(1,12);
SSIM_interictal = zeros(1,12);

for i = 3:12
    metrics = load("SARMetrics/"+int2str(i));
    RSNR_avg(i) = mean(mean(metrics.RSNR_out(:,1:100)));
    SSIM_avg(i) = mean(mean(metrics.ssim_out(:,1:100)));
    
    RSNR_normal(i) = mean(mean(metrics.RSNR_out(1:2,1:100)));
    RSNR_interictal(i) = mean(mean(metrics.RSNR_out(3:4,1:100)));
    RSNR_ictal(i) = mean(mean(metrics.RSNR_out(5,1:100)));

    
    SSIM_normal(i) = mean(mean(metrics.ssim_out(1:2,1:100)));
    SSIM_ictal(i) = mean(mean(metrics.ssim_out(5,1:100)));
    SSIM_interictal(i) = mean(mean(metrics.ssim_out(3:4,1:100)));
end

subplot(2,2,2)
hold on
yyaxis left
ylabel("RSNR (dB)")
plot(RSNR_avg)
plot_paper
yyaxis right
plot(100*SSIM_avg)
ylabel("SSIM (%)")
xlabel("ADC Resolution (bit)")
hold off
xlim([3,12])
plot_paper

subplot(2,2,3)
hold on
plot(RSNR_normal)
plot(RSNR_interictal)
plot(RSNR_ictal)
hold off
ylabel("RSNR (dB)")
xlabel("ADC resolution (bit)")
xlim([3,12])
legend(["Normal", "Interictal", "Ictal"], 'location', 'southeast')
plot_paper

subplot(2,2,4)
hold on
plot(100*SSIM_normal)
plot(100*SSIM_interictal)
plot(100*SSIM_ictal)
hold off
ylabel("SSIM (%)")
xlabel("ADC resolution (bit)")
legend(["Normal", "Interictal", "Ictal"], 'location', 'southeast')
xlim([3,12])
plot_paper



