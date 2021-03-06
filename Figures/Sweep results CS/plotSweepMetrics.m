clear;
normal = [1,2];
interictal = [3,4];
ictal = 5;

snr_ssim_index = 1;  %% Naming = metric_optimizedForMetric
ssim_ssim_index = 2;
snr_snr_index = 3;
ssim_snr_index = 4;
snr_smooth_index = 5;
ssim_smooth_index = 6;

avgSNRNormal = zeros(1,100);
avgSNRInterictal = zeros(1,100);
avgSNRIctal = zeros(1,100);

avgSSIMNormal = zeros(3,100);
avgSSIMInterictal = zeros(3,100);
avgSSIMIctal = zeros(3,100);

avgSNR = zeros(3,100);
avgSSIM = zeros(3,100);

%% get SSIM optimized data

for simPoint = 1:100
    data = load("Metrics/"+num2str(simPoint)+".mat").full_output;
    SNRdata = data{1, snr_ssim_index};
    SSIMdata = data{1, ssim_ssim_index};
    avgSNRNormal(1, simPoint) = mean(mean(SNRdata(1,normal)));
    avgSNRInterictal(1, simPoint) = mean(mean(SNRdata(1,interictal)));
    avgSNRIctal(1, simPoint) = mean(mean(SNRdata(1,ictal)));
    avgSNR(1,simPoint) = mean(mean(SNRdata));
    
    avgSSIMNormal(1, simPoint) = mean(mean(SSIMdata(1,normal)));
    avgSSIMInterictal(1, simPoint) = mean(mean(SSIMdata(1,interictal)));
    avgSSIMIctal(1, simPoint) = mean(mean(SSIMdata(1,ictal)));
    avgSSIM(1,simPoint) = mean(mean(SSIMdata));
end
% 
figure
subplot(2,1,1)
hold on
plot(avgSNRNormal)
plot(avgSNRInterictal)
plot(avgSNRIctal)
l=legend(["Normal", "Interictal", "Ictal"], 'Location', "northwest");
xlabel("Simulation settings index");
ylabel("RSNR (dB)");
hold off
plot_paper


subplot(2,1,2)
hold on
plot(100*avgSSIMNormal(1,:))
plot(100*avgSSIMInterictal(1,:))
plot(100*avgSSIMIctal(1,:))
l=legend(["Normal", "Interictal", "Ictal"], 'Location', "southeast");
xlabel("Simulation settings index");
ylabel("SSIM (%)");
hold off
plot_paper

%% get SNR optimized data

for simPoint = 1:100
    data = load("Metrics/"+num2str(simPoint)+".mat").full_output;
    SNRdata = data{1, snr_snr_index};
    SSIMdata = data{1, ssim_snr_index};
    avgSNRNormal(2, simPoint) = mean(mean(SNRdata(1,normal)));
    avgSNRInterictal(2, simPoint) = mean(mean(SNRdata(1,interictal)));
    avgSNRIctal(2, simPoint) = mean(mean(SNRdata(1,ictal)));
    avgSNR(2,simPoint) = mean(mean(SNRdata));
    
    avgSSIMNormal(2, simPoint) = mean(mean(SSIMdata(1,normal)));
    avgSSIMInterictal(2, simPoint) = mean(mean(SSIMdata(1,interictal)));
    avgSSIMIctal(2, simPoint) = mean(mean(SSIMdata(1,ictal)));
    avgSSIM(2,simPoint) = mean(mean(SSIMdata));
end

% figure
% hold on
% plot(avgSNRNormal(2,:))
% plot(avgSNRInterictal(2,:))
% plot(avgSNRIctal(2,:))
% l=legend(["Normal", "Interictal", "Ictal"]);
% title(l, "SNR optimized");
% xlabel("Simulation settings index");
% ylabel("RSNR (dB)");
% hold off
% plot_paper
% 
% figure
% hold on
% plot(100*avgSSIMNormal(2,:))
% plot(100*avgSSIMInterictal(2,:))
% plot(100*avgSSIMIctal(2,:))
% l=legend(["Normal", "Interictal", "Ictal"]);
% title(l, "SNR optimized");
% xlabel("Simulation settings index");
% ylabel("SSIM (%)");
% hold off
% plot_paper


%% get smoothness optimized data

for simPoint = 1:100
    data = load("Metrics/"+num2str(simPoint)+".mat").full_output;
    SNRdata = data{1, snr_smooth_index};
    SSIMdata = data{1, ssim_smooth_index};
    avgSNRNormal(3, simPoint) = mean(mean(SNRdata(1,normal)));
    avgSNRInterictal(3, simPoint) = mean(mean(SNRdata(1,interictal)));
    avgSNRIctal(3, simPoint) = mean(mean(SNRdata(1,ictal)));
    avgSNR(3,simPoint) = mean(mean(SNRdata));
    
    avgSSIMNormal(3, simPoint) = mean(mean(SSIMdata(1,normal)));
    avgSSIMInterictal(3, simPoint) = mean(mean(SSIMdata(1,interictal)));
    avgSSIMIctal(3, simPoint) = mean(mean(SSIMdata(1,ictal)));
    avgSSIM(3,simPoint) = mean(mean(SSIMdata));
end

% figure
% hold on
% plot(avgSNRNormal(3,:))
% plot(avgSNRInterictal(3,:))
% plot(avgSNRIctal(3,:))
% l=legend(["Normal", "Interictal", "Ictal"]);
% title(l, "Smoothness optimized");
% xlabel("Simulation settings index");
% ylabel("RSNR (dB)");
% hold off
% plot_paper
% 
% figure
% hold on
% plot(100*avgSSIMNormal(3,:))
% plot(100*avgSSIMInterictal(3,:))
% plot(100*avgSSIMIctal(3,:))
% l=legend(["Normal", "Interictal", "Ictal"]);
% title(l, "Smoothness optimized");
% xlabel("Simulation settings index");
% ylabel("SSIM (%)");
% hold off
% plot_paper

figure
subplot(2, 1, 1)
hold on
plot(avgSNR')
plot([39, 39], [0, 25], 'k')
plot([69, 69], [0, 25], 'k')
text(15, 23, "M = 75")
text(50, 23, "M = 150")
text(79.5, 23, "M = 192") 
l=legend(["SSIM optimized", "SNR optimized", "Smooth"], 'Location', 'northwest');
title(l, "Algorithm setting");
xlabel("Simulation settings index");
ylabel("RSNR (dB)");
hold off
plot_paper


subplot(2, 1, 2)
hold on
plot(100*avgSSIM')
plot([39, 39], [0, 100], 'k')
plot([69, 69], [0, 100], 'k')
text(15, 95, "M = 75")
text(50, 95, "M = 150")
text(79.5, 95, "M = 192") 
l=legend(["SSIM optimized", "SNR optimized", "Smooth"], 'Location', 'northwest');
title(l, "Algorithm setting");
xlabel("Simulation settings index");
ylabel("SSIM (%)");
hold off
plot_paper

powers=1e6*load("simPointPower.mat").simPointPower;
areas=load("simPointAreas.mat").simPointAreas;


figure
hold on
scatter(avgSNR(1,1:36), areas(1:36), 50, 'filled')
scatter(avgSNR(1,37:68), areas(37:68), 50, 'filled')
scatter(avgSNR(1,69:end), areas(69:end), 50, 'filled')
xlabel("RSNR (dB)")
ylabel("Capacitance (#technology minimum capacitances)")
set(gca, "YScale", "log")
rico = (11888-3846)/(21.7147-15.0446);
l= plot([15.0446, 24], [3846, 3846+rico*(24-15.0446)+15+0.446], ".-r");
legend(["M = 75", "M=150", "M=192"], 'Location', 'southeast');
hold off
plot_paper


figure
hold on
scatter(100*avgSSIM(1,1:36), areas(1:36), 50, 'filled')
scatter(100*avgSSIM(1,37:68), areas(37:68), 50, 'filled')
scatter(100*avgSSIM(1,69:end), areas(69:end), 50, 'filled')
xlabel("SSIM (%)")
ylabel("Capacitance (#technology minimum capacitances)")
% set(gca, "YScale", "log")
legend(["M = 75", "M=150", "M=192"], 'Location', 'southeast');
hold off
plot_paper


% figure
% hold on

scatter(avgSNR(1,:), powers, 50, 'filled')
xlabel("SNR (dB)")
ylabel("Power consumption (\muW)")
% set(gca, "YScale", "log")
% rico = (11888-3846)/(21.7147-15.0446);
% l= plot([15.0446, 24], [3846, 3846+rico*(24-15.0446)+15+0.446], ".-r");
hold off
plot_paper


% figure
% hold on
% scatter(avgSNR(1,1:36), powers(1:36).*areas(1:36), 50, 'filled')
% scatter(avgSNR(1,37:68), powers(37:68).*areas(37:68), 50, 'filled')
% scatter(avgSNR(1,69:end), powers(69:end).*areas(69:end), 50, 'filled')
% xlabel("SNR (dB)")
% ylabel("Power \cdot Capacitance ( \muW \cdot # technology minimum capacitances )")
% rico = (36115.7-4252.59)/(22.0575-14.7474);
% l= plot([14.7474, 24], [4252.59, 15000+rico*(24-14.7474)], ".-r");
% set(gca, "YScale", "log")
% legend(["M = 75", "M=150", "M=192"], 'Location', 'southeast');
% hold off
% plot_paper

figure
hold on
scatter(avgSSIM(1,1:36), powers(1:36).*areas(1:36), 50, 'filled')
scatter(avgSSIM(1,37:68), powers(37:68).*areas(37:68), 50, 'filled')
scatter(avgSSIM(1,69:end), powers(69:end).*areas(69:end), 50, 'filled')
xlabel("SSIM (%)")
ylabel("Power \cdot Capacitance ( \muW \cdot # technology minimum capacitances )")
set(gca, "YScale", "log")
rico = (24611.8-4252.59)/(-0.649112+0.900869);
l= plot([0.649112, 0.99], [4252.59, 16000+rico*(0.99-0.649112)], ".-r");
legend(["M = 75", "M=150", "M=192"], 'Location', 'southeast');
hold off
plot_paper
