clear;

areas = load('simPointAreas.mat').simPointAreas;
powers = load('simPointPower.mat').simPointPower;

%% Plot SSIM, Smooth and SNR together
% % Get SSIM
% load("Detection results/SSIM_Results.mat");
% 
% % Calculate average correctness over all signals
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i)+2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5;
% end
% % Split correct into separate categories
% 
% normalCorrect = zeros(100,1);
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     normalCorrect(l)=Correct(i)+Uncertain(i)*UncertainCorrect(i);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% figure
% subplot(3,2,1)
% hold on
% plot(100*correctPercent)
% plot(100*uncertainCorrectPercent)
% plot([39, 39], [0, 100], 'k')
% plot([69, 69], [0, 100], 'k')
% text(15, 67, "M = 75")
% text(50, 95, "M = 150")
% text(79.5, 95, "M = 192") 
% title("SSIM-optimized", "FontSize", 30, "FontWeight", "bold");
% % xlabel("Simulation point index")
% ylabel("Correct detection (%)")
% legend(["Total correct", "Uncertain but correct"], 'Location', 'northwest')
% plot_paper
% 
% subplot(3,2,2)
% hold on
% plot(100*normalCorrect)
% plot(100*interictalCorrect)
% plot(100*ictalCorrect)
% plot([39, 39], [0, 100], 'k')
% plot([69, 69], [0, 100], 'k')
% text(15, 45, "M = 75")
% text(50, 90, "M = 150")
% text(79.5, 90, "M = 192")
% l=legend(["Normal", "Interictal", "Ictal"], 'Location', 'northwest');
% title(l, "Signal type")
% % xlabel("Simulation point index")
% ylabel("Correct detection (%)")
% plot_paper
% 
% % Get SNR
% load("Detection results/SNR_Results.mat");
% UncertainCorrect=UncertainButCorrect;
% 
% % Calculate average correctness over all signals
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i)+2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5;
% end
% max(correctPercent)
% % Split correct into separate categories
% 
% normalCorrect = zeros(100,1);
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     normalCorrect(l)=Correct(i)+Uncertain(i)*UncertainCorrect(i);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% subplot(3,2,3)
% hold on
% plot(100*correctPercent)
% plot(100*uncertainCorrectPercent)
% title("RSNR-optimized", "FontSize", 30, "FontWeight", "bold");
% % xlabel("Simulation point index")
% plot([39, 39], [0, 100], 'k')
% plot([69, 69], [0, 100], 'k')
% text(15, 95, "M = 75")
% text(50, 95, "M = 150")
% text(79.5, 95, "M = 192") 
% ylabel("Correct detection (%)")
% % legend(["Total correct", "Uncertain but correct"], 'Location', 'northwest')
% plot_paper
% 
% subplot(3,2,4)
% hold on
% plot(100*normalCorrect)
% plot(100*interictalCorrect)
% plot(100*ictalCorrect)
% plot([39, 39], [0, 100], 'k')
% plot([69, 69], [0, 100], 'k')
% text(15, 5, "M = 75")
% text(50, 5, "M = 150")
% text(79.5, 5, "M = 192")
% % l=legend(["Normal", "Interictal", "Ictal"], 'Location', 'northwest');
% % xlabel("Simulation point index")
% ylabel("Correct detection (%)")
% plot_paper
% 
% 
% % Get Smooth
% load("Detection results/Smooth_Results.mat");
% 
% 
% % Calculate average correctness over all signals
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i)+2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5;
% end
% 
% % Split correct into separate categories
% 
% normalCorrect = zeros(100,1);
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     normalCorrect(l)=Correct(i)+Uncertain(i)*UncertainCorrect(i);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% subplot(3,2,5)
% hold on
% plot(100*correctPercent)
% plot(100*uncertainCorrectPercent)
% title("Smoothness-optimized", "FontSize", 30, "FontWeight", "bold");
% xlabel("Simulation point index")
% ylabel("Correct detection (%)")
% plot([39, 39], [0, 100], 'k')
% plot([69, 69], [0, 100], 'k')
% text(15, 95, "M = 75")
% text(50, 95, "M = 150")
% text(79.5, 95, "M = 192") 
% % legend(["Total correct", "Uncertain but correct"], 'Location', 'northwest')
% plot_paper
% 
% subplot(3,2,6)
% hold on
% plot(100*normalCorrect)
% plot(100*interictalCorrect)
% plot(100*ictalCorrect)
% plot([39, 39], [0, 100], 'k')
% plot([69, 69], [0, 100], 'k')
% text(15, 5, "M = 75")
% text(50, 5, "M = 150")
% text(79.5, 5, "M = 192")
% % l=legend(["Normal", "Interictal", "Ictal"], 'Location', 'northwest');
% xlabel("Simulation point index")
% ylabel("Correct detection (%)")
% plot_paper


%% Now plot interictal+ictal correctness only
% 
% load("Detection results/SSIM_Results.mat");
% 
% % Calculate average correctness over all signals
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3;
% end

% figure
% hold on
% plot(100*correctPercent)

%     
% load("Detection results/SNR_Results.mat");
% UncertainCorrect=UncertainButCorrect;
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3;
% end

% plot(100*correctPercent)

% load("Detection results/Smooth_Results.mat");
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3;
% end

% plot(100*correctPercent)
% xlabel("Simulation point index")
% ylabel("Correct detection (%)")
% plot_paper
% l=plot([0, 100], [98, 98], "-.r");
% l2=plot([0, 100], [90, 90], "-.r");
% l3=plot([0, 100], [95, 95], "-.r");
% legend(["SSIM", "RSNR", "Smooth"], "Location", "northwest")
% plot_paper

%% Calculate FOM's for points that make the cut

% load("Detection results/SNR_Results.mat");
% UncertainCorrect=UncertainButCorrect;
% ind = 1:3:length(Correct);
% % Split correct into separate categories
% 
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     percent(i) = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
%     if percent(i)>90
%         FOM(i) = percent(i)/(powers(i)*areas(i));
%     else
%         FOM(i) = nan;
%     end
% end
% 
% figure
% subplot(3,1,1)
% scatter(1:100, FOM, 50, "filled")
% 
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     percent(i) = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
%     if percent(i)>95
%         FOM(i) = percent(i)/(powers(i)*areas(i));
%     else
%         FOM(i) = nan;
%     end
% end
% subplot(3,1,2)
% scatter(1:100, FOM, 50, "filled")
% 
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     percent(i) = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
%     if percent(i)>98
%         FOM(i) = percent(i)/(powers(i)*areas(i));
%     else
%         FOM(i) = nan;
%     end
% end
% subplot(3,1,3)
% scatter(1:100, FOM, 50, "filled")
% 
% load("Detection results/Smooth_Results.mat");
% 
% % Split correct into separate categories
% 
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     percent(i) = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
%     if percent(i)>90
%         FOM(i) = percent(i)/(powers(i)*areas(i));
%     else
%         FOM(i) = nan;
%     end
% end
% subplot(3,1,1)
% hold on
% ylabel("FOM")
% scatter(1:100, FOM, 50, "filled","d")
% l=legend(["RSNR-optimized", "Smooth"]);
% title(l, "90%")
% plot_paper
% 
% 
% FOM = zeros(100,1);
% for i = 1:100
%     percent = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
%     if percent>95
%         FOM(i) = percent/(powers(i)*areas(i));
%     else
%         FOM(i) = nan;
%     end
% end
% subplot(3,1,2)
% hold on
% ylabel("FOM")
% scatter(1:100, FOM, 50,"filled", "d")
% l=legend(["RSNR-optimized", "Smooth"]);
% title(l, "95%")
% plot_paper
% 
% 
% FOM = zeros(100,1);
% for i = 1:100
%     percent = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
%     if percent>98
%         FOM(i) = percent/(powers(i)*areas(i));
%     else
%         FOM(i) = nan;
%     end
% end
% subplot(3,1,3)
% hold on
% scatter(1:100, FOM, 50, "filled","d")
% ylabel("FOM")
% xlabel("Simulation point index")
% l=legend(["RSNR-optimized", "Smooth"],'location','southeast');
% title(l, "98%")
% plot_paper


%% Plot FOM's for all points
% load("Detection results/SNR_Results.mat");
% UncertainCorrect=UncertainButCorrect;
% FOM_SAR = 8.8e-6/(98*512);
% % Split correct into separate categories
% 
% interictalCorrect = zeros(100,1);
% ind = 1:3:length(Correct);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     percent(i) = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
% %     FOM(i) = percent(i)/(powers(i)*log2(areas(i)));
%       FOM(i) = powers(i)/(percent(i)*512);
% end
% 
% figure
% hold on
% scatter(percent(1:36), FOM(1:36), 50, "filled")
% scatter(percent(37:68), FOM(37:68), 50, 'd', "filled")
% scatter(percent(69:end), FOM(69:end), 50, 's', "filled")
% set(gca, "YScale", "log")
% plot([0, 100], [FOM_SAR, FOM_SAR], 'k')
% plot_paper


% 
% load("Detection results/Smooth_Results.mat");
% 
% % Split correct into separate categories
% 
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     percent(i) = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
%     FOM(i) = percent(i)/(powers(i)*areas(i));
% end
% 
% scatter(percent, FOM, 50, "d", "filled")
% 
% load("Detection results/SSIM_Results.mat");
% 
% % Split correct into separate categories
% 
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     percent(i) = 100*(2/3*interictalCorrect(i)+1/3*ictalCorrect(i));
%     FOM(i) = percent(i)/(powers(i)*areas(i));
% end
% 
% scatter(percent, FOM, 50, "s", "filled")
% 

%% Make quad figure for SNR
% clear;
% areas = load('simPointAreas.mat').simPointAreas;
% powers = load('simPointPower.mat').simPointPower;
% 
% load("Detection results/SNR_Results.mat");
% UncertainCorrect=UncertainButCorrect;
% 
% %Split correct into separate categories
% 
% %Calculate average correctness over all signals
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5;
% end
% %Split correct into separate categories
% 
% normalCorrect = zeros(100,1);
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     normalCorrect(l)=Correct(i)+Uncertain(i)*UncertainCorrect(i);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     FOM(i) = 100*correctPercent(i)/(powers(i)*areas(i));
% end
% 
% avgSNR=zeros(1,100);
% 
% avgSSIM = zeros(1,100);
% %Get SNR & SSIM
% for simPoint = 1:100
%     data = load("Metrics/"+num2str(simPoint)+".mat").full_output;
%     SNRdata = data{1, 3};
%     SSIMdata = data{1, 4};
%     avgSNR(1,simPoint) = mean(mean(SNRdata));
%     avgSSIM(1,simPoint) = mean(mean(SSIMdata));
% end
% 
% figure
% subplot(2,2,1)
% yyaxis left
% hold on
% ylabel("RSNR (dB)")
% scatter(100*correctPercent(1:36), avgSNR(1:36), 50, 'filled')
% scatter(100*correctPercent(37:68), avgSNR(37:68), 50, 'd', 'filled')
% scatter(100*correctPercent(69:end), avgSNR(69:end), 50, '+')
% plot_paper
% yyaxis right
% hold on
% ylabel("SSIM (%)")
% scatter(100*correctPercent(1:36), 100*avgSSIM(1:36), 50, 'filled')
% scatter(100*correctPercent(37:68), 100*avgSSIM(37:68), 50, 'd', 'filled')
% scatter(100*correctPercent(69:end), 100*avgSSIM(69:end), 50, '+')
% l1=scatter(nan, nan, 50,'k', 'filled');
% l2=scatter(nan,nan, 50,'k', "d", 'filled');
% l3=scatter(nan,nan, 50,'k', '+');
% legend([l1,l2,l3], ["M = 75", "M = 150", "M = 192"], 'location', 'northwest');
% plot_paper
% xlabel("Correct detection (%)")
% 
% subplot(2,2,2)
% hold on
% 
% 
% subplot(2,2,3)
% hold on
% scatter(100*avgSSIM(1:36), areas(1:36).*powers(1:36), 50, 'filled')
% scatter(100*avgSSIM(37:68), areas(37:68).*powers(37:68), 50, 'filled')
% scatter(100*avgSSIM(69:end), areas(69:end).*powers(69:end), 50, 'filled')
% hold off
% set(gca,'Yscale', 'log')
% plot_paper
% xlabel("SSIM (%)")
% legend(["M = 75", "M = 150", "M = 192"], 'location', 'southeast');
% ylabel("#Capacitors \cdot Power (# Unit cap \cdot \muW) ")
% 
% subplot(2,2,4)
% hold on
% scatter(avgSNR(1:36),  areas(1:36).*powers(1:36), 50, 'filled');
% scatter(avgSNR(37:68), areas(37:68).*powers(37:68), 50, 'filled');
% scatter(avgSNR(69:end), areas(69:end).*powers(69:end), 50, 'filled');
% plot_paper
% set(gca, "YScale", "log")
% plot_paper
% ylabel("#Capacitors \cdot Power (# Unit cap \cdot \muW) ")
% xlabel("RSNR (dB)")
% legend(["M = 75", "M = 150", "M = 192"], 'location','southeast');
% hold off

%% Make quad figure for SSIM
% 
% clear;
% areas = load('simPointAreas.mat').simPointAreas;
% powers = load('simPointPower.mat').simPointPower;
% 
% load("Detection results/SSIM_Results.mat");
% 
% % Split correct into separate categories
% 
% % Calculate average correctness over all signals
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5;
% end
% 
% % Split correct into separate categories
% 
% normalCorrect = zeros(100,1);
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     normalCorrect(l)=Correct(i)+Uncertain(i)*UncertainCorrect(i);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     FOM(i) = 100*correctPercent(i)/(powers(i)*areas(i));
% end
% 
% avgSNR=zeros(1,100);
% 
% avgSSIM = zeros(1,100);
% %Get SNR & SSIM
% for simPoint = 1:100
%     data = load("Metrics/"+num2str(simPoint)+".mat").full_output;
%     SNRdata = data{1, 1};
%     SSIMdata = data{1, 2};
%     avgSNR(1,simPoint) = mean(mean(SNRdata));
%     avgSSIM(1,simPoint) = mean(mean(SSIMdata));
% end
% 
% figure
% subplot(2,2,1)
% yyaxis left
% hold on
% ylabel("RSNR (dB)")
% scatter(100*correctPercent(1:36), avgSNR(1:36), 50, 'filled')
% scatter(100*correctPercent(37:68), avgSNR(37:68), 50, 'd', 'filled')
% scatter(100*correctPercent(69:end), avgSNR(69:end), 50, '+')
% plot_paper
% yyaxis right
% hold on
% ylabel("SSIM (%)")
% scatter(100*correctPercent(1:36), 100*avgSSIM(1:36), 50, 'filled')
% scatter(100*correctPercent(37:68), 100*avgSSIM(37:68), 50, 'd', 'filled')
% scatter(100*correctPercent(69:end), 100*avgSSIM(69:end), 50, '+')
% l1=scatter(nan, nan, 50,'k', 'filled');
% l2=scatter(nan,nan, 50,'k', "d", 'filled');
% l3=scatter(nan,nan, 50,'k', '+');
% legend([l1,l2,l3], ["M = 75", "M = 150", "M = 192"], 'location','southeast');
% plot_paper
% xlabel("Correct detection (%)")
% 
% subplot(2,2,2)
% hold on
% 
% 
% subplot(2,2,3)
% hold on
% scatter(100*avgSSIM(1:36), areas(1:36).*powers(1:36), 50, 'filled')
% scatter(100*avgSSIM(37:68), areas(37:68).*powers(37:68), 50, 'filled')
% scatter(100*avgSSIM(69:end), areas(69:end).*powers(69:end), 50, 'filled')
% hold off
% set(gca,'Yscale', 'log')
% plot_paper
% xlabel("SSIM (%)")
% legend(["M = 75", "M = 150", "M = 192"], 'location', 'southeast');
% ylabel("#Capacitors \cdot Power (# Unit cap \cdot \muW) ")
% 
% subplot(2,2,4)
% hold on
% scatter(avgSNR(1:36),  areas(1:36).*powers(1:36), 50, 'filled');
% scatter(avgSNR(37:68), areas(37:68).*powers(37:68), 50, 'filled');
% scatter(avgSNR(69:end), areas(69:end).*powers(69:end), 50, 'filled');
% plot_paper
% set(gca, "YScale", "log")
% plot_paper
% ylabel("#Capacitors \cdot Power (# Unit cap \cdot \muW) ")
% xlabel("RSNR (dB)")
% legend(["M = 75", "M = 150", "M = 192"], 'location','southeast');
% hold off


 %% Make quad figure for Smooth


clear;
areas = load('simPointAreas.mat').simPointAreas;
powers = load('simPointPower.mat').simPointPower;

load("Detection results/Smooth_Results.mat");

% Split correct into separate categories

% Calculate average correctness over all signals
correctPercent = zeros(1,100);
uncertainCorrectPercent = zeros(1,100);
ind = 1:3:length(Correct);
for l=1:length(ind)
    i=ind(l);
    correctPercent(l) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
    uncertainCorrectPercent(l) = (2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5;
end
% Split correct into separate categories

normalCorrect = zeros(100,1);
interictalCorrect = zeros(100,1);
ictalCorrect = zeros(100,1);
for l=1:length(ind)
    i=ind(l);
    normalCorrect(l)=Correct(i)+Uncertain(i)*UncertainCorrect(i);
    interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
    ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
end


FOM = zeros(100,1);
percent = zeros(100,1);
for i = 1:100
    FOM(i) = 100*correctPercent(i)/(powers(i)*areas(i));
end

avgSNR=zeros(1,100);

avgSSIM = zeros(1,100);
 % Get SNR & SSIM
for simPoint = 1:100
    data = load("Metrics/"+num2str(simPoint)+".mat").full_output;
    SNRdata = data{1, 5};
    SSIMdata = data{1, 6};
    avgSNR(1,simPoint) = mean(mean(SNRdata));
    avgSSIM(1,simPoint) = mean(mean(SSIMdata));
end

% figure
% subplot(2,2,1)
% yyaxis left
% hold on
% ylabel("RSNR (dB)")
% scatter(100*correctPercent(1:36), avgSNR(1:36), 50, 'filled')
% scatter(100*correctPercent(37:68), avgSNR(37:68), 50, 'd', 'filled')
% scatter(100*correctPercent(69:end), avgSNR(69:end), 50, '+')
% plot_paper
% yyaxis right
% hold on
% ylabel("SSIM (%)")
% scatter(100*correctPercent(1:36), 100*avgSSIM(1:36), 50, 'filled')
% scatter(100*correctPercent(37:68), 100*avgSSIM(37:68), 50, 'd', 'filled')
% scatter(100*correctPercent(69:end), 100*avgSSIM(69:end), 50, '+')
% l1=scatter(nan, nan, 50,'k', 'filled');
% l2=scatter(nan,nan, 50,'k', "d", 'filled');
% l3=scatter(nan,nan, 50,'k', '+');
% legend([l1,l2,l3], ["M = 75", "M = 150", "M = 192"], 'location', 'southeast');
% plot_paper
% xlabel("Correct detection (%)")
% 
% subplot(2,2,2)
% hold on
% 
% load('simPointsSettings.mat')
% subplot(2,2,3)
% hold on
% scatter(100*avgSSIM(1:36), 1e6*areas(1:36).*powers(1:36), 50, 'filled')
% scatter(100*avgSSIM(37:68), 1e6*areas(37:68).*powers(37:68), 50, 'filled')
% scatter(100*avgSSIM(69:end), 1e6*areas(69:end).*powers(69:end), 50, 'filled')
% hold off
% set(gca,'Yscale', 'log')
% plot_paper
% xlabel("SSIM (%)")
% legend(["M = 75", "M = 150", "M = 192"], 'location', 'southeast');
% ylabel("#Capacitors \cdot Power (# Unit cap \cdot \muW) ")
% 
figure
hold on
scatter(1./(10.^(avgSNR(1:36)/10)),  1e6*areas(1:36).*powers(1:36), 50, 'filled');
scatter(1./(10.^(avgSNR(37:68)/10)), 1e6*areas(37:68).*powers(37:68), 50, 'filled');
scatter(1./(10.^(avgSNR(69:end)/10)), 1e6*areas(69:end).*powers(69:end), 50, 'filled');
plot_paper
% set(gca, "YScale", "log")
plot_paper
ylabel("#Capacitors \cdot Power (# Unit cap \cdot \muW) ")
xlabel("RSNR (dB)")
legend(["M = 75", "M = 150", "M = 192"], 'location','southeast');
hold off

% load("simPointsSettings.mat")
% figure
% hold on
% scatter(100*avgSSIM(1:36), 1e6*areas(1:36).*powers(1:36), 50, 'filled')
% for i=1:36
%     text(100*avgSSIM(i)+.2,  1e6*(areas(i).*powers(i)+.025*areas(i).*powers(i)), num2str(adcResolutions(i)));
%     text(100*avgSSIM(i)+.2,  1e6*(areas(i).*powers(i)), num2str(leakageBoosts(i)/2), 'Color', 'red');
%     text(100*avgSSIM(i)+.2,  1e6*(areas(i).*powers(i)-.025*areas(i).*powers(i)), num2str(samplingGains(i)), "Color", "blue");
% end
% set(gca,'Yscale', 'log')
% % scatter(100*avgSSIM(37:68), areas(37:68).*powers(37:68), 50, 'filled')
% % scatter(100*avgSSIM(69:end), areas(69:end).*powers(69:end), 50, 'filled')
% plot_paper
% xlabel("SSIM (%)")
% l1=plot([nan, nan], [nan,nan], "k");
% l2=plot([nan, nan], [nan,nan], "r");
% l3=plot([nan, nan], [nan,nan], "b");
% legend(["", "Resolution", "Leakage size boost", "Sampling gain"], 'location', 'northwest');
% hold off
% ylabel("#Capacitors \cdot Power (# Unit cap \cdot \muW) ")
% ylim([1.6e3, 9.8e3])
% 
% figure
% hold on
% scatter(100*avgSSIM(1:36), areas(1:36), 50, 'filled')
% for i=1:36
%     text(100*avgSSIM(i)+.2,  areas(i)+0.015*areas(i), num2str(adcResolutions(i)));
%     text(100*avgSSIM(i)+.2,  areas(i), num2str(leakageBoosts(i)/2), 'Color', 'red');
%     text(100*avgSSIM(i)+.2,  areas(i)-.015*areas(i), num2str(samplingGains(i)), "Color", "blue");
% end
% set(gca,'Yscale', 'log')
% % scatter(100*avgSSIM(37:68), areas(37:68).*powers(37:68), 50, 'filled')
% % scatter(100*avgSSIM(69:end), areas(69:end).*powers(69:end), 50, 'filled')
% plot_paper
% xlabel("SSIM (%)")
% l1=plot([nan, nan], [nan,nan], "k");
% l2=plot([nan, nan], [nan,nan], "r");
% l3=plot([nan, nan], [nan,nan], "b");
% legend(["", "Resolution", "Leakage size boost", "Sampling gain"], 'location', 'northwest');
% hold off
% ylabel("#Capacitors (# Unit cap) ")
% 
% figure
% hold on
% scatter(100*avgSSIM(1:36), 1e6*powers(1:36), 50, 'filled')
% for i=1:36
%     text(100*avgSSIM(i)-.1,  1e6*(powers(i)+0.06*powers(i)), num2str(adcResolutions(i)));
%     text(100*avgSSIM(i)-.1,  1e6*powers(i)+1e6*0.04*powers(i), num2str(leakageBoosts(i)/2), 'Color', 'red');
%     text(100*avgSSIM(i)-.1,  1e6*(powers(i)+0.02*powers(i)), num2str(samplingGains(i)), "Color", "blue");
% end
% set(gca,'Yscale', 'log')
% % scatter(100*avgSSIM(37:68), areas(37:68).*powers(37:68), 50, 'filled')
% % scatter(100*avgSSIM(69:end), areas(69:end).*powers(69:end), 50, 'filled')
% plot_paper
% xlabel("SSIM (%)")
% l1=plot([nan, nan], [nan,nan], "k");
% l2=plot([nan, nan], [nan,nan], "r");
% l3=plot([nan, nan], [nan,nan], "b");
% legend(["", "Resolution", "Leakage size boost", "Sampling gain"], 'location', 'northwest');
% hold off
% ylabel("Power (\muW) ")
% ylim([.95, 2.8])
% 


%% Test smth
% clear;
% load("simPointsSettings.mat");
% areas = load('simPointAreas.mat').simPointAreas;
% powers = load('simPointPower.mat').simPointPower;
% 
% load("Detection results/SNR_Results.mat");
% 
% UncertainCorrect=UncertainButCorrect;
% % Split correct into separate categories
% 
% % Calculate average correctness over all signals
% correctPercent = zeros(1,100);
% uncertainCorrectPercent = zeros(1,100);
% ind = 1:3:length(Correct);
% for l=1:length(ind)
%     i=ind(l);
%     correctPercent(l) = (2*Correct(i+1)+Correct(i+2) + 2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/3; %% Twice as many normal and interictal signals as ictal signals
%     uncertainCorrectPercent(l) = (2*Uncertain(i)*UncertainCorrect(i)+2*Uncertain(i+1)*UncertainCorrect(i+1)+Uncertain(i+2)*UncertainCorrect(i+2))/5;
% end
% % Split correct into separate categories
% 
% normalCorrect = zeros(100,1);
% interictalCorrect = zeros(100,1);
% ictalCorrect = zeros(100,1);
% for l=1:length(ind)
%     i=ind(l);
%     normalCorrect(l)=Correct(i)+Uncertain(i)*UncertainCorrect(i);
%     interictalCorrect(l)=Correct(i+1)+Uncertain(i+1)*UncertainCorrect(i+1);
%     ictalCorrect(l)=Correct(i+2)+Uncertain(i+2)*UncertainCorrect(i+2);
% end
% 
% 
% FOM = zeros(100,1);
% percent = zeros(100,1);
% for i = 1:100
%     FOM(i) = 100*correctPercent(i)/(powers(i)*areas(i));
% end
% 
% 
% figure
% hold on
% scatter(100*correctPercent(1:36), 1e6*areas(1:36).*powers(1:36), 50, 'filled')
% for i=1:36
%     text(100*correctPercent(i)+.2,  1e6*(areas(i).*powers(i)+.025*areas(i).*powers(i)), num2str(adcResolutions(i)));
%     text(100*correctPercent(i)+.2,  1e6*(areas(i).*powers(i)), num2str(leakageBoosts(i)/2), 'Color', 'red');
%     text(100*correctPercent(i)+.2,  1e6*(areas(i).*powers(i)-.025*areas(i).*powers(i)), num2str(samplingGains(i)), "Color", "blue");
% end
% set(gca,'Yscale', 'log')
% % scatter(100*avgSSIM(37:68), areas(37:68).*powers(37:68), 50, 'filled')
% % scatter(100*avgSSIM(69:end), areas(69:end).*powers(69:end), 50, 'filled')
% plot_paper
% xlabel("SSIM (%)")
% l1=plot([nan, nan], [nan,nan], "k");
% l2=plot([nan, nan], [nan,nan], "r");
% l3=plot([nan, nan], [nan,nan], "b");
% legend(["", "Resolution", "Leakage size boost", "Sampling gain"], 'location', 'northwest');
% hold off
% ylabel("#Capacitors \cdot Power (# Unit cap \cdot \muW) ")
% ylim([1.6e3, 9.8e3])
% 
% figure
% hold on
% scatter(100*correctPercent(1:36), areas(1:36), 50, 'filled')
% for i=1:36
%     text(100*correctPercent(i)+.2,  areas(i)+0.015*areas(i), num2str(adcResolutions(i)));
%     text(100*correctPercent(i)+.2,  areas(i), num2str(leakageBoosts(i)/2), 'Color', 'red');
%     text(100*correctPercent(i)+.2,  areas(i)-.015*areas(i), num2str(samplingGains(i)), "Color", "blue");
% end
% set(gca,'Yscale', 'log')
% % scatter(100*avgSSIM(37:68), areas(37:68).*powers(37:68), 50, 'filled')
% % scatter(100*avgSSIM(69:end), areas(69:end).*powers(69:end), 50, 'filled')
% plot_paper
% xlabel("SSIM (%)")
% l1=plot([nan, nan], [nan,nan], "k");
% l2=plot([nan, nan], [nan,nan], "r");
% l3=plot([nan, nan], [nan,nan], "b");
% legend(["", "Resolution", "Leakage size boost", "Sampling gain"], 'location', 'northwest');
% hold off
% ylabel("#Capacitors (# Unit cap) ")
% 
% figure
% hold on
% scatter(100*correctPercent(1:36), 1e6*powers(1:36), 50, 'filled')
% for i=1:36
%     text(100*correctPercent(i)-.1,  1e6*(powers(i)+0.06*powers(i)), num2str(adcResolutions(i)));
%     text(100*correctPercent(i)-.1,  1e6*powers(i)+1e6*0.04*powers(i), num2str(leakageBoosts(i)/2), 'Color', 'red');
%     text(100*correctPercent(i)-.1,  1e6*(powers(i)+0.02*powers(i)), num2str(samplingGains(i)), "Color", "blue");
% end
% set(gca,'Yscale', 'log')
% % scatter(100*avgSSIM(37:68), areas(37:68).*powers(37:68), 50, 'filled')
% % scatter(100*avgSSIM(69:end), areas(69:end).*powers(69:end), 50, 'filled')
% plot_paper
% xlabel("SSIM (%)")
% l1=plot([nan, nan], [nan,nan], "k");
% l2=plot([nan, nan], [nan,nan], "r");
% l3=plot([nan, nan], [nan,nan], "b");
% legend(["", "Resolution", "Leakage size boost", "Sampling gain"], 'location', 'northwest');
% hold off
% ylabel("Power (\muW) ")
% ylim([.95, 2.8])
