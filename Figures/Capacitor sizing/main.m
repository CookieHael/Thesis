clear;
F_per_M = 0.001025;  %F/m^2
C_mismatch_param_basic = 3.4878e-09;
k = 1.38064852e-23;
T = 293;
cap_gain_error = 1e-2;
cap_ratio = 15;
Vref = 2;

minimum_technology_cap = 1.995e-15;


DAC_caps1 = zeros(16,4);
C_mismatch_parameter = [1; 2; 10; 20]*C_mismatch_param_basic; %% Basic param is .25% mismatch, which is very little. Provided some larger ones for insights too

for i = 2:16
    for j=1:4
        
        C_mismatch = C_mismatch_parameter(j);
        number_of_bits = i;
        minimum_sampling_cap = 12*k*T*2^(2*number_of_bits)/Vref^2; 
        [minimum_match_cap_DAC, minimum_match_cap_MAC] = getCapSizes(number_of_bits,C_mismatch, minimum_technology_cap, cap_ratio, cap_gain_error, F_per_M);
        DAC_caps1(i,j) = max([ minimum_technology_cap, minimum_match_cap_DAC]);
        DAC_caps2(i,j) = max([ minimum_sampling_cap, minimum_technology_cap, minimum_match_cap_DAC]);
    end
end

% 
% figure('units','normalized','outerposition',[0 0 1 1])
% yyaxis left
% hold on
% plot(DAC_caps(:, 1))
% plot(DAC_caps(:,2), '-*')
% hold off
% ylabel("Unit capacitor size for linearity (F)")
% xlim([2,12])
% ylim([0.9*DAC_caps(2,1), 1.1*DAC_caps(12,2)])
% set(gca, "YScale", "log")
% plot_paper
% yyaxis right
% plot(DAC_caps(:,3:4))
% ylabel("Unit capacitor size for linearity (F)")
% ylim([0.9*DAC_caps(2,4), 1.1*DAC_caps(12,4)])
% set(gca, "YScale", "log")
% xlabel("Resolution in bits")
% legend(".25%", ".5%", "2.5%", "5%", "Location", 'northwest')
% plot_paper

fig=figure('units','normalized','outerposition',[0 0 1 1]);
colorord = get(gca, "colororder");
colors=colorord(1:4,:);
hold on
for i=1:4
plot(DAC_caps1(:,i), "-", "Color", colors(i,:))
end
plot_paper
ylabel("Unit capacitor size for matching (F)")
for i=1:4
plot(DAC_caps2(:,i), "--", "Color", colors(i,:))
end
set(gca, "YScale", "log")
plot_paper
hold off
xlim([3,16])
l=legend(".25%", ".5%", "2.5%", "5%","Location", "northwest");
title(l, "Mismatch")
xlabel("Resolution in bits")
ax=axes('Position',get(gca,'Position'),'Visible','Off');
line1 = line([0,0],[0,0],linestyle='-', color='k');
line2 = line([0,0],[0,0],linestyle='--', color='k');
l=legend(ax, [line1, line2], ["Mismatch only", "Mismatch and noise"], "Location", "northwest");
title(l, "Line type")
plot_paper


number_of_bits = 6;
cap_gain_errors = 1e-4:.5e-4:1e-2;
MAC_caps = zeros(length(cap_gain_errors),4);

for i=1:4
    
    C_mismatch = C_mismatch_parameter(i);
    for j=1:length(cap_gain_errors)
        cap_gain_error = cap_gain_errors(j);
        [minimum_match_cap_DAC, minimum_match_cap_MAC] = getCapSizes(number_of_bits,C_mismatch, minimum_technology_cap, cap_ratio, cap_gain_error, F_per_M);
        MAC_caps(j, i) = max([minimum_technology_cap, minimum_match_cap_MAC]);
    end
end
% 
figure('units','normalized','outerposition',[0 0 1 1])
MAC_caps = [fliplr(MAC_caps(:,1)),fliplr(MAC_caps(:,2)),fliplr(MAC_caps(:,3)),fliplr(MAC_caps(:,4))];
plot(100*cap_gain_errors,MAC_caps)
plot_paper
xlabel("Maximum capacitor gain error (%)")
ylabel("Capacitor size (F)")
set(gca, "YScale", "log")
legend(".25%", ".5%", "2.5%", "5%")





