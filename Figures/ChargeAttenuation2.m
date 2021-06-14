
attenuation = zeros(1,20);
capRatio = 20;
for i = 1:20
    attenuation(i) = (capRatio/(1+capRatio))^(20-i);
end

plot(attenuation);
plot_paper;
% set(gca, 'YScale', 'log')
xlabel("i")
ylabel("Attenuation factor due to additional accumulation")
title("Accumulation attenuation of different V_{in}[i] for C_1/C_2 = 20")