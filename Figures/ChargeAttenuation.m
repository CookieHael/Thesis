
attenuation = zeros(1,500);
for capRatio = 1:500
    attenuation(capRatio) = 1/(1+.2*capRatio)*(.2*capRatio/(1+.2*capRatio))^9;
end

plot(attenuation);
plot_paper
xlabel("Ratio of C_1/C_2")
ylabel("Attenuation factor")
title("Attenuation of V_{in}[1] during 10 accumulations")