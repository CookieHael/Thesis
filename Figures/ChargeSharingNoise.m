clear;
k = 1.38064852e-23;
T = 293;


nb_acc = 10;
C1=1.995e-15;
noise_sigma = zeros(30,10);
nb_acc_array = [5, 10, 15, 20, 25];

for z = 1:5
    for l=1:30
        nb_acc = nb_acc_array(z);
        C2=l*C1;

        sum=0;

        for i = 0:nb_acc
            sum = sum + (C2/(C1+C2))^(2*i);
        end


        noise_sigma(l,z) = ((k*T * C1/(C1+C2)^2) + (k*T*(C1*C2/(C1+C2))/C2^2)) * sum;
    end
end
C1=10*C1;
for z = 6:10
    for l=1:30
        nb_acc = nb_acc_array(z-5);
        C2=l*C1;

        sum=0;

        for i = 0:nb_acc
            sum = sum + (C2/(C1+C2))^(2*i);
        end


        noise_sigma(l,z) = ((k*T * C1/(C1+C2)^2) + (k*T*(C1*C2/(C1+C2))/C2^2)) * sum;
    end
end
figure
colorord = get(gca, "colororder");
colors=colorord(1:5,:);
hold on
for i = 1:5
plot(abs(10*log10(noise_sigma(:,i))),"Color", colors(i,:))
end
for i = 1:5
plot(abs(10*log10(noise_sigma(:,i+5))), "--","Color", colors(i,:))
end
plot_paper
xlabel("C2/C1")
l=legend(["5", "10", "15", "20","25"], "Location", "northwest");
title(l, "#Accumulations")
xlim([1, 30])
ylim([abs(10*log10(noise_sigma(1,1))), abs(10*log10(noise_sigma(30,6)))])
ylabel("Noise floor in SNR (dB)")
ax=axes('Position',get(gca,'Position'),'Visible','Off');
line1 = line([0,0],[0,0],linestyle='-', color='k');
line2 = line([0,0],[0,0],linestyle='--', color='k');
l=legend(ax, [line1, line2], ["C1 = 1.995 fF", "C1 = 19.95 fF"], "Location", "northwest");
title(l, "Line type")
plot_paper

