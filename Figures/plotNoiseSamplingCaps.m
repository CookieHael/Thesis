% for x = 10:50
%      numbers = zeros(1,x+1);
%     for i = 0:x
%         numbers(i+1) = (x/(1+x))^(2*i);
%     end
%     factor = (1/(1+x)+ (x/(1+x))/x^2);
%     output(x-9) = factor*sum(numbers);
% end
%plot(output);
   
figure
hold on
c = zeros(10,40);
for j=1:10
    C_1 = 1/j;
    for x = 10:50
        C_2 = 1/x;
        partials = zeros(1,30);
        for k = 0:31
            partials(k+1) = (C_2/(C_1+C_2))^(2*k);
        end
        factor = (C_1/(C_1 + C_2) + (C_1*C_2/(C_1+C_2))/(C_2)^2);
        output(j, x-9) = factor*sum(partials);
    end
end
surf(output)
