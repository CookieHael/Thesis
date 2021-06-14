function [correctPercent2,power2] = prepFit(correctPercent, power)
%Remove non-minimum data to fit trend
correctPercent2=[];
power2=[];
for i=1:length(correctPercent)
    currentP = power(i);
    currentC = correctPercent(i);
    flag = 0;
    for j=1:length(correctPercent)
        if (j~=i)
            if (power(j)<=currentP && correctPercent(j)>=currentC)
                flag = 1;
            end
        end
    end
    if (~flag)
        correctPercent2 = [correctPercent2; currentC];
        power2 = [power2; currentP];
    end
end

