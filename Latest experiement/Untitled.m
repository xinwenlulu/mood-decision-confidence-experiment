x = 0;
meanreward = []
for n = 1:size(rewardmHvH,2)
    x = x + rewardmHvH(n);
    meanreward(end + 1) = x/n;
end

plot(meanreward)