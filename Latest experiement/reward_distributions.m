function [reward_distribution]   = reward_distributions(mu,sigma,vectorLength)

%pd = makedist('Normal','mu',mu,'sigma',sigma); %input mean and standard deviation here
pd = makedist('Normal','mu',60,'sigma',25);
reward = [];

% randomly generate 10000 numbers from the normal distribution pd
for n = 1:10000 
    x = random(pd);
    reward(end+1) = x;
end

% sort the 10000 numbers in ascending orders
reward = sort(reward);


    reward_distribution = [];
    noise = randi([-50 50],1,15);
    while -100+645*15 + noise(end) > 9600
        noise = randi([-150 150],1,15);
    end
    
    for i = 1:15
        display(-100+645*i + noise(i))
            %take every 650th number from the 500th to the 9600th number in the distribution
           if x < 10000
              reward_distribution(end+1) = reward(-100+645*i + noise(i));
           end
    end
    
%round the numbers so they are integers
reward_distribution = round(reward_distribution);

%add a little noise to the lowest and highest value
reward_distribution(1) = reward_distribution(1) + randi([-1 1],1);
reward_distribution(end) = reward_distribution(end) + randi([-1 1],1);

%avoid values > 100 when mu = 60 and sigma = 25
for i = 1:15
    if reward_distribution(i) > 100
        reward_distribution(i) = 98 + randi([-1 1],1);
    end
end