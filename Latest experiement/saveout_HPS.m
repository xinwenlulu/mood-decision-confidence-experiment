clear
files=dir('AppyBrainPilot1*.mat')

for i=1:length(files)
    load(files(i).name)
    
    xxHPS(i,1) = data.id;
    xxRR(i,1)=data.id;
    if isfield(data, 'HPS')
        xxHPS(i,2) = data.HPS.HPS_score_mean;
        xxRR(i,2)= data.RR.RR_score_mean;
    else
        xxHPS(i,2) = nan;
        xxRR(i,2)=nan;
    end
end

for i=1: max(xxHPS(:,1))
    
   xHPS(i) = nanmean(xxHPS( xxHPS(:,1)==i, 2 ))
   xRR(i) = nanmean(xxRR( xxRR(:,1) ==i,2))
end

%%
x=xHPS*5;y=xRR;
x=x(~isnan(x)); y=y(~isnan(y));
figure;plot(x,y, '.','MarkerSize',12); hold on
coeffs = polyfit(x, y, 1);
fittedX = linspace(min(x), max(x), 200); fittedY = polyval(coeffs, fittedX);
plot(fittedX, fittedY, 'r-', 'LineWidth', 2);
[r, p]=corrcoef(x,y)