
function LL = QLearnDrift(Choices,LeftStim,RightStim,Outcome,Param)
% Choices are the choices made by participant 1,2 or 3
% LeftStim is which stimuli was on the left, RightStim is which stimuli was
% on the right
%Outcome is the reward at each timestep
% Param is a parameters vector [Beta,LR,Drift] to be estimated

Beta = Param(1);
LR =Param(2);
Drift = Param(3);

%initial parameters


Len=length(Choices);

Q = 50*ones(Len,3);
PChoice= zeros(Len,3);

MySoftMax = @(x) exp(x)./sum(exp(x));

for i = 1:Len
 
    % Calculate P values
    PChoice(i,[LeftStim(i),RightStim(i)]) = MySoftMax(Beta*Q(i,[LeftStim(i),RightStim(i)]));
     % Update Q for chosen option

    Q(i+1,Choices(i))=Q(i,Choices(i))+LR*(Outcome(i)-Q(i,Choices(i)));
     % Update Q for unchosen option 
    Unchosen = find([1,2,3]~=Choices(i));
    Q(i+1,Unchosen(1))=Q(i,Unchosen(1))+LR*(Drift-Q(i,Unchosen(1)));
    Q(i+1,Unchosen(2))=Q(i,Unchosen(2))+LR*(Drift-Q(i,Unchosen(2)));

end;
% Calculate log likelihood
PChoice(PChoice<0.001)=0.001;% avoid log(0)
LL=-sum(log(PChoice(Choices==1,1)))-sum(log(PChoice(Choices==2,2)))-sum(log(PChoice(Choices==3,3)));
