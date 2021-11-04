function Q = QLearn(Choices,LeftStim,RightStim,Outcome,Param)


Beta = Param(1);
LR =Param(2);
Drift = Param(3);


Len = length(Choices);

Q = 50*ones(Len,5);
PChoice= zeros(Len,3);

MySoftMax = @(x) exp(x)./sum(exp(x));

for i = 1:Len
PChoice(i,[LeftStim(i),RightStim(i)]) = MySoftMax(Beta*Q(i,[LeftStim(i),RightStim(i)]));   
Q(i+1,Choices(i))=Q(i,Choices(i))+LR*(Outcome(i)-Q(i,Choices(i)));
Unchosen = find([1,2,3]~=Choices(i));
Q(i+1,Unchosen(1))=Q(i,Unchosen(1))+LR*(Drift-Q(i,Unchosen(1)));
Q(i+1,Unchosen(2))=Q(i,Unchosen(2))+LR*(Drift-Q(i,Unchosen(2)));
Q(i+1,4) = Choices(i);
% Q(:,4) = chosen; Q(:,5) = unchosen
if LeftStim(i) == Choices(i)
    Q(i+1,5) = RightStim(i);
else
    Q(i+1,5) = LeftStim(i);
end
end