function var = getOtherBinProb(par,var,debug)
%  function var = getOtherBinProb(par,var,debug)
%  var to provide labels for the belief bins:
%       v.lab = {'little (5%)','25%','50%','75%','a lot (95%)'};
%  par must provide:
%       p.barArr: 3-D array with belief distros to choose from [:,peakIndex, confidenceIndex]
%       e.g. ready-made belief histograms:
%            load('BetaOp_5x5.mat');  p.barArr = Outp.pArr;  p.alpha  = Outp.abArr(:,:,1);   p.beta = Outp.abArr(:,:,2);
%

%% a little initialization:
catTotN = size(par.barArr,2); % Number of categories over which to express beliefs
var.belVec = -666*ones(1,catTotN);  % the output belief vector initialized to 'invalid' value
var.belPar = -666*ones(1,2);        % the outp. bel. vec. params
var.currBars = [2, 3];   % [ belief confidence level index, peak belief bin index]
maxInd = [size(par.barArr,2), size(par.barArr,3)];
ansStep = 1;
%%
try
    debug;
catch
    debug = 0;
end
if debug ~= 0 
  %% Don't use cogent - just the command line:
  %entryISlice=2;  tTrStart=100;
  % [entryISlice tTrStart] =  lastslice(par.port);
  
  while ProbAns < 0 || ProbAns > size(par.Cgx,2)
         
    % present the state of the system:
    %  stateStr = ['So far you have ',num2str(var.assets(1)),' points. '];

    qString = [stateStr,['Please select most probable level, from 1 (low) to '],size(par.Cgx,2),' >> '];
    disp ' ';
    % Now collect participant choice:
    Sp = input(qString);
    if isempty(Sp) && isempty(str2num(Sp)) % need to check this separately, and first.
        disp(' ');
        disp('No valid input detected - please try again:');
    else
        % Record participant valid choices. Invalid choices will simply
        % prolong the apparent response time (until a valid one is chosen!)
        ProbAns = Sp;
        
        %decisionSlice=10; tChoice=500;               
    end % end if a valid action entered
    
    error('Here we need to ask about confidence too ...');
    
  end % end while no valid action yet.
  %%
else % if not debugging, we have cogent and (emul)scan running ...
  %% Display all  in cogent:

  cgflip(0,0,0);                cgflip(0,0,0);
  % Announce the new round:
%   cgfont('Arial',20);   % Default font for this part
%   clearpict(1); 
%   preparestring(['How likey is it that you will get a shock?'],1,0,0);
%   drawpict(1);
%   wait(3000);
 
   q1 = sprintf(['Estimate the range of gems you can get from this flower ']);

%   q = ['What are the odds of you being shocked?'];
 

  % Record the time etc. of presentation just BEFORE FIRST state presentation. 
  %[entryISlice tTrStart] =  lastslice(par.port);
  %tTrStart = tTrStart / par.tUnitsPerSec - dat.expStartTime; % convert to seconds relative to start of expt.
  
  while (var.belVec(1) == -666) 
    key = -666;
    var.tBelQ  = [-666, -666];
    timestart = time;
    while key ~= par.confirmKey
      showBarQ( var, par,q1);      
      if var.tBelQ(1) == -666; var.tBelQ(1)=now; end; 
        
      [key, keyt, n] = waitkeydown(5000, [par.RKey, par.LKey, ...
                                         par.UKey, par.DKey, par.confirmKey]);
                                     
      if isempty(key)
          key = 59;
      end 
      key = key(end);                              
      % Adjust preference bin:
      if key == par.RKey && (var.currBars(2) +ansStep <= maxInd(2)  )
          var.currBars(2) = var.currBars(2)+ansStep; 
      elseif key == par.LKey && (var.currBars(2)-ansStep >= 1)
          var.currBars(2) = var.currBars(2)-ansStep;
      end
      
      % Adjust confidence level: 
      if key == par.UKey && (var.currBars(1) +ansStep <= maxInd(1)  )
          var.currBars(1) = var.currBars(1)+ansStep; 
      elseif key == par.DKey && (var.currBars(1)-ansStep >= 1)
          var.currBars(1) = var.currBars(1)-ansStep;
      end
      
      keytime =  time;

      %if you only want to display for a period of 4 seconds 
%        if keytime - timestart > 4000
%             key = par.confirmKey;
%        end


   end % whil  e MCQ ans not confirmed
         
    var.tBelQ(2) = now;  % record time belief MCQ confirmed
    % params of a beta pdf that 'bins' similar to pt. choice:
    var.belPar(1) = par.alpha(var.currBars(1),var.currBars(2));
    var.belPar(2) = par.beta(var.currBars(1),var.currBars(2));
    % the histogram / probability vector that the pt. chose:
    var.belVec = par.barArr(:,var.currBars(1),var.currBars(2));
    var.probestart = timestart;
    var.keytime = keytime; 

    %[decisionSlice, tChoice] =   lastslice(par.port);
    %tChoice = tChoice / par.tUnitsPerSec - dat.expStartTime; % convert to seconds relative to start of expt.

  %%
  end

cgflip(0,0,0);                    cgflip(0,0,0);
%% ********************************* End of function *******************************
end
%%
