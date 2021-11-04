function [ d v ] = guessOtherPolH( p, d, v )
%
%   Ask for inference on other-Honour policy only 

  cgfont('Arial',30);   % font for main question

  % first for CafterC (hon) :

  v.lab = {'little (5%)','25%','50%','75%','a lot (95%)'};
  v.hd  = 'How likely was your partner to HONOUR your offering?';
  v = getOtherBinProb(p, v);
  v.trialFee = v.belVec(v.PpolInd(1)) * p.OinfCoeff;
  d.evo(p.QRow(v.trial),14:15) = v.belPar;    % 'infOCCa',  'infOCCb' in prepDatIPD1

  % Add to actual game winnings:
  v.trialFee = v.trialFee + sum(p.roundEnd + v.trRets(:,1))*p.trWinCoeff;
  d.evo (p.QRow(v.trial),22) = v.trialFee;
  d.key_data(p.keySt(v.trial), 10)= v.trialFee;
  
end % of fn guessOtherPolH
