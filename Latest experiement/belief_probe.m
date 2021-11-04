% Demo for belief enquiry function
% p.display_type=0; debug = 1; configCogent(p,debug)

function [v] = belief_probe(debug, keys)
try 
    debug;
catch 
    debug = 1;
end 

%%  p structure has the parameters that the function will need.
% First, the keys we'll need to move around must be provided:

 p.UKey = 95;           % usu. up arrow.
 p.DKey = 100;          % usu. down arrow
 p.LKey = 97;
 p.RKey = 98;
 p.confirmKey   = keyCode('RETURN');

 % Then, a few parameters related to the graphics:

p.image_dir = 'images\'; 
p.imMCQbar = ['images\column.bmp'];
p.imMCQh = 400;
p.imMCQw = 120;
p.resolution = [1024 768]; 
p.HonDefForRetPic = [];        % sorry, this is relic ...

% Last, the option array - these are parametrized as beta beliefs.
load('BetaOp_5x5.mat');  p.barArr = Outp.pArr;  p.alpha  = Outp.abArr(:,:,1);   p.beta = Outp.abArr(:,:,2);

%% Now to use the function:

  % Labels for the belief bins
  v.lab = {'10-20','20-40','40-60','60-80','80-100'};
%   v.lab = {'<5/100','<15/100','<25/100%','75%','very likely(95%)'};

%   cgfont('Arial',20);   % Default font for this part
%   preparestring(['How likey is it that you will get a shock?'],1,0,200);
%   drawpict(1);
%   wait(3000);

%change true value depending on condition

  % this is the true index that the pt. is trying to guess
                           % Set this from 1 to 5 here.

  v = getOtherBinProb(p, v);            % Ask for index and confidence
%   belief_vector = v.belVec;
%   belief_parameters =  v.belPar;
%   belief_bars = v.currBars;
%   belief_time = v.tBelQ;
%   belief_fee = v.belVec(trueValue); 
 
  
  
 
  