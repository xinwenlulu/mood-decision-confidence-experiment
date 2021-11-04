function tBarsShown = showBarQ( v, p,q1)
%
% v must contain:
%    v.currBars is a vector j, k with the coordinates of
%          p.barArr(:,j,k) e.g. j=confidence level, k=prefered bin
%          p.imMCQbar =  bmp with graphic to use for the bar
%    v.hd = main heading. Default: 'How ''forgiving'' was your partner this time ?'
%    v.lab = e.g. {'5%','25%','50%','75%','95%'}
% Optional to ensure non-clashing with existing stuff etc:
%          p.barTran  = transparent colour to use
%          p.barSpSt = starting index for Sprites for the bars


% SEE advAcorns.m and TST_mcq1.m for hints of how I did this sort of thing
% before ...

try % Sprite Start number for Bars
    p.BarSpSt;
catch
    p.BarSpSt = 600;
end
try 
    p.barTran;
catch
    p.barTran= 'm'; % ie [1, 0, 1] ie magenda to signify 'transparent'
end
try
   v.hd;
catch
   v.hd  = 'How ''forgiving'' was your partner this time?';
%    v.lab = {'less than 5%','15%','25%','35%','more than 50%'};
   v.lab = {'0-20','20-40','40-60','60-80','80-100'};

%    v.lab = {'unlikely <=1:20','1:8','1:4','1:3','likey >=1:2'};

end
try
   p.black;
catch
   p.black = [ 0 0 0];   p.dbrown=  [0.6  0.2  0];
end
try % see if honour (and so assume defect, forgive, retaliate) pictures provided
   p.HonDefForRetPic;
catch   
   p.HonDefForRetPic = [];
end
try
   p.image_dir;
catch
   p.image_dir = [];
end

try
   q;
catch
   q = [];
end

barTotN = size(p.barArr,2);
levTotN = size(p.barArr,3);
sp0 = p.BarSpSt ;   % shorthand for number of sprite zero, which doesn't exist
tableSp = sp0+1;     % next sprite for the table / base of the bars.

base = -100; 
colh = p.imMCQh;        colw = p.imMCQw;

% Make sprites. REM instructions are at page B26 on of the f*g c*t
% manual G2UsrManv131, e.g. in ...\Cogent\Cogent2000v1.31\Documents
% Remember the f*g cgsetsprite(0) before you draw ...

cgmakesprite(sp0,p.resolution(1),p.resolution(2),p.black);
cgmakesprite(tableSp,p.resolution(1),0.29*p.resolution(1),p.black);
cgtrncol(sp0,p.barTran);
cgsetsprite(sp0);
% cgrect(0,0,p.resolution(1),p.resolution(2)); % ,p.barTran);
cgloadbmp(sp0,p.imMCQbar,colw,colh); % ,(i-floor(barTotN/2))*(colw+20),z );

cgsetsprite(0)
% Draw bars
for i=1:barTotN
  h = p.barArr(i,v.currBars(1),v.currBars(2))*colh;
  z = base+h- round(colh/2);
  cgdrawsprite(sp0, (i-floor(barTotN/2)-1)*(colw+20),z );
end

% Draw table on which they rest
if isempty(p.image_dir)  % i.e. no image expected for table.
  cgpencol(p.black);  
  cgrect(0,base-p.resolution(2)/4-1,(barTotN+1)*(colw+20),p.resolution(2)/2);
  cgpencol(p.dbrown);  
  cgrect(0,base-21,(barTotN+1)*(colw+20),40);
  cgrect(0,base-256,(barTotN+1)*colw,500);
  cgpencol(p.black);  
  cgrect(0,base-256,(barTotN+1)*(colw-20),300);
  
else  % 'table.bmp' suitable image expected in p.image_dir
  try
      cgloadbmp(tableSp,[p.image_dir 'table.bmp']);
      cgdrawsprite(tableSp,0,base-0.145*p.resolution(1));
      cgloadbmp(tableSp+1,'2.bmp',150,150);
      cgdrawsprite(tableSp+1, 0, -250);
  catch
      error('No table.bmp found in p.image_dir');
  end
end
cgpencol(1,1,1);  

% Now draw the writing - main header w. its bits and labels of categories:
if ~isempty(p.HonDefForRetPic)  
    RpicSp = 102; LpicSp = 101; 
    cgmakesprite(RpicSp,p.resolution(1),p.resolution(2),p.black);
    cgmakesprite(LpicSp,p.resolution(1),p.resolution(2),p.black);

    %if v.trial > 1  % if at least one trial has been played before, 
       % show how participant 2 responded in the previous trial to use
       % in reminding the player. Rem prevResps is of the form
       % [honour forgive betray retaliate] x 2 rows
       NoD = v.trResps(2,2) + v.trResps(2, 4);
       NoC = v.trResps(2,1) + v.trResps(2, 3);
    %end
    
    if strcmp(v.hd,'How likely was your partner to HONOUR your offering?')
         %if v.trial > 1 
            remstr = ['Your partner honoured ' num2str(v.trResps(2,1)) ...
                  ' times after ' num2str(NoC) ' offers you made'];      
         %end     
        cgloadbmp(RpicSp, p.HonDefForRetPic{1}); 
        cgloadbmp(LpicSp, p.HonDefForRetPic{2}); 
        cgdrawsprite(RpicSp, 440, 280); cgdrawsprite(LpicSp, -440, 280); 
    elseif strcmp(v.hd,'How likely was your partner to RETALIATE your NOT offering?')
        %if v.trial > 1 
           remstr = ['Your partner retaliated ' num2str(v.trResps(2,4)) ...
                  ' times after ' num2str(NoD) ' NO-offers from you.'];    
        %end
        cgloadbmp(RpicSp, p.HonDefForRetPic{4}); 
        cgloadbmp(LpicSp, p.HonDefForRetPic{3});     
        cgdrawsprite(RpicSp, 440, 280); cgdrawsprite(LpicSp, -440, 280); 
    end
    
end

if ~ isempty(q1)
  fontSize = 40;
  cgfont('Arial',fontSize);
  %cgalign('l','t')
  cgtext(q1, 0, 200);
%   cgtext(v.hd, 0, base+colh+fontSize);
end 
   

try  % if there are fancy information / question strings, show them.
  fontSize = 30;
  cgfont('Arial',fontSize);
  cgtext(remstr, 0, base+colh+2*fontSize);
  cgtext(v.hd, 0, base+colh+fontSize);
end

for i=1:barTotN
      cgtext(v.lab{i}, (i-floor(barTotN/2)-1)*(colw+20),base-fontSize+5 );
end


tBarsShown = cgflip(p.black);

% Clear up a bit:
if ~isempty(p.HonDefForRetPic) 
    cgfreesprite(RpicSp); cgfreesprite(LpicSp); cgfreesprite(tableSp); 
end;

return;

