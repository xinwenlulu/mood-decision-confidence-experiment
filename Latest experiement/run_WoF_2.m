function [outcome] = run_WoF_Sascha(data, block)


%%

outcome = data.wofdraws
background = [0, 0, 0];
cgflip(data.background);

% %% Preamble
if ~strcmp( data.initials, 'debugging' )
    cgfont('Arial',data.FontSizeTask);
    cgpencol(1,1,1); cgtext('+', 0, 0);
    cgflip(data.background);
    %wait(7000)
    cgtext('Hold on.....', 0, 0);
    cgflip(data.background);
    wait(2000)
    if block==1
        cgtext(['You', '''re in with a huge chance...'], 0, 0);
    else
        cgtext(['You', '''re in with another chance...'], 0, 0);
    end
    cgflip(data.background);
    wait(3000)
end


if outcome == -1
    outcome = -6;
%     if block==2
%         outcome=-999;
%     end
elseif outcome == 1
    outcome = 5;
%     if block==2
%         outcome=999;
%     end
end

data.wofvalue = outcome;

%cgloadbmp(5,'WoF_01.bmp')
% cgloadbmp(10,'WoF_01.bmp'); w=500; h=500;
%
% cgdrawsprite(10,0,0)
% cgtext('|', 0, 270);
% cgtext('Get ready', 0, -300);
% cgflip(background);

%wait(3000)
%
%cgsetsprite(10)                          % = Cogent buffer 5
%[w,h]=loadpng('WoF_01.png',10,'n');    % 'g'=green is transparent
cgloadbmp(5,'Test_Wof_money.bmp'); w=500; h=500;

% disp('Generating rotated sprites...')
% tic
% for ang=1:359
%     buff=10+ang;
%     cgmakesprite(buff, w,h, [0 0 0])
%     cgsetsprite(buff)
%     cgrotatesprite(10,0,0,-ang)
% end
% toc
%

%Positions for gems at end
xgem = -200;
ygem = 310;


cgfont('Arial',data.FontSizeTask);
cgpencol(1,1,1);


cgdrawsprite(5,0,-50) %cgdrawsprite(5,0,0)
cgtext('|', 0, 170); %Originally 270
cgflip(background);


switch outcome
    case -6; dang = 15.69;
    case 5; dang = 15.83;
    case -5; dang = 14.5;
%     case 999; dang = 14.65;
end

ang = 2;


cgsetsprite(0)              % make sure to draw in sprite 0
cgdrawsprite(5,0,-100) %cgdrawsprite(5,0,0)
t=cgflip(background);
waituntil(t+5000)

cgdrawsprite(5,0,-100) %cgdrawsprite(5,0,0)
cgpencol(1,1,1);cgtext('|', 0, 170); %Originally 270
cgfont('Arial',data.FontSizeText);
cgtext('Press the space bar to spin the wheel !', 0, 350);
cgtext('You can win real money, but',0,300);
cgtext('you can also lose part of your payout.',0,250);
%cgpencol(0,0,0);
%cgrect(420,350,200,100);
%cgpencol(1,1,1)
%cgfont('Arial',data.FontSizeText);
%cgtext('Your Gems',420,365);
%cgtext(num2str(data.currentscore),420,330);

        % show extra points box
       if data.wofscore >= 0  
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('  £',400,330);        
        cgalign('l','c')
        cgtext(num2str(data.wofscore),428,330);
        cgalign('c','c');
      else
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('- £',400,330);        
        cgalign('l','c')
        cgtext(num2str(abs(data.wofscore)),428,330);
        cgalign('c','c');
      end

cgflip(background);
pause;
% waitkeydown(100000);

cgfont('Arial',80);
cgpencol(1,1,1);


while dang>0
    
    ang = ang + dang;
    if dang>10, x=0.1;
    elseif dang<=10 && dang>2, x=dang/100;
    elseif dang<=2 && dang >.05, x=dang/125;
    else x=0.005;
    end;
       
%     if                dang>8; x = 0.1;
%     elseif dang<=8 && dang>3; x = 0.05;
%     elseif dang<=3 && dang>1; x = 0.01;
%     else                      x = 0.005;
%     end
    dang = dang - x;
    
    cgsetsprite(0)              % make sure to draw in sprite 0
    %cgdrawsprite(5,0,0)
    cgdrawsprite([floor(rem(ang,360)) + double(~mod(floor(rem(ang,360)),2))]+10, 0,-100)
    cgfont('Arial',data.FontSizeTask)
    cgpencol(1,1,1); cgtext('|', 0, 170); %changed from 270
    cgflip(background);
end



%%

%wait(500);

%%
%cgdrawsprite(5,0,0)


cgfont('Arial',data.FontSizeTask);
if outcome > 0
    cgdrawsprite( floor(rem(ang,360))+10, 0,-100) %cgdrawsprite( floor(rem(ang,360))+10, 0,0)
    cgpencol(1,1,1); cgtext('|', 0, 170); %changed from 270
    cgpencol(0,1,0); % Green text
    cgtext('YOU WON !!!',0,235); %originally ,0,50
%     cgpencol(0,0,0);
%     cgrect(420,350,200,100);
%     cgpencol(1,1,1)
%     cgfont('Arial',data.FontSizeText)
%     cgtext('Your Gems',420,365);
%     cgtext(num2str(data.currentscore),420,330);
    
            % show Bonus box
      if data.wofscore >= 0  
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('  £',400,330);        
        cgalign('l','c')
        cgtext(num2str(data.wofscore),428,330);
        cgalign('c','c');
      else
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('- £',400,330);        
        cgalign('l','c')
        cgtext(num2str(abs(data.wofscore)),428,330);
        cgalign('c','c');
      end
    
    wait(750)
    cgflip(background)
    wait(750)
    playsound(15);
    playsound(16);
    
    for j=1:outcome/0.5,
        cgdrawsprite( floor(rem(ang,360))+10, 0,-100) %cgdrawsprite( floor(rem(ang,360))+10, 0,0)
        cgfont('Arial',data.FontSizeTask);
        cgpencol(1,1,1);cgtext('|', 0, 170); %changed from 270
        cgpencol(0,1,0); % Green text
        cgtext('YOU WON !!!',0,235); %originally ,0,50
        %%% Tally
%         cgpencol(0,0,0);
%         cgrect(420,350,200,100);
%         cgpencol(1,1,1)
%         cgfont('Arial',data.FontSizeText)
%         cgtext('Your Gems',420,365);
%         cgtext(num2str((data.currentscore)),420,330);
        data.wofscore = data.wofscore+0.5;
        
                % show extra points box
      if data.wofscore >= 0  
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('  £',400,330);        
        cgalign('l','c')
        cgtext(num2str(data.wofscore),428,330);
        cgalign('c','c');
      else
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('- £',400,330);        
        cgalign('l','c')
        cgtext(num2str(abs(data.wofscore)),428,330);
        cgalign('c','c');
      end
        
%         data.currentscore = data.currentscore+111;
        %%%
        img = data.winwofstim.picture; % Load reward stimuli (dimensions of which will be used for background of loss stimuli too)
        rewardstimDimensions = data.winwofstim.stimDimensions;
%         cgloadbmp(3,img,rewardstimDimensions(1),rewardstimDimensions(2))
%         cgdrawsprite(3,120,-128) %cgdrawsprite(3,120,-100)
%         cgtext(['+',num2str(outcome)],-70,-130); %originally -70,-30
        %draw the appropriate number of gems
        smallrewardDimensions = [40 40];
        cgloadbmp(4,img,smallrewardDimensions(1),smallrewardDimensions(2))
        for m=1:j, cgdrawsprite(4,xgem+(m-1)*45,ygem); end;
        cgflip(background);
        wait(200);
    end;
    


else
    cgdrawsprite( floor(rem(ang,360))+10, 0,-100) %cgdrawsprite( floor(rem(ang,360))+10, 0,0)
        cgpencol(1,1,1);
        cgtext('|', 0, 170); %changed from 270
        cgpencol(1,0,0); % Red
        cgtext('YOU LOST !! ',0,235); %Originally ,0,50
%         cgpencol(0,0,0);
%         cgrect(420,350,200,100);
%         cgpencol(1,1,1)
%         cgfont('Arial',data.FontSizeText)
%         cgtext('Your Gems',420,365);
%         cgtext(num2str(data.currentscore),420,330);
        
                % show extra points box
      if data.wofscore >= 0  
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('  £',400,330);        
        cgalign('l','c')
        cgtext(num2str(data.wofscore),428,330);
        cgalign('c','c');
      else
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('- £',400,330);        
        cgalign('l','c')
        cgtext(num2str(abs(data.wofscore)),428,330);
        cgalign('c','c');
      end
        
        wait(750)
        cgflip(background)
        wait(750)
        playsound(17);
        %playsound(18);
    for j=1:(-outcome)/0.5
        cgdrawsprite( floor(rem(ang,360))+10, 0,-100) %cgdrawsprite( floor(rem(ang,360))+10, 0,0)
        cgfont('Arial',data.FontSizeTask);
        cgpencol(1,1,1)
        cgtext('|', 0, 170); %changed from 270
        cgpencol(1,0,0); % Red
        cgtext('YOU LOST !! ',0,235); %Originally ,0,50
        %%% Tally
%         cgpencol(0,0,0);
%         cgrect(420,350,200,100);
%         cgpencol(1,1,1)
%         cgfont('Arial',data.FontSizeText)
%         cgtext('Your Gems',420,365);
%         cgtext(num2str((data.currentscore)),425,330);
        data.wofscore = data.wofscore-0.5;
        
                % show extra points box
      if data.wofscore >= 0  
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('  £',400,330);        
        cgalign('l','c')
        cgtext(num2str(data.wofscore),428,330);
        cgalign('c','c');
      else
        cgpencol(0,0,0);
        cgrect(420,350,200,100);
        cgpencol(1,1,1)
        cgfont('Arial',data.FontSizeText)
        cgtext('Bonus',420,365);
        cgtext('- £',400,330);        
        cgalign('l','c')
        cgtext(num2str(abs(data.wofscore)),428,330);
        cgalign('c','c');
      end
       
        %%%
        smallrewardDimensions = [40 40];
        img2 = data.losewofstim.picture;
        cgloadbmp(7,img2,smallrewardDimensions(1),smallrewardDimensions(2));
        for m=1:j, cgdrawsprite(7,xgem+(m-1)*45,ygem); end;
        %playsound(15);
        cgfont('Arial',data.FontSizeTask);
        cgflip(background);
        wait(200)
        end    
end
wait(5000);

cgfont('Arial',data.FontSizeText);
if data.wofscore == 5
    cgpencol(1,0.8,0); 
    cgtext('Congratulations! You just won £5!', 0, 150);
    cgtext('This bonus will be added to your payment at the end of the experiment!', 0, 0);
    cgflip(background);
    wait(data.time.instructiontime);
    cgtext('Congratulations! You just won £5!', 0, 150);
    cgtext('This bonus will be added to your payment at the end of the experiment!', 0, 0);
    cgtext('Click to continue',0,-200);
else
    cgpencol(0.6,0,0); 
    cgtext('Oops! You just lost £6!', 0, 0);
    %cgtext("This bonus will be added to your payment at the end of the experiment!", 0, 0);
    cgflip(background);
    wait(data.time.instructiontime);
    cgtext('Oops! You just lost £6!', 0, 0);
    cgtext('Click to continue',0,-200);
end
cgflip(background);
[~,~,~,mp]=cgmouse; mp=0; while mp==0, [~,~,~,mp]=cgmouse; end;
%%
cgflip(data.background);
wait(2000);