%% user input
clc; clear; close all;

%---ial Conditions-----------------------------------------------------
mE = 5.972e24; % Earth mass
m = [mE 0.7*mE 0.9*mE 0.8*mE, 0.75*mE, 0.85*mE];
r1 = [0;0;0];           v1 = [0;0;0];
r2 = [0;-9e8;9e8];      v2 = [100;-300;-900];
r3 = [13e8;-4e8;1e8];   v3 = [-300;-200;300];
r4 = [-12e8;5e7;-9e8];  v4 = [-800;200;400];
r5 = [1e9;1.2e9;-6e8];  v5 = [80;-900;50];
r6 = [-1e8;7e8;1.2e9];    v6 = [400;400;-500];
sv = [r1;v1;r2;v2;r3;v3;r4;v4;r5;v5;r6;v6];
% you can hardcode this. SI base units

%---Sim Times--------------------------------------------------------------
tEnd = seconds(years(2));           % End time
tScale = seconds(years(1/12));      % aim for 1 month per four second
dt = seconds(days(1/2));            % 2 frame per day
fps = 30;                     % framerate to aim for in animation
trailLifetime = seconds(years(10));   % How long trails should live. Set to inf to not erase them
rotrate = 15; % deg/sec rotate

useCoM = true;  % whether to use CoM as origin
useTex = true;  % whether to use tex labels (otherwise MATLAB labels)
useDarkMode = true; % whether to use dark mode
% Note: search for ANIMATIONCOMMENT for stuff to toggle based on whether
% animating or not

%% run sim

% came to me in a fever dream. Prevents the spacing between frames from
% being a non-integer
% Comes from finding the dt required to make spacing an integer. Done so by
% taking the expression that goes into sp and setting it equal to itself
% rounded
dt = (tEnd/(fps*tEnd/tScale))/ceil((tEnd/(fps*tEnd/tScale))/dt);

numbods = length(m);
if(numbods ~= length(sv)/6)
    error("Improper inputs. Check to make sure your state vector has 4x as many elements as your masses");
end

ops = odeset('RelTol', 1e-10, 'AbsTol', 1e-2);
[t,SV] = ode45(@(t,sv) dynamics(t,sv, m, numbods), 0:dt:tEnd, sv, ops);

%% Prepare plots
x = zeros(length(t),numbods);
y = zeros(length(t),numbods);
z = zeros(length(t),numbods);
vx = zeros(length(t),numbods);
vy = zeros(length(t),numbods);
vz = zeros(length(t),numbods);
for idx = 1:numbods
    x(:,idx) = SV(:,idx*6-5);
    y(:,idx) = SV(:,idx*6-4);
    z(:,idx) = SV(:,idx*6-3);
    vx(:,idx) = SV(:,idx*6-2);
    vy(:,idx) = SV(:,idx*6-1);
    vz(:,idx) = SV(:,idx*6);
end

if(useCoM)
    M=sum(m);
    xc=sum(x.*m,2)/M; yc=sum(y.*m,2)/M; zc=sum(z.*m,2)/M; vxc=sum(vx.*m,2)/M; vyc=sum(vy.*m,2)/M; vzc=sum(vz.*m,2)/M;
    x = x-xc;   y = y-yc;   z=z-zc; vx = vx-vxc;    vy = vy-vyc; vz=vz-vzc;
end
if(useTex)
    interp = "latex";
else
    interp = "tex";
end
H = linspace(0,1-1/numbods,numbods)';   % hue
colors = hsv2rgb([H, ones(size(H)), (.75+0.25*useDarkMode)*ones(size(H))]); %color order
%% Paths
figure('color','w');
colororder(colors);
plot3(x, y, z, '-'); ax=gca;
lgtxt = {};
for idx=1:numbods
    mExp = floor(log10(m(idx))); mCoef = round(m(idx)/10^mExp, 3);
    lgtxt{end+1} = "$$"+mCoef+"\times 10^{"+mExp+"}$$";
end
lg=legend(lgtxt,Interpreter=interp,Location="eastoutside", box="off");
lg.ItemTokenSize=[10 5]; lg.Title.String="Masses [kg]"; 
set(gca,'TickLabelInterpreter',interp)
title(numbods+"-Body Sim", Interpreter=interp); axis equal;grid on;

if(useDarkMode)
    set(gcf, "Color", 'k'); set(gca,'Color','k');
    set(gca,'GridColor','w');set(gca,'XColor','w');set(gca,'YColor','w');set(gca,'ZColor','w');
    lg.Title.Color='w'; ax.Title.Color='w'; set(lg,'textcolor','w')
end

[az, el] = view; xl = xlim; yl=ylim; zl=zlim;
%% Animation
figure; set(gcf, "color", 'w');
colororder(colors);
objs = plot3(xl',yl',zl', 'ko'); animAx = gca;
axis equal; grid;
title(numbods+"-Body Sim", Interpreter=interp); subtitle("$$t=0$$ years", Interpreter=interp);
set(gca,'TickLabelInterpreter',interp)
for idx=1:numbods
    al(idx) = animatedline("LineWidth", 0.1,"Color",colors(idx,:), "MaximumNumPoints",trailLifetime/dt);
end
lg=legend(["", lgtxt],Interpreter=interp,Location="eastoutside", box="off");
lg.ItemTokenSize=[10 5]; lg.Title.String="Masses [kg]";
view(az, el);
if(useDarkMode)
    set(gcf, "Color", 'k'); set(gca,'Color','k');
    set(gca,'GridColor','w');set(gca,'XColor','w');set(gca,'YColor','w');set(gca,'ZColor','w');
    lg.Title.Color='w'; animAx.Title.Color='w'; animAx.Subtitle.Color='w';
    set(lg,'textcolor','w'); set(objs,'color','w');
end

pause; tic; xlim(xl);ylim(yl);zlim(zl);
animation = VideoWriter("animation3d_"+string(datetime, "yyyy-MM-dd-hhmmss")+".mp4",'MPEG-4');      % ANIMATIONCOMMENT
animation.FrameRate = fps; animation.Quality=100; open(animation);      % ANIMATIONCOMMENT
sp = (length(t)-1)/(fps*tEnd/tScale); % frames (sim) per frame (display)
for idx=sp:sp:length(t)
    az = az + rotrate/fps;
    view(az, el);
    set(objs, "XData", x(idx,:), "YData", y(idx,:), "ZData", z(idx,:));
    for alIdx=1:length(al)
        addpoints(al(alIdx), x(idx+(1-sp:0),alIdx), y(idx+(1-sp:0),alIdx),  z(idx+(1-sp:0),alIdx));
    end
    animAx.Subtitle.String = sprintf("$$t=%.2f$$ years", years(seconds(t(idx))));
    %pause(max([1/fps-toc, 0]));tic;   % ANIMATIONCOMMENT
    writeVideo(animation, getframe(gcf));   % ANIMATIONCOMMENT
end
close(animation);   % ANIMATIONCOMMENT


%% Helper function
% it works, don't touch
function dsv = dynamics(t, sv, m, num)
G = 6.6743e-11;
dynMtx = zeros(6*num,6*num);
for targetnum = 1:num
    dynMtx((targetnum-1)*6+[1:3],(targetnum-1)*6+[4:6]) = eye(3);
    for sourcenum = 1:num
        d = norm(sv((1:3)+(sourcenum-1)*6) - sv((1:3)+(targetnum-1)*6));
        gmd3 = (G*m(sourcenum)/d^3)*eye(3);
        if(targetnum ~= sourcenum)
            dynMtx((targetnum-1)*6+[4:6],(sourcenum-1)*6+[1:3]) = gmd3;
            dynMtx((targetnum-1)*6+[4:6],(targetnum-1)*6+[1:3]) = dynMtx((targetnum-1)*6+[4:6],(targetnum-1)*6+[1:3])-gmd3;
        end
    end
end
dsv = dynMtx * sv;
end