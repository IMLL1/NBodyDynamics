%% user input
clc; clear; close all;

mSun = 1.989e30;
m = [mSun/3 mSun/3 mSun/3 mSun/3];
r1init = [1e9;0]; r2init = [0;1e8]; r3init = [-1e9;0];  r4init = [0;-1e8];
v1init = [0;-5]; v2init = [4.6;0];  v3init = [0;5]; v4init = [-4.6;0];
svinit = [r1init;v1init;r2init;v2init;r3init;v3init;r4init;v4init];
% you can hardcode this. Base distance unit is km

tEnd = seconds(years(50));    % sim for 50 years
tScale = seconds(years(2));   % aim for a year per second
dt = seconds(years(1/200));   % 200 frames per year
fps = 30;

useCoM = true;  % whether to use CoM as origin
useTex = true;  % whether to use tex labels (otherwise MATLAB labels)
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
if(numbods ~= length(svinit)/4)
    error("Improper inputs. Check to make sure your state vector has 4x as many elements as your masses");
end

ops = odeset('RelTol', 1e-10, 'AbsTol', 1e-5);
[t,SV] = ode45(@(t,sv) dynamics(t,sv, m, numbods), 0:dt:tEnd, svinit, ops);

%% Prepare plots
x = zeros(length(t),numbods);
y = zeros(length(t),numbods);
vx = zeros(length(t),numbods);
vy = zeros(length(t),numbods);
for idx = 1:numbods
    x(:,idx) = SV(:,idx*4-3);
    y(:,idx) = SV(:,idx*4-2);
    vx(:,idx) = SV(:,idx*4-1);
    vy(:,idx) = SV(:,idx*4);
end

if(useCoM)
    M=sum(m);
    xc=sum(x.*m,2)/M; yc=sum(y.*m,2)/M; vxc=sum(vx.*m,2)/M; vyc=sum(vy.*m,2)/M;
    x = x-xc;   y = y-yc;   vx = vx-vxc;    vy = vy-vyc;
end
if(useTex)
    interp = "latex";
else
    interp = "tex";
end
H = linspace(0,1-1/numbods,numbods)';   % hue
colors = hsv2rgb([H, ones(size(H)), .75*ones(size(H))]); %color order
%% Paths
figure;
colororder(colors)
plot(x, y, '-');
lgtxt = [];
for idx=1:numbods
    lgtxt = [lgtxt, "$$r_{"+idx+"}(t)$$"];
end
legend(lgtxt, Interpreter=interp); set(gca,'TickLabelInterpreter',interp)
title("Body Locations", Interpreter=interp); axis equal;
xlabel("X Position ($$x$$) [km]", Interpreter=interp); ylabel("Y Position ($$y$$) [km]", Interpreter=interp); grid on;
xl = xlim; yl=ylim;
%% Animation
figure;
colororder(colors)
objs = plot(nan,nan, 'ko'); animAx = gca;
xlim(xl); ylim(yl);
title("Body Locations", Interpreter=interp); subtitle("$$t=0$$ years", Interpreter=interp);
xlabel("X Position ($$x$$) [km]", Interpreter=interp); ylabel("Y Position ($$y$$) [km]", Interpreter=interp);
set(gca,'TickLabelInterpreter',interp)
for idx=1:numbods
    al(idx) = animatedline("LineWidth", 0.1,"Color",colors(idx,:));
end
pause; tic;

animation = VideoWriter("animation"+string(datetime, "yyyy-MM-dd-hhmmss")+".mp4",'MPEG-4');      % ANIMATIONCOMMENT
animation.FrameRate = fps; animation.Quality=100; open(animation);      % ANIMATIONCOMMENT
sp = (length(t)-1)/(fps*tEnd/tScale); % frames (sim) per frame (display)
for idx=sp:sp:length(t)
    set(objs, "XData", x(idx,:), "YData", y(idx,:));
    for alIdx=1:length(al)
        addpoints(al(alIdx), x(idx+(1-sp:0),alIdx), y(idx+(1-sp:0),alIdx));
    end
    animAx.Subtitle.String = sprintf("$$t=%.2f$$ years", years(seconds(t(idx))));
    %pause(max([1/fps-toc, 0]));tic;   % ANIMATIONCOMMENT
    writeVideo(animation, getframe(gcf));   % ANIMATIONCOMMENT
end
close(animation);   % ANIMATIONCOMMENT


%% Helper function
% it works, don't touch
function dsv = dynamics(t, sv, m, num)
G = 6.6743e-20; % division to convert to km as base distance unit
dynMtx = zeros(4*num,4*num);
for targetnum = 1:num
    dynMtx((targetnum-1)*4+[1,2],(targetnum-1)*4+[3,4]) = eye(2);
    for sourcenum = 1:num
        d = norm(sv((1:2)+(sourcenum-1)*4) - sv((1:2)+(targetnum-1)*4));
        gmd3 = (G*m(sourcenum)/d^3)*eye(2);
        if(targetnum ~= sourcenum)
            dynMtx((targetnum-1)*4+[3,4],(sourcenum-1)*4+[1,2]) = gmd3;
            dynMtx((targetnum-1)*4+[3,4],(targetnum-1)*4+[1,2]) = dynMtx((targetnum-1)*4+[3,4],(targetnum-1)*4+[1,2])-gmd3;
        end
    end
end
dsv = dynMtx * sv;
end