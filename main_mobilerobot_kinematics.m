clear all
close all
clc

%% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% test platform %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %mass
% m1 = 150;
% m2 = 45;
% m = m1 + m2;
% g = 9.81;
% 
% %wheel radius
% r = 0.1;
% %wheel base
% D = 1.2;
% %width
% L = 0.8;
% 
% %z value from global origin to link1 mass center
% h1 = 0.45;
% %z value from global origin to link2 mass center
% h2 = 1.3;
% %z value from joint 0 to link1 mass centor
% h0c1 = h1;
% %z value from joint 1 to link2 mass centor
% h1c2 = 0.7;
% 
% %x value from global origin to link2 mass center
% a = 0.4;
% %y value from global origin to link2 mass center
% b = 0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% prototype? %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mass
m1 = 50;
m2 = 40;
m = m1 + m2;
g = 9.81;

%wheel radius
r = 0.1;
%wheel base
D = 0.6;
%width
L = 0.4;

%z value from global origin to link1 mass center
h1 = 0.7;
%z value from global origin to link2 mass center
h2 = 1.5;
%z value from joint 0 to link1 mass centor
h0c1 = h1;
%z value from joint 1 to link2 mass centor
h1c2 = 0.5;

%x value from global origin to link2 mass center
a = 0.2;
%y value from global origin to link2 mass center
b = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%link1 mass center location from origin
rc1 = [0 0 h1]';
%link2 mass center location from origin
rc2 = [a b h2];
%link1 mass centor location from joint0
r0c1 = [0 0 h1]';
%link2 mass centor location from joint1
r1c2 = [a b h1c2]';

%% kinematics
dt = 0.01;

%trajectory generation
% max velocity
v_max = 1;
% first, define the width of road
road_width = 2;
%and then, find the radius of circle
circle_r = (road_width - L*cos(pi/4)) / (1-cos(pi/4));
%and then, find the time to draw circle for pi/2
t = 0 : dt : pi *circle_r/(2*v_max) +0.5;
%and then, find the angular velocity, (linear velocity is set to 1)
v = v_max*[-cos(pi .* (t(1:1/dt)))/2+0.5 1*ones(1,length(t) - 2/dt) +cos(pi .* (t(1:1/dt)))/2+0.5];
w = 1*(v_max/circle_r)*[-cos(pi .* (t(1:1/dt)))/2+0.5 1*ones(1,length(t) - 2/dt) +cos(pi .* (t(1:1/dt)))/2+0.5];

% v = v_max*[1*ones(1,length(t) - 1/dt) +cos(pi .* (t(1:1/dt)))/2+0.5];
% w = (v_max/circle_r)*[1*ones(1,length(t) - 1/dt) +cos(pi .* (t(1:1/dt)))/2+0.5];


% v = v_max*[-cos(pi .* (t(1:1/dt)))/2+0.5 1*ones(1,length(t) - 1/dt) ];
% w = (v_max/circle_r)*[-cos(pi .* (t(1:1/dt)))/2+0.5 1*ones(1,length(t) - 1/dt) ];

% v = v_max*[1*ones(1,length(t))];
% w = 5*(v_max/circle_r)*[1*ones(1,length(t))];


theta(1) = 0;
for i = 2 : length(t)
    theta(i) = theta(i-1) + w(i-1)*dt;
end

wheel2body = [r/2 r/2; -r/L r/L];
[w_wheel] = inv(wheel2body) * [v; w];
wl = w_wheel(1,:);
wr = w_wheel(2,:);

for i = 1 : length(t)
    velocity_global(:,i) = [cos(theta(i)) 0 ; sin(theta(i)) 0; 0 1] * [v(i); w(i)];
end

position_global(:,1) = [0;0;0];
for i = 2 : length(t)
    position_global(:,i) = position_global(:,i-1) + velocity_global(:,i-1)*dt;
end

acceleration_global(:,1) = [0;0;0];
for i = 2 : length(t)
    acceleration_global(:,i) = (velocity_global(:,i) - velocity_global(:,i-1))/dt;
end

%% ZMP calculation
for i = 1 : length(t)
    dHx(i) = -(m1*h1 + m2*h2)*acceleration_global(2,i) - m2*h2*(a*acceleration_global(3,i) - b*velocity_global(3,i)^2);
    dHy(i) = (m1*h1 + m2*h2)*acceleration_global(1,i) - m2*h2*(b*acceleration_global(3,i) + a*velocity_global(3,i)^2);
    
    rpx(i) = (-dHy(i) + a*m2*g) / (m*g);
    rpy(i) = (dHx(i) + b*m2*g) / (m*g);
end

%
figure;plot(t,rpx,'b')
hold on
plot(t,rpy,'r')
legend('rpx', 'rpy')
grid on
title('ZMP point')
xlabel('t(sec)'); ylabel('y(m)')
% set(gcf,'Position', [10 550 560 420])
set(gcf,'Position', [10 460 448 336])


%
for i = 1 : length(theta)
    rpx_body(:,i) = inv([cos(theta(i)) -sin(theta(i));sin(theta(i)) cos(theta(i))])*[rpx(i);rpy(i)];
end
figure;plot(rpx_body(1,:),rpx_body(2,:),'LineWidth',2)
hold on
% plot(rpx,rpy,'--r','Linewidth',2)
plot([-D/2 D/2], [-L/2 -L/2],'k')
plot([-D/2 D/2], [L/2 L/2],'k')
plot([-D/2 -D/2], [-L/2 L/2],'k')
plot([D/2 D/2], [-L/2 L/2],'k')
t_r = [0 : 0.01 : 2*pi];
plot(L/2*cos(t_r), L/2*sin(t_r))
axis equal;
grid on
title('ZMP pint w/ polygon')
xlabel('x(m)'); ylabel('y(m)')
% set(gcf,'Position', [580 550 560 420])
set(gcf,'Position', [480 460 448 336])

%
% figure;plot(t,sqrt(rpx.^2 + rpy.^2),'b')
% hold on
% plot(t,L/2*ones(length(t)),'r')


%
figure;plot(t,wl,'b')
hold on
plot(t,wr,'r')
grid on
title('wheel velocity')
xlabel('t(sec)'); ylabel('velocity(rad/sec)')
% set(gcf,'Position', [10 50 560 420])
set(gcf,'Position', [10 40 448 336])


%
figure;plot(position_global(1,:), position_global(2,:), 'LineWidth', 2)
grid on
axis equal;
title('X-Y plane')
xlabel('x(m)'); ylabel('y(m)')
% set(gcf,'Position', [580 50 560 420])
set(gcf,'Position', [480 40 448 336])
hold on
axis equal;
plot([-circle_r+road_width circle_r], [0 0],'k', 'LineWidth',3)
plot([circle_r circle_r], [0 2*circle_r-road_width],'k', 'LineWidth',3)
plot([-circle_r+road_width circle_r-road_width], [road_width road_width],'k', 'LineWidth',3)
plot([circle_r-road_width circle_r-road_width], [road_width 2*circle_r-road_width],'k', 'LineWidth',3)

plot([0 circle_r*cos(-pi/2)], [circle_r circle_r+circle_r*sin(-pi/2)],'c')
box_x = [-D/2 D/2 D/2 -D/2]';
box_y = [0 0 L L]';
box1 = fill(box_x, box_y,'r');

plot([0 circle_r*cos(-pi/4)], [circle_r circle_r+circle_r*sin(-pi/4)],'c')
rot_45 = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)];
box_45 = rot_45 * [box_x box_y]';
box2 = fill(box_45(1,:)'+circle_r*cos(-pi/4), box_45(2,:)'+circle_r+circle_r*sin(-pi/4),'r');

plot([0 circle_r*cos(0)], [circle_r circle_r+circle_r*sin(0)],'c')
rot_90 = [cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
box_90 = rot_90 * [box_x box_y]';
box3 = fill(box_90(1,:)'+circle_r*cos(0), box_90(2,:)'+circle_r+circle_r*sin(0),'r');



% figure;
% subplot(3,1,1)
% plot(t,position_global(1,:))
% subplot(3,1,2)
% plot(t,velocity_global(1,:))
% subplot(3,1,3)
% plot(t,acceleration_global(1,:))





