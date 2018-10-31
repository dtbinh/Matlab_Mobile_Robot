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
% %moment of inertia
% I = 1/12*m1*(D^2 + L^2) + m2*(a^2+b^2); %m2는 mulped mass로 가정
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
b = 0.0;
%moment of inertia
I = 1/12*m1*(D^2 + L^2) + m2*(a^2+b^2); %m2는 mulped mass로 가정
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% torque input to force input
B = [1/r 1/r; -L/2/r L/2/r];

%each wheel velocity to body velocity
C = [r/2 r/2; -r/L r/L];

%% Dynamics
dt = 0.01;

v_max = 1;
road_width = 2;
circle_r = (road_width - L*cos(pi/4)) / (1-cos(pi/4));

t = 0 : dt : pi *circle_r/(2*v_max) +0.5;

v_ref = v_max*[-cos(pi .* (t(1:1/dt)))/2+0.5 1*ones(1,length(t) - 2/dt) +cos(pi .* (t(1:1/dt)))/2+0.5];
w_ref = 1*(v_max/circle_r)*[-cos(pi .* (t(1:1/dt)))/2+0.5 1*ones(1,length(t) - 2/dt) +cos(pi .* (t(1:1/dt)))/2+0.5];

s_ref(1) = 0;
t_ref(1) = 0;
for i = 1 : length(t)-1
    s_ref(i+1) = s_ref(i) + v_ref(i) * dt;
    t_ref(i+1) = t_ref(i) + w_ref(i) * dt;
end

s_local(1) = 0;t_local(1) = 0;
v_local(1) = 0;w_local(1) = 0;
dv_local(1) = 0;dw_local(1) = 0;
TauL(1) = 0;TauR(1) = 0;
for i = 1 : length(t)-1
    
    % computed torque control
    ddq_ctc(:,i) = -50*( [s_local(i); t_local(i)] - [s_ref(i); t_ref(i)] ) -20*( [v_local(i); w_local(i)] - [v_ref(i); w_ref(i)] );
%     ddq_ctc(:,i) = -50*( [s_local(i); t_local(i)] - [s_ref(i); t_ref(i)] ) -2*( [v_local(i); w_local(i)] - [v_ref(i); w_ref(i)] );
    TAU(:,i) = inv(B)*([m 0;0 I]*ddq_ctc(:,i) + [-m*a*w_local(i)^2; m*a*w_local(i)*v_local(i)]);
    TauL(1,i) = TAU(1,i);
    TauR(1,i) = TAU(2,i);
    
    % Forward Dynamics
    dv_local(i+1) = (1/r*(TauR(i) + TauL(i)) + m*a*w_local(i)^2) / m;
    dw_local(i+1) = (L/2/r*(TauR(i)-TauL(i)) - m*a*w_local(i)*v_local(i)) / I;
     
    v_local(i+1) = v_local(i) + dv_local(i) * dt;
    w_local(i+1) = w_local(i) + dw_local(i) * dt;
    
    s_local(i+1) = s_local(i) + v_local(i) * dt;
    t_local(i+1) = t_local(i) + w_local(i) * dt;
    
end

wheel_velocity = inv(C) * [v_local;w_local];
wl = wheel_velocity(1,:);
wr = wheel_velocity(2,:);

%% ZMP calculation
for i = 1 : length(t)
    velocity_global(:,i) = [cos(t_local(i)) 0; sin(t_local(i)) 0; 0 1] * [v_local(i);w_local(i)];
end

position_global(:,1) = [0;0;0];
for i = 2 : length(t)
    position_global(:,i) = position_global(:,i-1) + velocity_global(:,i-1)*dt;
end

acceleration_global(:,1) = [0;0;0];
for i = 2 : length(t)
    acceleration_global(:,i) = (velocity_global(:,i) - velocity_global(:,i-1))/dt;
end


for i = 1 : length(t)
    dHx(i) = -(m1*h1 + m2*h2)*acceleration_global(2,i) - m2*h2*(a*acceleration_global(3,i) - b*velocity_global(3,i)^2);
    dHy(i) = (m1*h1 + m2*h2)*acceleration_global(1,i) - m2*h2*(b*acceleration_global(3,i) + a*velocity_global(3,i)^2);
    
    rpx(i) = (-dHy(i) + a*m2*g) / (m*g);
    rpy(i) = (dHx(i) + b*m2*g) / (m*g);
end

%% Plot
%ZMP point
figure;plot(t,rpx,'b')
hold on
plot(t,rpy,'r')
legend('rpx', 'rpy')
grid on
title('ZMP point')
xlabel('t(sec)'); ylabel('y(m)')
set(gcf,'Position', [10 460 448 336])

%ZMP pint w/ polygon
for i = 1 : length(s_local)
    rpx_body(:,i) = inv([cos(s_local(i)) -sin(s_local(i));sin(s_local(i)) cos(s_local(i))])*[rpx(i);rpy(i)];
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
set(gcf,'Position', [480 460 448 336])

%wheel velocity
figure;plot(t,wl,'b')
hold on
plot(t,wr,'r')
grid on
title('wheel velocity')
xlabel('t(sec)'); ylabel('velocity(rad/sec)')
set(gcf,'Position', [10 40 448 336])

%X-Y plane
figure;plot(position_global(1,:), position_global(2,:), 'LineWidth', 2)
grid on
axis equal;
title('X-Y plane')
xlabel('x(m)'); ylabel('y(m)')
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


%local velocity
figure;
subplot(2,1,1)
plot(t,v_ref,'r','LineWidth',2)
hold on
plot(t,v_local,'--b','LineWidth',2)
grid on
title('v local')
xlabel('t(sec)'); ylabel('velocity(m/sec)')
legend('reference','controlled','Location','South')
subplot(2,1,2)
plot(t,w_ref,'r','LineWidth',2)
hold on
plot(t,w_local,'--b','LineWidth',2)
grid on
title('w local')
xlabel('t(sec)'); ylabel('velocity(rad/sec)')
legend('reference','controlled','Location','South')