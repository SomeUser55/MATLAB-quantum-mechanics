% Physical system is a simple pendulum. Initial speed of bob = 0; 
% Differential equation: (d/dt)phi_dot + omega_0^2*sin(phi) = 0, "phi" - an angle of deviation,
% omega_0^2 = g/L, "g" is a gravitational acceleration, "L" is a length of the pendulum; 
% This function plots a dependence of period from an initial angle of deviation.

function SokolovIgor_SimplePendulum_6
global omega_0;
omega_0=1;

Period_Ideal = 2*pi/omega_0;
N_Points=1000;
phi_0=linspace(0, pi, N_Points);
phi_dot_0=0.;
t_max=1000*Period_Ideal;
Time_Events=zeros(1, N_Points);
t_span=linspace(0, t_max, N_Points);

% y(1)=phi y(2)=phi_dot
options=odeset('Events', @Period);
for i=2:N_Points
    y_start=[phi_0(i); phi_dot_0];
    [t, Y, TE]=ode45(@RightSide, t_span, y_start, options);
Time_Events(i)=TE;
end %for i=1:N_Points

figure;
hold on;
xlabel('phi start  axis');
ylabel('period axis');
plot(phi_0, Time_Events, '-r');
Log_Time_Events(numel(Time_Events))=0;
Log_Time_Events(:)=log10(Time_Events(:));
figure;
hold on;
xlabel('phi start  axis');
ylabel(' log10 period axis');
set(gca,'Xlim',[3.138, pi]);
plot(phi_0, Log_Time_Events, '-r');
end %SokolovIgor_SimplePendulum_1

function y_prime=RightSide(t,y)
global omega_0;
% y(1)=phi y(2)=phi_dot
y_prime = [y(2); -(omega_0^2)*sin(y(1))];
end %function y_prime=RightSide(t,y)

function [value, isterminal, direction]=Period(t, y)
% y(1)=phi y(2)=phi_dot
%value=zeros(2,1);
value = y(1);
isterminal = 1;
direction = 0;
end %function [value, isterminal, direction]=Period(t, y)