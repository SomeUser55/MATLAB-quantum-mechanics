% This programm solves an equation: (d/dx)^2*psi(x) + (E/E_0  - U(x)/E_0)*psi(x)=0  
% System is a symmetric potential well, determined from arrays: U_y_arr, U_x_arr.
% We search for two first states. 


function SokolovIgor_FindEnergy_5

clear all

global U_x_arr U_y_arr x_arr k_wave L m h_bar E_0 sol pp i
L=10^(-10);  %meters
m=9.109382*(10^(-31)); %kg
h_bar= 1.054571726*10^(-34); %J*sec  6.58211928*10^(-16); %eV*sec
%E_0=h_D^2/(2*m*L^2);
(2*m*(L^2)/(h_bar^2)) * (1.602176565*10^(-19)) 
k_wave = pi/2/25.000;
 U_y_arr=[-20.000 -21.363 -23.757 -24.476 -22.572 -15.935  -4.825  -0.734  -0.098  -0.012]; %eV
 U_x_arr=[0.000   0.556   1.111   1.667   2.222   2.778   3.333   3.889   4.444   5.000]; % Angstrom

 pp=interp1(U_x_arr, U_y_arr, 'spline','pp');

x_arr = linspace(0, 25, 300);

%figure;
%hold on;
%set(gca, 'YLim', [-25, 1]);
plot(x_arr, U(x_arr));

Guess_Energy_i = [-24.000, -22.000]; %-20.; % -100/E_0; %-E_0; %-24.476/E_0; %-20.000;  %U(0.); % 3./2.;
RelErr = 1.e-6;  % -6   -9
AbsErr = 1.e-7; % -7  -10
options = bvpset('RelTol', RelErr, 'AbsTol', AbsErr, 'Vectorized', 'on');
Energy_i(2)=0;
%psi_arr(2, numel(x_arr))=0;
psi_arr_i(2, numel(x_arr))=0;
norm_i(2)=0;
Y=ones(1, length(x_arr));

figure;
hold on;
xlabel('x axis');
zlabel('psi axis');
ylabel('State axis');
grid on;

for i=1:2
Guess_Energy=Guess_Energy_i(i);
solinit = bvpinit(x_arr, @psi_guess, Guess_Energy);
sol = bvp4c(@RightSide, @Boundary_Cond, solinit, options);
Energy_i(i)=sol.parameters;
psi_arr=deval(sol, x_arr);
psi_arr_i(i, :)=psi_arr(1, :);
norm_i(i)=sqrt(2*quadl(@Psi_Func, 0, 25));
psi_arr_i(i, :)=psi_arr_i(i, :)./norm_i(i);
fprintf('Energy of %1.0f state = %12.6f eV\n', i-1, Energy_i(i));
Z=psi_arr_i(i, :);
plot3(x_arr, (i-1)*Y, Z, '.-g', 'LineWidth', 3);
end %for i=0:1

set(gca, 'XLim', [0, 5]);
view(27,30);
end %SokolovIgor_FindEnergy_1
%===============================================================================

function y = psi_guess(x)
global k_wave i

if i==1
    y = [cos(k_wave*x); - k_wave*sin(k_wave*x)];
else 
    y=[sin(2*k_wave*x); cos(k_wave*x)];
end %if else

end % function psi_guess
%=================================================================

function res = Boundary_Cond(ya, yb, Energy)
global i

if i==1
   res=[ya(1) - 1.0; ya(2); yb(1)];
else
   res=[ya(1); ya(2) - 1.0; yb(1)]; 
end %if else

end % function Boundary_Cond
%=================================================================

function y_prime = RightSide(x, y, Energy)
global m h_bar L
y_prime = [y(2,:); -(2*m*(L^2)/(h_bar^2)) * (1.602176565*10^(-19)) * ( Energy - U(x)).*y(1,:) ];
end
%==========================================
         
function u_arr=U(x)
global pp 
n=0;

for i=1:numel(x)
    if x(i)<=5.000
    n=n+1;
    else
        break;
    end %if else
end %for i=1:numel(x

u_arr(1:n)=ppval(pp, x(1:n));
u_arr(n+1:numel(x))=0.;
end 
%======================================


function psi=Psi_Func(x)
global sol 
psi=deval(sol, x);
psi=psi(1, :).^2;
end
%========================================================