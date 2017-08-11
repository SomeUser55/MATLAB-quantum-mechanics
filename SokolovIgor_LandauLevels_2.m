% System is a quantum oscillator with U=x^2. Search for n-th level.
% This program solves shroedinge equation: 
% psi^2 + (2*Eps - ksi^2)*psi =0 
% where all values are nondimensional: Eps=E/(h_prime*omega),
% ksi=x/sqrt(h_prime/m*omega)
 
function SokolovIgor_LandauLevels_2()

global k_wave sol n p_x y_0;

n=5; %level. n=1,3,5,7,9...
%y(1) = psi
%y(2) = psi_prime

c = 2.99792458 * 10 ^ 10; %cm/s
e = -4.803204673 * 10 ^ (-10); %Franklin
m = 9.10938291 * 10 ^ (-28); %g
H = 0.5; 
omega_H = abs(e) * H / m / c;
h_bar = 1.054571726 * 10 ^ (-27); %erg*second
sigma = 0.5;
p_x = m * c/100; %g * cm / s     1E-3;
y_0_CGS = - c * p_x / e / H;
y_0=y_0_CGS/sqrt(h_bar*omega_H); % non-dimensional y_0

ksi_max = 15.0;
N_Points = 4000;
ksi_span = linspace(0, ksi_max, N_Points);
k_wave = pi/ksi_max;

xlabel('ksi'); 
ylabel('psi');

Guess_Energy =n +0.5; %3./2.; % E/(hbar*omega);
plot(ksi_span, psi_guess(ksi_span)); % guess function
hold on;
solinit = bvpinit(ksi_span, @psi_guess, Guess_Energy);
 
sol = bvp4c(@RightSide, @Boundary_Cond, solinit);
Energy=sol.parameters;
fprintf('Energy = %12.6f\n', Energy)

psi=deval(sol, ksi_span, 1);
norm=quadl(@psi_2, 0, ksi_max);
plot(ksi_span, psi/sqrt(norm), 'g');

end % Schroedinger_Harmonic_bvp
%=================================================================

function y = psi_guess(ksi)
global k_wave n
y = [sin(k_wave*ksi); k_wave*cos(k_wave*ksi)];
end % function psi_guess
%=================================================================

function y = psi_2(ksi)

global sol

psi = deval(sol, ksi);
y = psi(1, :).^2;

end % function psi_2
%=================================================================

function res = Boundary_Cond(ya, yb, Energy)

res = [ya(1);    ya(2)-1;       yb(1)];

end % function Boundary_Cond
%=================================================================

function y_prime = RightSide(ksi, y, Energy)
global y_0
%y(1) = psi
%y(2) = psi_prime

y_prime = [y(2, :);
          - (2*Energy - (ksi-y_0).^2).*y(1, :)];

end %  function y_prime = RightSide(x, y)
