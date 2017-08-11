% System is a quantum oscillator in a ground state with U=x^2 + alpha*x^4, where 0<=alpha<=1.
% This program solves 100 shroedinger equations: 
% psi^2 + (2*Eps - ksi^2 - alpha*ksi^4 =0  alpha=0:0.01:1
% where all values are nondimensional: Eps=E/(h_prime*omega),
% ksi=x/sqrt(h_prime/m*omega)
 
function SokolovIgor_QuantumOscillator_bvp_6()

global k_wave sol alpha_n  alpha psi_n norm_n;
alpha_n=linspace(0,1,100);
psi_n=zeros(100,2,200);
norm_n=zeros(1,100);
psi=zeros(1,200);
%y(1) = psi
%y(2) = psi_prime

Fig_Psi_h = figure;
set(Fig_Psi_h, 'Position', [10 60 1120 720]); %@1200x800
set(Fig_Psi_h, 'Name', 'Ground state wavefunction', 'NumberTitle', 'off');
Axes_Psi_h = axes('Parent', Fig_Psi_h, 'Color', [1.0 1.0 1.0], ...
             'Units', 'points', 'Position', [50 80 600 400], ...
             'Fontsize', 10);

ksi_max = 15.0;
N_Points = 200;
ksi_span = linspace(0, ksi_max, N_Points);
k_wave = pi/2/ksi_max;

plot(Axes_Psi_h, [0, ksi_max], [0, 0], '-k');
xlabel('ksi'); 
ylabel('psi');
set(gca, 'XLim', [0, 3.7]);
set(gca, 'YLim', [-0.1, 1.4]);
hold on;

Guess_Energy = 3./2.; % E/(hbar*omega);
%solinit = bvpinit(x_span, @psi_guess, Guess_Energy);

RelErr = 1.e-6;  % -6   -9
AbsErr = 1.e-7; % -7  -10
options = bvpset('RelTol', RelErr, 'AbsTol', AbsErr, 'Vectorized', 'on');
Energy(100)=0;

for i=1:100
    
    alpha = alpha_n(i);
    if i==1 solinit = bvpinit(ksi_span, @psi_guess, Guess_Energy);
    else solinit=sol;
    end %if
 
sol = bvp4c(@RightSide, @Boundary_Cond, solinit, options);
Energy(i)=sol.parameters;

psi_n(i,:,:)=deval(sol, ksi_span);
norm_n(i)=quadl(@psi_2, 0, ksi_max);

psi=squeeze(psi_n(i,1,:));
plot( ksi_span, psi(:)/sqrt(norm_n(i)), 'g');
hold on;
drawnow;
end %for i=1:100


Energy_FirstOrderCorrection(100)=0;
Energy_SecondOrderCorrection(100)=0;
Energy_FirstOrderCorrection(:)=0.5 + (3/8)*alpha_n(:);
Energy_SecondOrderCorrection(:)=0.5 + (3/8)*alpha_n(:) -(21/32) * alpha_n(:).^2;

Fig_Energy_h=figure; 
set(Fig_Energy_h, 'Name', 'Ground state energy');
hold on;
xlabel('alpha axis');
ylabel('Energy axis, h bar*omega');

set(gca, 'YLim', [0, 1]);
plot(alpha_n, Energy, 'r', 'LineWidth', 2);
plot(alpha_n, Energy_FirstOrderCorrection, 'b', 'LineWidth', 2);
plot(alpha_n, Energy_SecondOrderCorrection, 'g', 'LineWidth', 2);
legend('Exact solution', 'First order correction', 'Second order correction', 'Location', 'North');
%Energy = sol.parameters;
    

%fprintf('i = %4d    Energy = %12.6f\n', i, Energy)

%fprintf('The lowest energy level is %7.3f\n', sol.parameters)

end % Schroedinger_Harmonic_bvp
%=================================================================

function y = psi_guess(ksi)
global k_wave
y = [cos(k_wave*ksi); - k_wave*sin(k_wave*ksi)];
end % function psi_guess
%=================================================================

function y = psi_2(ksi)

global sol

psi = deval(sol, ksi);
y = psi(1, :).^2;

end % function psi_2
%=================================================================

function res = Boundary_Cond(ya, yb, Energy)

res = [ya(1) - 1.0;    ya(2);       yb(1)];

end % function Boundary_Cond
%=================================================================

function y_prime = RightSide(ksi, y, Energy)
global alpha
%y(1) = psi
%y(2) = psi_prime

y_prime = [y(2, :);
          - (2*Energy - ksi.^2 - alpha*ksi.^4).*y(1, :)];

end %  function y_prime = RightSide(x, y)
