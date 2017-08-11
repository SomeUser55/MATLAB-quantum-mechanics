% This function searches for permitted levels of Energy.
% Psi function is symmetric and tends to 0 when x tends to inf.
% System is: a harmonic oscillator. Potential energy acts like x^2.
% Shrodinger equation: (d/d(ksi))^2*psi + (Eps - ksi^2)*psi = 0.
% Eps and ksi are nondimensional energy and coordinate respectively: Eps=E/(h/(2*pi)*omega/2),
% ksi=x/L; L is scale of oscillator: L=sqrt(h/2*pi*m*omega) [L]=meter; h is a Planck constant; m is a mass
% and omega is a frequency.

function SokolovIgor_Harmonic_Oscillator_3

global L;
L=53*10^(-12);
global accuracy; %It is considered that psi=0 with this accuracy.
accuracy=1.e-6; 
global Eps;
Eps=linspace(0,10,101);
global Eps_n;
Eps_n=0;
N_Points=1000;


ksi=linspace(0, 15, N_Points);

y_start=[1; 0];
% y(1)=phi y(2)=phi_dot
array=zeros(101,2,1000);
for i=1:101
Eps_n=Eps(i);
    [Ksi, Y]=ode45(@RightSide, ksi, y_start);
  array(i,1,:)=Ksi;
array(i,2,:)=Y(:,1);
end %for i=1:N_Points

figure;
hold on;
xlabel('ksi  axis');
ylabel('psi axis');
set(gca, 'XLim', [0, 15]);
set(gca, 'YLim', [-1.5, 1.5]);
drawnow;
for i=1:9
  %plot(array(i,1,:),array(i,2,:), '-r');  
 comet(array(i,1,:),array(i,2,:));
end 

end %SokolovIgor_Harmonic_Oscillator_1

function y_prime=RightSide(ksi,y)
global Eps_n;
% y(1)=psi y(2)=psi_prime
y_prime = [y(2); -(Eps_n-ksi^2)*y(1)];
end %function y_prime=RightSide(ksi,y)

