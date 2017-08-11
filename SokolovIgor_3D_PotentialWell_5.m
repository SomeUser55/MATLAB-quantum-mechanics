% Three dimensional potentional well. Potential V(r) is isotropic. V(r)=0 r<a, V(r) = V_0 r>a.
% Hamiltonian H=(p^2)/2m + V(r)
% Stationary Shroedinger`s eq.: H*psi(r)=E*psi(r)
% x=r*sin(eta)*cos(phi), y=r*sin(teta)*sin(phi), z=r*cos(eta)
% Angular momentum operator L=(L_x, L_y, L_x) L=[r p]
% l - eigenvalue of L
% [L_x, H]=[L_y, H]=[L_z, H]=0 => [L^2, H]=0
% Radial momentum p_r=(h_bar/i)*(1/r)*(d/dr)*r=(h_bar/i)*(d/dr + 1/r)
% p^2=p_r^2 + L^2/r^2, r>0.
% H=(p_r^2)/2m + L^2/(2*m*r^2) + V(r), (p_r^2)/2m - radial kinetic energy,
% L^2/(2*m*r^2) - rotational kinetic energy, V(r) - potential energy.
% [(p_r^2)/2m + L^2/(2*m*r^2) + V(r)] * psi(r,teta,phi) = E * psi(r,teta,phi)
% At first we find eigenvalues of L^2, then we find eigenfunctions of L^2, which
% obeys Shroed. eq.
%
% Spherical functions
% [L^2, L_x]=[L^2, L_y]=[L^2, L_z]=0.
% Spherical functions Y_l_m(teta,phi) -  eigenfunctions of L^2 and L_x
% (L^2) * Y_l_m(teta,phi) = l(l+1)*h_bar^2 * Y_l_m(teta,phi)
% L_z * Y_l_m(teta,phi)= m*h_bar * Y_l_m(teta,phi), l=0,1,2...; m=-l,...,+l
% 
% Radial equation
% We find eigenfunction  of H, L^2, L_z:
% psi_l_m(r,teta,phi)=Y_l_m(teta,phi) * chi_l(r)
% psi_l_m(r,teta,phi) is eigenfunction of H, Y_l_m(teta,phi) eigenfunction of L^2 
% => [(p_r^2)/2m + L^2/(2*m*r^2) + V(r) - E] = chi_l(r) = 0
% y_l(r)=r * chi_l(r)
% [(-h_bar^2/(2*m))*(d/dr)^2 + l(l+1)h_bar^2/(2*m*r^2) + V(r) - E]*y_l(r)=0
% solution: y_l(r)~ r^(l+1) r->0 => chi_l(r)~ r^l r->0
%
% Thus and so 3-D Shroed. eq. has reduced itself to 1-D one (28). 
% Effective potential U(r)= l(l+1)h_bar^2/(2*m*r^2) + V(r) 
% We consider a diffuse state so l=2
% If E<U an energy spectrum is discrete. 
% If E>U an energy spectrum is continuum.

function SokolovIgor_3D_PotentialWell_5
clear all;
global angstrom U_0 U_0_SI a a_SI m h_bar e_charge_SI sol k_wave r_min L
angstrom=10^(-10);  %meters
e_charge_SI = 1.602176565*10^(-19);
U_0=30; %eV
U_0_SI=U_0*1.602176*10^(-19); % Joule
a=5; %angstroms
a_SI=a*10^(-10); %meters
m=9.109382*(10^(-31)); %kg
h_bar = 1.054571*(10^(-34)); %Joule*sec
L=2; % *h_bar; L is an angular momentum

r_max=25.; %angstorms
r_min=0.5*10^(-3); %1; %3; %0.5; %1; %0.5*10^(-3); %angstorms
k_wave = pi/2/r_max;
r_arr = linspace(r_min, r_max, 300);

Guess_Energy=0.5;%  1; %eV
RelErr = 1.e-6;  % -6   -9
AbsErr = 1.e-7; % -7  -10
options = bvpset('RelTol', RelErr, 'AbsTol', AbsErr, 'Vectorized', 'on');

solinit = bvpinit(r_arr, @psi_guess, Guess_Energy);
sol = bvp4c(@RightSide, @Boundary_Cond, solinit, options);
Energy=sol.parameters;

norm=sqrt(2*integral(@chi_squared, 0 , r_max)); % norm of radial wave function
r_arr2=[linspace(0, r_min, 10) r_arr];
chi_arr=chi_function(r_arr2, norm);

h=figure;
hold on; 
xlabel('radius, angstroms');
ylabel('chi axis');
plot(r_arr2, chi_arr, '.-g', 'LineWidth', 3);
set(gca, 'XLim', [0, 10]);
set(h, 'Name', 'Bound diffuse state: radial wave function');

fprintf('a = %6.4f angstrom, U = %6.4f eV\n', a, U_0);
fprintf('Energy of d-state = %6.6f eV\n\n', Energy);

%k=figure;
%hold on; 
%set(gca, 'XLim', [0, 10]);
%set(gca, 'YLim', [0, 40]);
%set(k, 'Name', 'Potential well');
%label('radius, angstroms');
%ylabel('U_eff axis');
%L=0;
%lot(r_arr2, U_eff(r_arr2), '.-r', 'LineWidth', 1);
%=1;
%plot(r_arr2, U_eff(r_arr2), '.-g', 'LineWidth', 1);
%L=2;
%plot(r_arr2, U_eff(r_arr2), '.-b', 'LineWidth', 1);
%legend('L=0', 'L=1', 'L=2');
plot(r_arr2, psi_guess(r_arr2), '.-b', 'LineWidth', 1);

end %function SokolovIgor_3D_PotentialWell
%================================================================================================

function y = psi_guess(r)
global k_wave 

y=[sin(2*k_wave*r); cos(k_wave*r)];

end % function psi_guess
%==============================================================================================

function res = Boundary_Cond(ya, yb,Energy)
global r_min L;

res=[ya(1) - r_min^(L+1); ya(2) - (L+1)*r_min^L ; yb(1)]; 

end % function Boundary_Cond
%======================================================================================================

function y_prime = RightSide(r, y, Energy)
global m h_bar angstrom e_charge_SI 

y_prime = [y(2,:);  -(2*m*(angstrom^2)/(h_bar^2)) * e_charge_SI * ( Energy - U_eff(r)).*y(1,:) ];
end %function y_prime = RightSide(x, y, Energy)
%======================================================================================================

function res=U(r) % Returns an original potential well: U(r)=0 r<a, U(r) = U_0 r>a
global a U_0;

res(numel(r))=0.;

switch int2str([isempty(r(r<a)) isempty(r(r>=a))])
    case '0  1'
        res(r<a)=0.;
    case '1  0'
        res(r>=a)=U_0;
    case '0  0'
        res(r<a)=0.;
        res(r>=a)=U_0;
end %switch
end %function res=U(x)
%======================================================================================================

function res=U_eff(r) % Returns effective pontial well: U_eff(r) = U(r) + h_bar^2*L*(L + 1)/(2mr^2); r_min<r<r_max.
global m h_bar e_charge_SI angstrom L;

res=U(r) + L*(L+1)*(h_bar)^2./(2*m*(r*angstrom).^2)/e_charge_SI;
end %function res=U_eff(r)
%======================================================================================================


function res=chi_squared(r) % Returns unnormalized squared radial wave function within the interval 0<r<r_max.
global sol r_min L

switch int2str([isempty(r(r<r_min)) isempty(r(r>=r_min))])
    case '0  1'
        res=(r(r<r_min)).^(L);
    case '1  0'
        res=deval(sol, r(r>=r_min), 1)./r(r>=r_min);
    case '0  0'
        res=[(r(r<r_min)).^L  deval(sol, r(r>=r_min), 1)./r(r>=r_min)];
end %switch

res=res.^2;
end %function f=chi_squared(r)
%======================================================================================================

function res=chi_function(r, norm) % Returns normalized radial wave function within the interval 0<r<r_max
global sol r_min L

switch int2str([isempty(r(r<r_min)) isempty(r(r>=r_min))])
    case '0  1'
        res=(r(r<r_min)).^(L);
    case '1  0'
        res=deval(sol, r(r>=r_min), 1)./r(r>=r_min);
    case '0  0'
        res=[(r(r<r_min)).^L  deval(sol, r(r>=r_min), 1)./r(r>=r_min)];
end %switch

res=res./norm;
end %function res=chi_function(r)
%======================================================================================================



















