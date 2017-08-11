% Our system is a one-dimensional rectangular potential energy well.
% Parameters "a > 0" and "U_0 > 0" are a half-width and a depth of the well respectively. [a]=angstrem, [U_0]=eV 
% Box:
%   U = U_0     |x|>a,
%   U = 0       |x|<a. 
% Schrodinger equations: 
%   (d/dx)^2 * psi(x) + 2*m/h_bar^2 * E * psi(x) = 0            |x|<a,
%   (d/dx)^2 * psi(x) + 2*m/h_bar^2 * (E - U_0) * psi(x) = 0    |x|>a.
% Substitution:
%   k^2 = 2*m/h_bar^2 * E              |x|<a
%   alpha^2 = (2*m/h_bar) * (U_0 - E)  |x|>a
% Schrodinger equations:
%   (d/dx)^2 * psi(x) = -k^2 * psi(x)       |x|<a
%   (d/dx)^2 * psi(x) = alpha^2 * psi(x)    |x|>a 
% Even solutions:
%   psi_in=cos(kx)                            |x|<a,
%   psi_out=B*exp(-alpha*x) B=exp(eta)*cos(z)  x>a.
%   z = a*k, eta = alpha*a
% System: 
%   eta = sqrt(C^2 - z^2) 
%   eta = z*tan(z),
%   C=sqrt( 2*m*U_0_SI*a_SI^2/(h_bar)^2 )
%   eta = z*tan(z) where z = a*k, eta = alpha*a, h_bar = h/2pi, [h_bar]=Joule*sec;
%   sqrt(const-z^2)=f1; z*tan(z)=f2; 
% Odd solutions:
%   psi_in=sin(kx)  |x|<a
%   psi_out=B*exp(-alpha*x) B=exp(eta)*cos(z)  x>a.
% System:
%   eta = sqrt(C^2 - z^2) 
%   eta = - z*cot(z)
%   sqrt(const-z^2)=f1; -z*cot(z)=f2; 
% root(1,i) is an array of z components; root(2,i) is an array of eta components. 

function SokolovIgor_PotentialWell_16
clear all;
global eps U_0 U_0_SI a a_SI m h_bar C e_charge_SI root ifEven1_ifOdd2
e_charge_SI = 1.602176565*10^(-19);
eps=2.2204*(10^(-16));
U_0=30; %eV
U_0_SI=U_0*1.602176*10^(-19); % Joule
a=5; %angstroms
a_SI=a*10^(-10); %meters
m=9.109382*(10^(-31)); %kg
h_bar = 1.054571*(10^(-34)); %Joule*sec
C=sqrt(2*m*U_0_SI*a_SI^2/(h_bar)^2);
fprintf('a = %6.4f angstrom, U = %6.4f eV, C = %6.4f\n', a, U_0, C);

for ifEven1_ifOdd2=1:2
z=linspace(0, C, 3000);
f1=sqrt(C^2-z.^2);
f2=Fun_f2(z);
NumberOfRoots=Fun_NumberOfRoots();
intervals=Fun_intervals(NumberOfRoots);
root=Fun_root(NumberOfRoots, intervals); 
Energy=(root(1,:)/a_SI).^2 *(h_bar^2)/(2*m)/e_charge_SI;
Norm_arr=Fun_norm();
B=Fun_B; 
P_in=2*integral(@(x)psi_squared(x), 0, a, 'ArrayValued', true)./(Norm_arr.^2);
x_squared_avrg=2*integral(@Fun_x_squared_times_psi_squared, 0, inf, 'ArrayValued', true)./(Norm_arr.^2);
PlotAndPrint(z, f1, f2, P_in, Norm_arr, Energy, B, NumberOfRoots, x_squared_avrg);
end %for ifEven1_ifOdd2=1:2

end %function SokolovIgor_Potential_Well
%==================================================================================================

function PlotAndPrint(z, f1, f2, P_in, Norm_arr, Energy, B, NumberOfRoots, x_squared_avrg)
global ifEven1_ifOdd2 C root a 
Odd_or_Even={'Even levels:', 'Odd levels:'};
Odd_or_Even_legend={'z*tan(z)', '-z*cot(z)'};
h(ifEven1_ifOdd2)=subplot(2,2,ifEven1_ifOdd2);
axes(h(ifEven1_ifOdd2));
hold on;
plot(z,f1,'.-r', z,f2,'.-b');
title(Odd_or_Even{ifEven1_ifOdd2});
set(gca,'Ylim',[-pi/2, C+pi/2], 'Xlim',[-pi/2, C+pi/2]);
xlabel('eta axis');
ylabel('eta axis');
legend('sqrt(C^2 - z^2)', Odd_or_Even_legend{ifEven1_ifOdd2}, 'Location', 'NorthOutside');
hold off;

fprintf('%s\n', Odd_or_Even{ifEven1_ifOdd2});
fprintf('   ½tate   Energy,eV    z=k*a        P_in         norm^2       B, 10^5      x^2_avrg, A^2\n');
fprintf('   %1.0f. %12.5f %12.5f %12.5f %12.5f %12.5f %12.4f\n', [(1:numel(Energy)); Energy; (root(1, :)); P_in'; (Norm_arr.^2)'; B/(10^5); x_squared_avrg']);

x=linspace(0, 2*a, 2000);
psi=Wave_Function(x, Norm_arr);
h(ifEven1_ifOdd2+2)=subplot(2,2,ifEven1_ifOdd2+2);
axes(h(ifEven1_ifOdd2+2));
hold on;
plot3((kron(ones(NumberOfRoots, 1), x))', (kron((linspace(1, NumberOfRoots, NumberOfRoots))', ones(1, numel(x))))', psi', '.-g', 'LineWidth', 3);
xlabel('x axis');
zlabel('psi axis');
ylabel('State axis');
grid on;
view(64,52);
hold off;
end %function PlotAndPrint(z, f1, f2, P_in, Norm_arr, Energy, B, NumberOfRoots)
%==================================================================================================

function res=Fun_x_squared_times_psi_squared(x) % This function returns the product of x^2 times psi^2. Psi is unnormalized wave function. 
res=psi_squared(x).*x^2;
end %function res=Fun_x_squared_avrg()
%==================================================================================================


function root=Fun_root(NumberOfRoots, intervals)
global C ifEven1_ifOdd2
root=zeros(2, NumberOfRoots);
if ifEven1_ifOdd2==1
    j=0;
    for i=1:2:NumberOfRoots*2-1
        j=j+1;
        a_left = intervals(i);
        b_right = intervals(i + 1) - 1.e-10;
        root(1,j)=fzero(@diff, [a_left, b_right]);
    end %for i=1:2:NumberOfRoots*2-1
    for i=1:NumberOfRoots
        root(2,i)=sqrt(C^2-(root(1,i))^2);
    end %for i=1:NumberOfRoots
else 
    j=0;
    for i=2:2:NumberOfRoots*2
        j=j+1;
        a_left = intervals(i);
        b_right = intervals(i + 1) - 1.e-10;
        root(1,j)=fzero(@diff, [a_left, b_right]);
    end %for i=1:2:NumberOfRoots*2
    for i=1:NumberOfRoots
    root(2,i)=sqrt(C^2-(root(1,i))^2);
    end %for i=1:NumberOfRoots
end %if else
end %function root=Fun_root(NumberOfRoots)
%==========================================================================

function B=Fun_B()
global root ifEven1_ifOdd2
if ifEven1_ifOdd2==1
    B=exp(root(2,:)).*cos(root(1,:));
else 
    B=exp(root(2,:)).*sin(root(1,:));
end %if else
end %function B=Fun_B()
%==========================================================================

function f2=Fun_f2(z)
global ifEven1_ifOdd2
if ifEven1_ifOdd2==1
    f2=z.*tan(z);
else 
    f2=-z.*cot(z);
end %if else
end %function f2=Fun_f2(z)
%===================================================================================================

function NumberOfRoots=Fun_NumberOfRoots()
global ifEven1_ifOdd2 C
if ifEven1_ifOdd2==1
    NumberOfRoots=floor(C/pi)+1;
else 
    NumberOfRoots=floor(C/pi - 0.5) + 1;
end %if else
end %function Fun_NumberOfRoots()
%===================================================================================================

function intervals=Fun_intervals(NumberOfRoots)
global ifEven1_ifOdd2 C
if ifEven1_ifOdd2==1
    intervals=zeros(1,NumberOfRoots*2);
    for i=1:(NumberOfRoots*2)-1
        intervals(i+1)=intervals(i)+pi/2;
    end %for i=1:(NumberOfRoots*2)-1
    intervals(NumberOfRoots*2)=min(intervals(NumberOfRoots*2), C );
else 
    intervals=zeros(1,NumberOfRoots*2 + 1);
    for i=1:(NumberOfRoots*2)
        intervals(i+1)=intervals(i)+pi/2;
    end %for i=1:(NumberOfRoots*2)-1
    intervals(NumberOfRoots*2 + 1)=min(intervals(NumberOfRoots*2 + 1), C );
end %if else
end %function res=Fun_intervals()
%===================================================================================================


function d=diff(z) %This function computes the difference between f1 and f2
global ifEven1_ifOdd2 C;
if ifEven1_ifOdd2==1
    d=(z*tan(z))-(sqrt(C^2-z^2));
else 
    d=(-z*cot(z))-(sqrt(C^2-z^2));
end %if else
end %function d=diff(ksi)
%====================================================================================================

% This function calculates an integral N=(psi^2)dx   -inf<x<+inf
function Norm_arr=Fun_norm()
Norm_arr=sqrt(2*integral(@(x)psi_squared(x), 0, inf, 'ArrayValued', true) );
end %function N=norm(root,  NumberOfRoots)
%========================================================================================================

function res=psi_squared(x) %z=k*a; This function returns unnormalized squared wave functions within the interval 0<x<inf.
global ifEven1_ifOdd2 a root
indices_in=find(x<a);
indices_out=find(x>=a);

if indices_in
    if ifEven1_ifOdd2==1
        var_psi_in_squared=(cos(kron((root(1,:))', x(indices_in))/a)).^2;
    else
        var_psi_in_squared=(sin(kron((root(1,:))', x(indices_in))/a)).^2;
    end %if else
end %if indices_in

if indices_out
    B_matrix=kron((Fun_B)', ones(1, numel(indices_out) ) );
    var_psi_out_squared=(B_matrix.* (exp(-kron((root(2,:))', x(indices_out))/a)) ).^2;
end %if indices_out

if indices_in & indices_out
    res=[var_psi_in_squared var_psi_out_squared];
elseif indices_in
    res=var_psi_in_squared;
elseif indices_out
    res=var_psi_out_squared;
end %if indices_in & indices_out
end %function f=psi_in_squared(x)
%======================================================================================================

function res=Wave_Function(x, Norm_arr) %z=k*a; This function returns normalized wave functions within the interval 0<x<+inf.
global ifEven1_ifOdd2 root a ;
indices_in=find(x<a);
indices_out=find(x>=a);

if indices_in
    Norm_Matrix=kron(Norm_arr, ones(1, numel(indices_in) ) );
    if ifEven1_ifOdd2==1
        psi_in=cos(kron((root(1, 1:numel(Norm_arr)))', x(indices_in))/a);
    else
        psi_in=sin(kron((root(1, 1:numel(Norm_arr)))', x(indices_in))/a);
    end %if else
    psi_in=psi_in./Norm_Matrix;
end %if indices_in

if indices_out
    Norm_Matrix=kron(Norm_arr, ones(1, numel(indices_in) ) );
    B_matrix=kron((Fun_B)', ones(1, numel(indices_out) ) );
    psi_out=B_matrix.* (exp(-kron((root(2, 1:numel(Norm_arr)))', x(indices_out))/a));
    psi_out=psi_out./Norm_Matrix;
end

if indices_in & indices_out
    res=[psi_in psi_out];
elseif indices_in
    res=psi_in;
elseif indices_out
    res=psi_out;
end %if indices_in & indices_out
end %function psi=Wave_Function(x, Norm_arr)
%=======================================================================================================