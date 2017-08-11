% n is an index of the n-th row of Hermite_Matrix and correspomds to the (n-th)-1 condition. n=1, 2, ...
% N is an index of the N-th condition. N=0, 1, 2, ...
% Function "H=Hermite_Matrix()" calculates coefficiens of 0-th, 1-st, ..., N-th Hermite polynomes. 
% The N-th polynome H_N is defined by the (N+1)-th row of (n)x(n) matrix,
% so that H_N=H(N+1,1)*x^0 + H(N+1, 2)*x^1 + ... + H(N+1, N+1)*x^N.

% Wave function psi_n=1/sqrt(2^n * n! * pi^0.5 * L) * exp(-0.5*ksi^2) * H_2(ksi),
% where ksi is a nondimensional quantity: ksi=x/L; L is a scale of
% oscillator: L=sqrt(h/2*pi*m*omega) [L]=meter; h is a Planck constant; m is a mass
% and omega is a frequency.
function SokolovIgor_Hermite_4()
global L;
global n; % Don`t touch it! 
global N; % Touch it! 
global ksi; % ksi is an array of arguments.
L=53*10^(-12);
N=10;
n=N+1;
ksi=linspace(-10,10,10000);

Wave_Function();

end %function SokolovIgor_Hermite_1
%============================================================================================================

function H=Hermite_Matrix()
global n;
H=zeros(n, n);
H(1,1)=1;
H(2,2)=2;

for k=3:n
    for m=1:k
        H(k,m)=-2*(k-2)*H(k-2,m);
    end %m=1:k
    
    for m=2:k
        H(k,m)=H(k,m) + 2*H(k-1,m-1);
    end %m=2:k

end %for k=3:n+1

end %function H=Hermite_Matrix(n)
%========================================================================================================

function Herm_Values=Hermite_Values()
global n ksi;
H=zeros(n, n);
H=Hermite_Matrix();
Herm_Values=zeros(n,numel(ksi));
for i=1:n
    
for k=1:i
    Herm_Values(i, :)=Herm_Values(i, :)+H(i,k).*ksi.^(k-1);
end %k=1:i
    
end %i=1:n

end %function Hermite_Values=Hermite(ksi, H_n)
%========================================================================================================


function psi=Wave_Function()
global L n ksi;
psi=zeros(n,numel(ksi));
Herm_Values=zeros(n,numel(ksi));
Herm_Values=Hermite_Values();
 
  for i=1:n
psi(i, :)= exp(-0.5.*ksi.^2) .*Herm_Values(i, :) .* (1/sqrt(2^i * factorial(i) * pi^0.5 * L));
  end %i=1:n
  
figure;
hold on;
xlabel('ksi axis');
zlabel('psi axis');
ylabel('State axis');
grid on;
Y=ones(1, length(ksi));
for i=1:n
    Z=psi(i,:);
    plot3(ksi, (i-1)*Y, Z, '.-g', 'LineWidth', 3);
end %i=1:n

view(72,54);
end %function psi=Wave_Function(ksi, H_n)
%========================================================================================================