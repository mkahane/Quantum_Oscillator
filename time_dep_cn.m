clc
clear all
close all

nx = 200;
n = 5;
a = -20;
b = -a;
dx = (b-a)/nx;
x = linspace(a,b,nx+1);


V = x.^2;

tfinal = 100;
dt = 0.1;
nt = tfinal/dt;



psi = (1/sqrt((2^n)*factorial(n))).*hermite(n,x)...
    .*(1/pi)^(.25).*exp(-.5*x.^2);
psi = psi';



H1 = zeros(nx+1,nx+1);
H1(1,1) = -2;

for i = 2:nx+1
    H1(i,i) = -2;
    H1(i-1,i) = 1;
    H1(i,i-1) = 1;
end


H2 = zeros(nx+1,nx+1);
for i = 1:nx+1
    H2(i,i) = V(i);
end


H = (1/(dx^2)).*(H1)-H2;
%H = -H;


num = eye(nx+1) + 1i.*H.*(dt/2);
denom = eye(nx+1) - 1i.*H.*(dt/2);

for t = 1:nt
    subplot(1,2,1)
    plot(x,real(psi));
    axis([a b -1 1]);
    title(strcat('real part | n = ',num2str(n), '| t = ' , num2str(t)));
    subplot(1,2,2)
    plot(x,imag(psi),'r');
    axis([a b -1 1]);
    title(strcat('imag part | n = ',num2str(n), '| t = ' , num2str(t)));
    psi = denom\(num*psi);
    pause(.05)
end

