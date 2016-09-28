# MicroWave
# SomeMatlabSimualtion
close all;
clear all;
clc;

beta = 0*pi;
M = 360;
theta = linspace(0,pi,M/2);
phi = linspace(0,2*pi,M);
[tt,pp] = meshgrid(theta,phi);

ttM = acos(-sin(tt).*cos(pp));
ppM = angle(-sin(tt).*sin(pp)+1i*cos(tt));
%电流1
E1tm = (-exp(1i*0.25*pi*cos(ttM))./sin(ttM) - cos(ttM).*exp(-1i*0.25*pi*cos(ttM))./(1i*sin(ttM))).*exp(1i*0.25*pi*sin(ttM));
E1t = E1tm .* (-cos(ttM).*cos(ppM).*cos(tt).*sin(pp)-cos(ttM).*sin(ppM).*sin(tt)+sin(ttM).*cos(tt).*cos(pp));
E1p = E1tm .* (-cos(ttM).*cos(ppM).*cos(pp)-sin(ttM).*sin(pp));

% Ett = abs(abs(E1t)+1i*abs(E1p));
% [x1,y1,z1] = sph2cart(pp,pi/2-tt,Ett);
% mesh(x1,y1,z1);
% diff = 20*log10(max(max(Ett))/min(min(Ett)))

%电流2
E2t = 2*tan(tt).*sin(0.25*pi*cos(tt));
E2p = 0;

%电流3
E3tm = (exp(1i*0.25*pi*cos(ttM))./sin(ttM) + cos(ttM).*exp(-1i*0.25*pi*cos(ttM))./(1i*sin(ttM))).*exp(-1i*0.25*pi*sin(ttM));
E3t = E3tm .* (-cos(ttM).*cos(ppM).*cos(tt).*sin(pp)-cos(ttM).*sin(ppM).*sin(tt)+sin(ttM).*cos(tt).*cos(pp));
E3p = E3tm .* (-cos(ttM).*cos(ppM).*cos(pp)-sin(ttM).*sin(pp));

Et = E1t + E2t*exp(1i*beta) + E3t;
% Et = E1t + E2t.*exp(1i*0.25*pi*sin(tt)) + E3t;
Ep = E1p + E2p + E3p;

E13t = E1t+E3t;
E13p = E1p+E3p;
E13 = sqrt(abs(E13t).^2+abs(E13p).^2);
% E13 = 20*log10(E13);
E2 = abs(E2t);
% E2 = 20*log10(E2);
% E = sqrt(abs(Et).^2+abs(Ep).^2);
% E = 20*log10(E);

% E = abs(E13t+E13p+E2t*exp(1i*beta));
%求椭圆极化总场的振幅
alpha = angle(Et-Ep);
EtA = abs(Et);
EpA = abs(Ep);
A = 1./(-EtA.^2.*sin(alpha).^2);
C = 1./(-EpA.^2.*sin(alpha).^2);
B = 2*cos(alpha)./(EtA.*EpA.*sin(alpha).^2);
E = sqrt(-2./(A+C+sqrt((A-C).^2+B.^2)));
Etest = abs(Ep);
[i,j] = find(Etest <= 0.005);
E(i,j) = Et(i,j);
E = abs(E);

[x1,y1,z1] = sph2cart(pp,pi/2-tt,E13);
mesh(x1,y1,z1);

[x2,y2,z2] = sph2cart(pp,pi/2-tt,E2);
figure
mesh(x2,y2,z2);

[xx,yy,zz] = sph2cart(pp,pi/2-tt,E);
figure
mesh(xx,yy,zz);
diff = 20*log10(max(max(E))/min(min(E)))
