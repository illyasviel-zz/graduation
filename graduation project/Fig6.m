
f0=24;%Grad/s
t0=67.5;%pico-second
FWHM=350*10^-3;%pico-second
tao0=1.2011*FWHM;%full 1/e maximum
Omega1=0;%ps
Omega2=446;%ps^2
Omega3=0.8;%ps^3

t=-5000:1:5000;%ps

%(a) Envelope of the generated pulse. 
z=0.5;%km
A=4*(t-Omega1*z)./tao0;
B=32*Omega3*z./tao0^3;
C=8*Omega2*z./tao0^2;
P=A.^2-1/(1+C^2).*B*(1-C^2-((1-3*C^2)*A.^2)./(6*(1+C^2))).*A-1/(2*(1+C^2)^3)*B^2*(1-6*C^2+C^4-(1-10*C^2+5*C^4)/(8*(1+C^2)).*A.^2).*A.^2;
r=(1+C^2)^-0.25*exp(-0.25/(1+C^2).*P);

normalized_r=r/max(r);
%(b) Spectrum (inset: zoom-in display).
f=0:0.01:30;%GHz
t1=t0^2./(8*pi^2*Omega3.*(f./1000).^2)-Omega2^2/(2*Omega3);
A1=4*(t1-Omega1*z)./tao0;
P=A1.^2-1/(1+C^2).*B*(1-C^2-((1-3*C^2)*A1.^2)./(6*(1+C^2))).*A1-1/(2*(1+C^2)^3)*B^2*(1-6*C^2+C^4-(1-10*C^2+5*C^4)/(8*(1+C^2)).*A1.^2).*A1.^2;
r1=(1+C^2)^-0.25*exp(-0.25/(1+C^2).*P);
normalized_r1=r1*2.5*10;
%(c) Instantaneous frequency versus time (solid line: predicted by (6), circle: ob-tained from numerical result). 
wrf=(1000*(Omega2^2+2*Omega3.*t).^-0.5*t0)/(2*pi);%GHz
%wrf=(1000*(1/Omega2-Omega3/(Omega2)^3.*t)*t0)/(2*pi);%GHz
%(d) Compressed pulse obtained by autocorrelation(inset: zoom-in display).
TBMP=Omega2*f0*(1/tao0-2*pi*f0);


i=r.*(1+cos((Omega2^2+2*Omega3.*t).^0.5*t0/Omega3+(Omega2*t0/Omega3)));
%i=r.*r;

normalized_i=i/max(i);

subplot(321);

%plot(t,normalized_r);
plot(t,normalized_i);

subplot(322);
plot(f,normalized_r1);
ylim([0,1]);

subplot(323);

plot(t,wrf);
ylim([23 25]);
%t,wrflinear);
%subplot(324);
%plot(t,normalized_i);
%plot(t,i);

subplot(325);
plot(f,normalized_r1);
xlim([22,26]);
ylim([0,0.3]);

%subplot(326);


