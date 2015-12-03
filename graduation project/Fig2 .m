
f0=29;%Grad/s
t0=80.9;%pico-second
FWHM=350*10^-3;%pico-second
tao0=1.2011*FWHM;%full 1/e maximum
Omega1=0;%ps
Omega2=440;%ps^2
Omega3=14;%ps^3
z=20;%km
t2=-5:0.001:5;%ns
%t2=-0.4:0.001:0.4;%ns
t=-5000:1:5000;%ps

%(a) Envelope of the generated pulse. 
z=z/40;%normalize Z
A=4*(t-Omega1*z)./tao0;
B=32*Omega3*z./tao0^3;
C=8*Omega2*z./tao0^2;
P=A.^2-1/(1+C^2).*B*(1-C^2-((1-3*C^2)*A.^2)./(6*(1+C^2))).*A-1/(2*(1+C^2)^3)*B^2*(1-6*C^2+C^4-(1-10*C^2+5*C^4)/(8*(1+C^2)).*A.^2).*A.^2;
r=(1+C^2)^-0.25*exp(-0.25/(1+C^2).*P);
normalized_r=r/max(r);
%(b) Spectrum (inset: zoom-in display).
f=0:0.1:50;%GHz
t1=t0^2./(8*pi^2*Omega3.*(f./1000).^2)-Omega2^2/(2*Omega3);
A1=4*(t1-Omega1*z)./tao0;
P=A1.^2-1/(1+C^2).*B*(1-C^2-((1-3*C^2)*A1.^2)./(6*(1+C^2))).*A1-1/(2*(1+C^2)^3)*B^2*(1-6*C^2+C^4-(1-10*C^2+5*C^4)/(8*(1+C^2)).*A1.^2).*A1.^2;
r1=(1+C^2)^-0.25*exp(-0.25/(1+C^2).*P);
normalized_r1=r1*9;
%(c) Instantaneous frequency versus time (solid line: predicted by (6), circle: ob-tained from numerical result). 
wrf=(1000*(Omega2^2+2*Omega3.*t).^-0.5*t0)/(2*pi);%GHz
%wrflinear=(1000*(1/Omega2-Omega3/(Omega2)^3.*t3)*t0)/(2*pi);%GHz
%(d) Compressed pulse obtained by autocorrelation(inset: zoom-in display).
TBWP=Omega2*(f0/1000)*(1/tao0-2*pi*(f0/1000));%Time-Bandwidth Product
%i=r.*(1+cos((Omega2^2+2*Omega3.*t).^0.5*t0/Omega3+(Omega2*t0/Omega3)));
r2=0.25*exp(-2.*(t2./(tao0/2)).^2)+0.25*exp(-2.*((t2-t0)./(tao0/2)).^2)+0.5*exp((-t2.^2-(t2-t0).^2)./((tao0/2))^2).*cos(t0/Omega2.*t2);
r2=r2./max(r2);
i=r2.*(1+cos((Omega2^2+2*Omega3.*(t2*1000)).^0.5*t0/Omega3+(Omega2*t0/Omega3)));
normalized_i=i/max(i);

%(e) zoom-in display of Spectrum

%(f) zoom-in display of compressed pulse

%i=r2.*(1+cos((Omega2^2+2*Omega3.*t).^0.5*t0/Omega3+(Omega2*t0/Omega3)));


%i=r.*r;

%normalized_i=i/max(i);

subplot(321);
plot(t,normalized_r);
%plot(1000*t2,normalized_i);
title('Envelope of the generated pulse');
xlabel('time(ps)');
ylabel('Normalized Amplitude');

subplot(322);
plot(f,normalized_r1);
ylim([0,1]);
title('Spectrum');
xlabel('frequency(GHz)');
ylabel('Normalized Amplitude');

subplot(323);
plot(t,wrf);%t,wrflinear);
title('Instantaneous frequency versus time');
xlabel('time(ps)');
ylabel('frequency(GHz)');

subplot(324);
plot(1000*t2,normalized_i);
title('Compressed pulse obtained by autocorrelation');
%plot(t,i);
xlabel('time(ps)');
ylabel('Normalized Amplitude');

subplot(325);
plot(f,normalized_r1);
xlim([20,45]);
ylim([0,0.1]);
title('Instantaneous frequency versus time (zoom-in display)');
xlabel('frequency(GHz)');
ylabel('Normalized Amplitude');

subplot(326);
plot(1000*t2,normalized_i);
xlim([-400 400]);
title('Compressed pulse obtained by autocorrelation (zoom-in display)');
xlabel('time(ps)');
ylabel('Normalized Amplitude');

fprintf('     The Time-Bandwidth Product is %8.5f\n', TBWP)


