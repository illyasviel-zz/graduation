function non_linear_chirped()
f0=29;%Grad/s
t0=80.9;%pico-second
FWHM=350*10^-3;%pico-second
tao0=1.2011*FWHM;%full 1/e maximum
Omega1=0;%ps
Omega2=450;%ps^2
Omega3=-140;%ps^3
z=20;%km
t2=-50:0.001:50;%ns
%t2=-0.4:0.001:0.4;%ns
t=-5000:1:5000;%ps
wrf=(1000*(Omega2^2+2*Omega3.*t).^-0.5*t0)/(2*pi);
r2=0.25*exp(-2.*(t2./(tao0/2)).^2)+0.25*exp(-2.*((t2-t0)./(tao0/2)).^2)+0.5*exp((-t2.^2-(t2-t0).^2)./((tao0/2))^2).*cos(t0/Omega2.*t2);
r2=r2./max(r2);
i=r2.*(1+cos((Omega2^2+2*Omega3.*(t2*1000)).^0.5*t0/Omega3+(Omega2*t0/Omega3)));
%subplot(2,2,1);
%plot(t,wrf);
%title('Instantaneous frequency versus time');
%xlabel('time(ps)');
%ylabel('frequency(GHz)');
plot(1000*t2,i);
%title('Compressed pulse obtained by autocorrelation');
%plot(t,i);
xlabel('time(ps)');
ylabel('Normalized Amplitude');
xlim([-0.04e4,0.04e4])
%non_linear_chirped1()
%non_linear_chirped2()
%non_linear_chirped3()
%Omega3=
