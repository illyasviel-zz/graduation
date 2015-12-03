%Photonic Generation
j=(-1)^0.5;
Omega2=87.4;%ps^2
tao0=210*10^-3;%ps
deltat=9.88;%ps
Psi=0;
tao=abs(Omega2)/tao0;%ps
fc=deltat/(2*pi*abs(Omega2))*1000;%GHz
TBWP=abs(Omega2)*fc/1000*(1/tao0-2*pi*fc/1000);
t=-1:0.004:1;%ns
t=1000.*t;
Ed=exp(t.^2/(tao)).*exp(-j*t.^2/(2*Omega2));
Eout=(exp(-(t-deltat).^2/tao^2).*exp(-j*(t-deltat).^2/(2*Omega2))+exp(-t.^2/tao^2).*exp(-j*t.^2/(2*Omega2)).*exp(j*Psi));
I=abs(Eout).^2;
plot(t/1000,I);
