function linear_nonchirped()
%Photonic Generation
j=(-1)^0.5;
Omega1=0;
Omega2=80;%ps^2
Omega3=0;
z=65;%km
tao0=210*10^-3;%ps
deltat=9.88;%ps
t=-1000:0.004:1000;%ns
Psi=0;
tao=abs(Omega2)/tao0;%ps
fc=deltat/(2*pi*abs(Omega2))*1000;%GHz
TBWP=abs(Omega2)*fc/1000*(1/tao0-2*pi*fc/1000);

%£¨a£©envelope of the generated pulse
z=z/40;%normalize Z
A=4*(t-Omega1*z)./tao0;
B=32*Omega3*z./tao0^3;
C=8*Omega2*z./tao0^2;
P=A.^2-1/(1+C^2).*B*(1-C^2-((1-3*C^2)*A.^2)./(6*(1+C^2))).*A-1/(2*(1+C^2)^3)*B^2*(1-6*C^2+C^4-(1-10*C^2+5*C^4)/(8*(1+C^2)).*A.^2).*A.^2;r=(1+C^2)^-0.25*exp(-0.25/(1+C^2).*P);
normalized_r=r/max(r);
%
Ed=exp(t.^2/(tao));%.*exp(-j*t.^2/(2*Omega2));
Eout=(exp(-(t-deltat).^2/tao^2).*exp(-j*(t-deltat).^2/(2*Omega2))+exp(-t.^2/tao^2).*exp(-j*t.^2/(2*Omega2)).*exp(j*Psi));

I=r.*abs(Eout).^2;
subplot(2,2,1);

plot(t/1000,I);
xlabel('a.time');
ylabel('current')
linear_wav1()
linear_wav2()
linear_wav3()
