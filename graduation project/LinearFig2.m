FWHM=500*10^-3;%pico-second
tao0=1.2011*FWHM;%full 1/e maximum
lamda0=1558;%nm
deltalamda=2;%nm
l=1*10^3;%um
Ch=2.0;%nm/cm
DeltaL0=0;%mm
ne=2.0;
delta=0;
%still need:
c=3*10^5;%nm/ps
Omega1=0;%ps
%Omega2=440;%ps^2


Omega2_lamda=538;%ps/nm
Omega2=Omega2_lamda*lamda0^2/c;
Omega3=0;
%Omega3=Ch*Omega2^3/tao0;%ps^3


%Omega2_lamda=Omega2*c./lamda^2

z=30.8;%km

detune=-1:0.01:1;
lamda=lamda0+400.*detune;%nm
lamda1=lamda0+2.*detune;%nm

eta=0.1;deltan=0.15;N=1;

t1=-1:0.001:1;%ns
t=1000.*t1;%ps

%(a) Envelope of the generated pulse. 
z=z/40;%normalize Z
A=4*(t-Omega1*z)./tao0;
B=32*Omega3*z./tao0^3;
C=8*Omega2*z./tao0^2;
P=A.^2-1/(1+C^2).*B*(1-C^2-((1-3*C^2)*A.^2)./(6*(1+C^2))).*A-1/(2*(1+C^2)^3)*B^2*(1-6*C^2+C^4-(1-10*C^2+5*C^4)/(8*(1+C^2)).*A.^2).*A.^2;
r=(1+C^2)^-0.25*exp(-0.25/(1+C^2).*P);
normalized_r=r/max(r);
%
GP=lamda0/(2*ne);
gamma=1/(eta*deltan).*(lamda./lamda0-1);
phase_term=eta*deltan.*(1000*N*GP./lamda).*(1-gamma.^2).^0.5;
w=abs((sinh(phase_term).^2)./((cosh(phase_term)).^2-gamma.^2));
normalized_w=w/max(w);
w1=(1+cos(c*tao0./(lamda.*ne)));

%w=exp(-(tao0.*pi*c./(ne.*lamda)).^2).*(1+cos(c*tao0./(lamda.*ne)));

%T=0.5.*w.*(1+cos(4*pi*ne/lamda0^2.*lamda.*(DeltaL0+delta.*lamda/Ch)));
%T=0.5*1.*(exp(-(t./(tao0/2)).^2)+exp(-((t-tao0/2)./(tao0/2)).^2));
T1=0.5.*(1+cos(4*pi*ne/lamda0^2.*lamda1.*10^6*(DeltaL0+10*deltalamda/Ch)));%exp(-(lamda.*Omega2./(tao0/2)).^2);
normalized_T1=T1/max(T1);
T=w.*1;%T1;
%T=exp(-(lamda.*Omega2./(tao0/2)).^2).*(1+cos(4*pi*ne/lamda0^2.*lamda.*(DeltaL0+delta.*lamda/Ch)));
%T=normalized_w.*normalized_T1;
normalized_T=T./max(T);
%y=0.5.*r.*w.*(1+cos((4*pi*ne.*t/(lamda0^2*Omega2_lamda)).*(DeltaL0+delta.*t/(Ch*Omega2_lamda))));
%normalized_y=y/max(y);

TBWP=4*ne*Ch*l^2/lamda0^2;%Time-Bandwidth Product
%plot(lamda1,w);
%plot(lamda1,normalized_w);
plot(lamda1,normalized_T);
%plot(t1,normalized_y);
%fprintf('     The Time-Bandwidth Product is %8.5f\n', TBWP)