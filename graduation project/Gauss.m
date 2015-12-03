Seq=[1];

Ts=0.4e-9;         
num=2^8;
dt=Ts/num;        
fs=1/dt;
nt=Ts/2;
m=1;
CC=0;
FWHM=0.13e-9;
T0=FWHM/(2*sqrt(log(2)));
Number=length(Seq)*num;

t=zeros(1,Number);
Es=zeros(1,Number);

for k=1:1:Number
    t(k)=dt*(k-1);
    Es(k)=Seq(fix(t(k)/Ts)+1)*exp(-(1+i*CC)/2*((t(k)-(fix(t(k)/Ts)*Ts)-nt)/T0)^(2*m));
end
figure 
plot(t*1e9,Es/max(Es),'k')
