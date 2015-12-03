clc;
clear all;
dt=1e-12;
Num=100000;
t=(1:1:Num)*dt;
T=500*dt;
tau=100*dt;
P=exp(-((t-T)/tau).^2);                   %��˹����ʱ����ʽ

df=1/dt/Num;
f=(1:1:Num)*df;
Pf=fft(P);
Pf=fftshift(P);

plot(f,Pf);xlim([5e11,5.1e11])


