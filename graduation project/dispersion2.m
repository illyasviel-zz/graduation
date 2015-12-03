%用傅立叶变换法进行数据处理的Matlab计算程序：
clear;
clc;
clf;
global Co
Co=299792458;   %light velocity constant

%INPUT %
widt=1400;      %sample points
win_ct=790;     %filter center
win_wd=30;       %filter half width

x1=1260;         %start wavelength
x2=1400;         %end wavelength


L=750;         %fiber length unit:m
order=3;        %Order of polynomial Fit
nFit=20;      %fit points
sign=-1;        %sign of dispersion

%Main Program %
 [filename, pathname] = uigetfile('*.txt', 'Pick an Data file');
        if isequal(filename,0)
           disp('User selected Cancel')
        else
           disp(['User selected', fullfile(pathname, filename)])
        end
FilePath=[pathname filename];
fid = fopen(FilePath,'r');
[A,count] = fscanf(fid,'%f');
t=length(A)/2;

wavelength=p;
omega1=2*pi*Co./wavelength*10^9;
signal1=q;
l1=length(signal1);

figure(1)
plot(wavelength,signal);
xlabel('\fontsize{12}\bfWavelength(nm)');
ylabel('\fontsize{12}\bfIntensity(a.u)');

figure(2)
plot(omega,signal);
xlabel('\fontsize{12}\bfomega(nm)');
ylabel('\fontsize{12}\bfIntensity(a.u)');

figure(3)
fsi=fft(signal);
fsi=fftshift(fsi);
l2=length(fsi);
plot(abs(fsi),'b');
xlabel('\fontsize{12}\bfTime');
ylabel('\fontsize{12}\bfIntensity(a.u)');
grid on 

figure(4)
win=zeros(1,widt);
win(round(win_ct-win_wd):round(win_ct+win_wd))=1;
w_fsi=fsi.*win;           %truncate the correlation part
plot(abs(w_fsi));
xlabel('\fontsize{12}\bfTime');
ylabel('\fontsize{12}\bfIntensity(a.u)');

figure(5)
w_fsi=ifftshift(w_fsi);
fsi=ifft(w_fsi);
wphase=-phase(fsi);
wphase=wphase-wphase(round(widt/2));
b=1;
for a=1:50
    subphase(a)=wphase(b);
    subomega(a)=omega(b);
    b=b+19;
end
plot(subomega,subphase,'o');
xlabel('\fontsize{12}\bfOmega(rad/s)');
ylabel('\fontsize{12}\bfPhase(rad)');

figure(9)
set(gca,'FontSize',10);
hold on;
plot(lambda_Fit,beta2*(10^24),'b');

xlabel('wavelength\lambda(nm)','FontSize',12);
ylabel('beta2(ps^2/m)]','FontSize',12);
box on
%grid on
hold off;

%SAVE %       %save to file
Dispersion=[lambda_Fit;D];
savefile=strcat('D_','.TXT');
save(savefile,'Dispersion','-ascii');
Beta2=[lambda_Fit;beta2*(10^24)];
savefile=strcat('beta2_','.TXT');
save(savefile,'Beta2','-ascii');

