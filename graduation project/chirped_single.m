
function chirped_single()


% sta=1544;
% en=1554;
sta=1554;
en=1563;
v=1;



wav_len=sta:0.01:en;
wav_len=wav_len.*(10^(-9));
n_eff=1.485;

% FBG_periods=521.4*10^(-9);
FBG_periods=524.42*1e-9;

FBG_wavlen=FBG_periods*(2*n_eff)
% FBG_wavlen_1=FBG_periods_1*(2*n_eff);



pp=192*FBG_periods*10;
L=pp*10;%光栅长度
N=500;
aa=L/N;%曝光区长度
z=[-N/2:N/2-1]/N*L;

Cg=2e-1;
% Cg=2e-3;
a=2;
b=1;




FBG_periods_z=FBG_periods*(1-Cg*z);
fai=2*pi*(z(a)/FBG_periods_z(a)-z(b)/FBG_periods_z(b))/2

beita=0.8;
alfa=2;

% %有切趾的时候
%dn=5*10^(-4).*(cos(pi*z/L)).^2;%exp(-10*z.^2/L^2);%(1+tanh(beita*(1-2*(2*z/L).^alfa)));%(1-(2*z/L).^2)./(1-(2*0.3*z/L).^2);%ones(1,N);%.*exp(-10*z1.^2/L^2);%(cos(pi*z1/L)).^2;%.*exp(-10*z1.^2/L^2);%
% 
% 无切趾的时候
dn=5*10^(-4).*ones(1,N);


plot(z,dn)
dn_neff=5*10^(-4);
k=pi*v*dn(1)./wav_len;
p=2*pi*n_eff./wav_len*aa;
% dp=p-pi/FBG_periods+2*pi*dn_neff./wav_len;

dfai=(4*pi*n_eff)*Cg*z*FBG_wavlen./(2*n_eff*FBG_periods_z).^2;
dp=2*pi*n_eff./wav_len-pi/FBG_periods+2*pi*dn_neff./wav_len-dfai(1);

s=sqrt(k.^2-dp.^2);
sa=s*aa;


% dn=5*10^(-4);
% k=pi*v*dn./wav_len;
% p=2*pi*n_eff./wav_len;
% dp=p-pi/FBG_periods+2*pi*dn./wav_len;
% s=sqrt(k.^2-dp.^2);
% sa=s*aa;


% A=(cosh(sa)-i*dp.*sinh(sa)./s)*exp(-i*fai).*exp(i*p*(aa));
% B=(-i*k.*sinh(sa)./s)*exp(-i*fai).*exp(-i*p*(aa));
% C=(i*k.*sinh(sa)./s)*exp(i*fai).*exp(i*p*(aa));
% D=(cosh(sa)+i*dp.*sinh(sa)./s)*exp(i*fai).*exp(-i*p*(aa));
A=(cosh(sa)-i*dp.*sinh(sa)./s).*exp(-i*p*(aa));
B=(-i*k.*sinh(sa)./s).*exp(i*p*(aa));
C=(i*k.*sinh(sa)./s).*exp(-i*p*(aa));
D=(cosh(sa)+i*dp.*sinh(sa)./s).*exp(i*p*(aa));

a1=A;
b1=B;
c1=C;
d1=D;

for x=2:N

%       k=pi*v*dn(x)./wav_len;
%       dp=p-pi/FBG_periods+2*pi*dn_neff./wav_len;
%       s=sqrt(k.^2-dp.^2);
%       sa=s*aa;    
% %     
% fai=2*pi*(z(x)/FBG_periods_z(x)-z(x-1)/FBG_periods_z(x-1))/2;

% A=(cosh(sa)-i*dp.*sinh(sa)./s)*exp(-i*fai).*exp(i*p*(aa));
% B=(-i*k.*sinh(sa)./s)*exp(-i*fai).*exp(-i*p*(aa));
% C=(i*k.*sinh(sa)./s)*exp(i*fai).*exp(i*p*(aa));
% D=(cosh(sa)+i*dp.*sinh(sa)./s)*exp(i*fai).*exp(-i*p*(aa));

k=pi*v*dn(x)./wav_len;
dp=2*pi*n_eff./wav_len-pi/FBG_periods+2*pi*dn_neff./wav_len-dfai(x);
s=sqrt(k.^2-dp.^2);
sa=s*aa;

 A=(cosh(sa)-i*dp.*sinh(sa)./s).*exp(-i*p*(aa));
B=(-i*k.*sinh(sa)./s).*exp(i*p*(aa));
C=(i*k.*sinh(sa)./s).*exp(-i*p*(aa));
D=(cosh(sa)+i*dp.*sinh(sa)./s).*exp(i*p*(aa));
    
    
    A1=a1.*A+b1.*C;
    B1=a1.*B+b1.*D;
    C1=c1.*A+d1.*C;
    D1=c1.*B+d1.*D;
    
    a1=A1;
    b1=B1;
    c1=C1;
    d1=D1;
end

r=abs(C1./A1);
R=r.^2;

% 
Q=phase(-C1./A1);
   Y(1)=Q(1);Y(2)=Q(2);Y(3)=Q(3);
   for ii=4:length(wav_len)
       if(abs(Q(ii-1)-Q(ii))<=1)
           Y(ii)=(sta+ii*0.001)^2*1e-18/(2*pi*3e-4)*((Q(ii-1)-Q(ii))/(0.001e-9));
       else
           Y(ii)=(sta+ii*0.001)^2*1e-18/(2*pi*3e-4)*((Q(ii-3)-Q(ii-2))/(0.001e-9));
       end
   end  


% hold on
% R=10*log10(R);
wav_len=wav_len*1e9;
subplot(2,2,1)
plot(wav_len,R);
axis([1550 1565 0 1]);
xlabel('(a）Wavelength（nm）');
ylabel('Reflectivity');
hold on
chirped_single1();
chirped_single2();
chirped_single3();





%figure(2)
%plot(wav_len,Y)
%xlabel('Wavelength（nm）');
%ylabel('Group delay (ps)');
   

% lambda0=FBG_wavlen%nm
% %c=2.0e-9/(1e-2)%nm/cm
% delta_lambda=wav_len-lambda0;
% delta_L0=0;
% c=Cg;
% lambda0=FBG_wavlen;
% T_lambda=1/2*R.*(1+cos(4*pi*n_eff.*wav_len/(lambda0^2).*(delta_L0+abs(delta_lambda)/c)))
% %figure
% plot(wav_len,T_lambda)
% 




