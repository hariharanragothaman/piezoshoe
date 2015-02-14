clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%  MATERIAL PROPERTIES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Setting Values for PZT Materail 
hpzt=10e-6;
rho_PZT=7600;
wpzt=5e-3;
%w=0.1*10^-3:1*10^-6:20*10^-3;
E=5.2e10;
lpzt=14.7e-3;
%l=1e-3:1e-7:20e-3;
%l=15e-3;
mpzt=hpzt*rho_PZT.*wpzt.*lpzt;

% Parameters of Auln Bonding Layer
lauln=lpzt;
hauln=5e-6;
wauln=wpzt;
%rho_AULN
%Eauln=;
%mauln=rho_AULN.*lauln.*hauln.*wauln;

%Si Substrate layer
lsi = lpzt;
hsi = 10e-6;
wsi=wpzt;
rho_SI=2.33e3;
Esi=1.7e11;
msi=rho_SI.*lsi.*hsi.*wsi;

%Proof Mass
lm=6e-3;
hm=4.5e-4;
wm=wpzt;
rho_MASS=1.93e4;
mmass=rho_MASS.*lm.*wm.*hm;

m_total = msi+msi+mmass;
%m_total = mmass;

h_total = hpzt+hsi+hauln;
%h_total = 4.5e-4;



w = wpzt;
l = lpzt-(lm/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% K CALCULATION
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%k2 calculation - for walking frequency
f=10; 
wn=2*pi*f;
k2=(wn*wn).*m_total

%k1 calculation - for beam 
I=w.*(h_total.^3)/12
k1=(3*Esi*I)./(l.^3)
%/(l.^3);

%PLOTTING GRAPH 
figure;
grid on;
%xlim([h(8800) h(end)]); 
plot(wpzt,k1,'r');
hold on;
grid on;
%xlim([h(8800) h(end)]); 
plot(wpzt,k2,'b');
%hold on;
%grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%Power calculations
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k31 = 3.5e-1 %memsnet
S = 0.03 %assumed as 0.01 in paper
Y=E % youngs modulus of PZT
k=1-k31^2%1-k31^2
Et=1700 %dielectric constant
Eo = 8.85*10^-12;
e = Et*Eo
Apzt = wpzt.*lpzt
C = e.*Apzt/hpzt;
w_resonance=wn
Ropt = (1./(wn.*C)).*((2.*S)./sqrt((4.*(S.^2))+(k31.^4)))
d=1.8*10^-10%memsnet
g=9.8
A = 2.5
b=S.*(2.*sqrt(m_total.*k1)) % S * (2*sqrt(meff*ksp))

P = (((Ropt.*C.^2).*((Y.*d.*g)./(k.*Et)).^2)/((2.*S.*(w_resonance).^2.*(Ropt.*C)).^2+(2.*S.*w_resonance+(1-k).*Ropt.*C.*w_resonance.^2).^2)).*(3*b/(2.*lpzt.^2)).*(A).^2


figure;
grid on;
plot(wpzt,P,'g');
