clc
% constants
L=15e-2; %length of beam
b=5.2e-2; %width of beam
h=20e-5; %height of beam
E1=1.9e11; %young modulus of Si
E2=5.2e10; %young modulus of PZT
Db=2329; %density of the Si
Dm= 19250; %density of proof mass Tungsten
Dp=7.5e3; %density of PZT
Xm=0.005; %length of proof mass
Ym=5e-4;
Zm=1.5e-2;
t=1e-5; %thickness of PZT
e33=3500*8.85e-12;
e31=250;
g=9.81;
m=Dm*Ym*Zm*Xm; %mass of proof mass
m2=Db*L*b*h+Dp*L*b*t; % mass of the beam
meff=m+0.25*m2;
a=g;%acceleration for walking used for power calculation
a2=2.6*g; %acceleration of jumping used for calculation of safety factor
F=meff*a; %effective force
Fend=m*a2; %force at the end of the beam

I=b*h^3/12;
%part1-calculation of the natural frequency

%u=F*L^3/(3*E1*I); %maximum deflection at the end of the beam
ybar=(b*h^2*E1/2+(h+t/2)*b*t*E2)/(b*h*E1+b*t*E2); %center of gravity
EIeff=E1*(((b*h^3)/12)+b*h*(ybar-h/2)^2)+E2*(((b*t^3)/12)+b*t*(t/2+h-ybar)^2); %effective EI
K=3*EIeff/L^3;
W=sqrt(K/meff); %natural frequency rad/s
fr=W/(2*pi) %Hz

%stress distribution
[x,y]=meshgrid(0:0.003:L,-h/2:0.000004:h/2);
sigdis=-y*(L-x)*F/I;
mesh(x,y,sigdis)
xlabel 'X-axis'
ylabel 'Y-axis'
zlabel 'Stress distribution'


%part2- power calculation
A=(F*L^2)/(2*EIeff); %slope at the end of the beam
R=t/(b*L*e33*W); %Resistance
C=(W*h*e31*A*b)/(2*(1+b*L*e33*W*R/t)); %Current
V=R*C %Voltage
P=((W*b*h*e31*A)^2/(8*(1+b*L*e33*W*R/t))^2)*R %power


%calculation of maximum stress
rou=Fend*L/EIeff;
sigmax=E1*h/2*rou;
sigallow=75*10^6;
safetyfactor=sigallow/sigmax

