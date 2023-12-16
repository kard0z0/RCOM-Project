clear al; close all; clc;


Pe=10;                                  % Emisson Power
he=40;                                  % Height of emitter
hr=100;                                 % Height of receptor
f=10*10^9;                              % Frequence in Hertz
fmhz=f*10^-6                            % Frequence in MegaHertz
d=45000;                                % Distance in meters
dkm=d*10^-3;                            % Distance in kilometers
comp_onda=3*10^8/f;                     % Wavelength
l=3;                                    % Diameter of the antenna's dish
rend=0.5;                               % Efficiency
Ag=44.3;                                % Atenuation dos guias por km
Dg=60;                                  % Length dos guias em m
r = 6370*10^3;                          % Earth's ray

Pre=1017;                               %Atmospheric Pressure
Temp=18;                                %Temperature in Celsius Degrees
no=1000;
dndh=-43*10^-6;

Ar=-(10*log10(10^-((Ag*Dg/1000)/10)))               
Ae=-Ar

A0=32.4+20*log10(dkm)+20*log10(fmhz)

ang_f = atan((hr-he)/d)

Gr=20*log10((pi*l)/comp_onda)+10*log10(rend) 
Ge=Gr

Pd=Pe-Ae-Ar+Ge+Gr-A0
Pdw=10^(Pd/10)

x=0:45000;
y=-((x-(d/2)).^2)/(2*r);

a=1/2*sqrt((d-0)^2+(hr-he)^2);
b=sqrt(comp_onda*a)/2;

y0=-((0-(d/2)).^2)/(2*r)+he;
y1=-((d-(d/2)).^2)/(2*r)+hr;

t=linspace(0,2*pi,300);
X=a*cos(t);
Y=b*sin(t);
ang=atan2((y1-y0),d);

m=(0+d)/2+X*cos(ang)-Y*sin(ang);
n=(y0+y1)/2+X*sin(ang)+Y*cos(ang);

tanE = d*x/(2*r) + he - d^2/(8*r);

tanR = -d*x/(2*r) + (3/8)* d^2/r + hr;

syms a;
S=double(real(solve(a.^3-1.5*d*a.^2+(0.5.*d^2-r*(he+hr))*a+he*d*r==0,a,'MaxDegree',3)))

for i=1:length(S)
    if S(i)>0 && S(i)<d
        P_especular=S(i)
    end
end

d1=P_especular;
d2=d-d1;

he_equi=he-d1^2/(2*r)
hr_equi=hr-d2^2/(2*r)

y01=-((0-(d/2)).^2)/(2*r)+he_equi;
y11=-((d-(d/2)).^2)/(2*r)+hr_equi;

ang_incidencia=atan(he_equi/d1);

D=1/(sqrt(1+(2*d1*d2)/(r*d*sin(ang_incidencia))))

%figuras
figure()
plot(x,y);
xlim([-5*10^3 50*10^3]);
hold on
plot(tanE);
plot(tanR);
plot([0 0],[-((0-(d/2)).^2)/(2*r) -((0-(d/2)).^2)/(2*r)+he]);
plot([d d],[-((d-(d/2)).^2)/(2*r) -((d-(d/2)).^2)/(2*r)+hr]);       
plot([0 d],[-((0-(d/2)).^2)/(2*r)+he -((d-(d/2)).^2)/(2*r)+hr]);
plot(m,n,'-k');
plot(P_especular,-((P_especular-(d/2)).^2)/(2*r),'.k','MarkerSize',15);
grid on;
hold off

a2=1/2*sqrt((d-0)^2+(hr_equi-he_equi)^2);
b2=sqrt(comp_onda*a2)/2;

t1=linspace(0,2*pi,300);
X1=a2*cos(t1);
Y1=b2*sin(t1);
ang_i=atan((((y11)-(y01)))/d)

m_equi=(d)/2+X*cos(ang_i)-Y1*sin(ang_i);
n_equi=(he_equi+hr_equi)/2+X1*sin(ang_i)+Y1*cos(ang_i);

figure()
xlim([-5*10^3 50*10^3]);
hold on
plot([0 0],[0 he_equi],'-r');
plot([d d],[0 hr_equi],'-r');       
plot([0 d],[he_equi hr_equi]);
plot(m_equi,n_equi,'-k');
plot(P_especular,0,'.k','MarkerSize',15);
grid on;    
hold off

rh=(sin(ang_i)-sqrt(no.^2-(cos(ang_i)).^2))/(sin(ang_i)+sqrt(no.^2-(cos(ang_i)).^2))
rv=(no.^2*sin(ang_i)-sqrt(no.^2-(cos(ang_i)).^2))/(no.^2*sin(ang_i)+sqrt(no.^2-(cos(ang_i)).^2))

rte=2*(he_equi*hr_equi)/d
rtp=2*(he*hr)/d

Rh=angle(rh);
Rv=angle(rv);

Coef_Fres_hor=Rh*exp(atan((he+hr)/d))
Coef_Fres_vert=Rv*exp(atan((he+hr)/d))

dif_fase=2*pi/Rh*(d^4/16)*(1/(r*d))*(he/(d^2/4))

Pt_plano=Pd+20*log10(abs(2*sin((2*pi*he*he)/(comp_onda*d))));
Pt_esf=Pd+20*log10(abs(2*sin((2*pi*he_equi*he_equi)/(comp_onda*d))));

Rp=Pre/1013;
Rt=288/(273+Temp);

Ao57=((7.27*Rt)/((57)^2+0.351*Rp^2)+7.5/((57-57)^2+2.44*Rp^2*Rt^5))*57^2*Rp^2*Rt^2*10^-3;
Ao63=((7.27*Rt)/((63)^2+0.351*Rp^2)+7.5/((63-57)^2+2.44*Rp^2*Rt^5))*63^2*Rp^2*Rt^2*10^-3;

if (f/10^9)<=57
    Ao=((7.27*Rt)/((f/10^9)^2+0.351*Rp^2)+7.5/(((f/10^9)-57)^2+2.44*Rp^2*Rt^5))*(f/10^9)^2*Rp^2*Rt^2*10^-3
else
    if 57<(f/10^9) && (f/10^9)<=63
        Ao=(((f/10^9)-60)*((f/10^9)-63)/18)*Ao57-1.66*Rp^2*Rt^8.5*((f/10^9)-57)*((f/10^9)-63)+(((f/10^9)-57)*((f/10^9)-60)/18)*Ao63
    else
       Ao=((2*10^-4*Rt^1.5*(1-1.2*10^-5*(f/10^9)^1.5))+(4/(((f/10^9)-63)^2+1.5*Rp^2*Rt^5))+(0.28*Rt^2)/(((f/10^9)-118.75)^2)+2.84*Rp^2*Rt^2)*(f/10^9)^2*Rp^2*Rt^2*10^-3
    end
end

Aw=((3.24*10^-2*Rt+1.67^-3*Pre*Rt^7/Rp+7.7*10^-4*((f/10^9)/10^9)^0.5)+...
    (3.79/(((f/10^9)-22.235)^2+9.81*Rp^2*Rt))+((11.73*Rt)/((((f/10^9)/10^9)-183.31)^2+11.85*Rp^2*Rt))+...
    ((4.01*Rt)/((((f/10^9)/10^9)-325.153)^2+10.44*Rp^2*Rt)))*((f/10^9)/10^9)^2*Pre*Rt*Rp*10^-4

At=(Ao+Aw)*d

Ri=42;   %Parameters from ITU-R       
p=0.01;  %Parameters from ITU-R
r1=0;       
Kh=0.0101;
Kv=0.00887;
Ah=1.276;
Av=1.264;

k=(Kh+Kv+(Kh-Kv)*cos(ang)^2*cos(2*r1))/2;
a=(Kh*Ah+Kv*Av+(Kh*Ah-Kv*Av)*cos(ang)^2*cos(2*r1))/(2*k);

Yr=k*Ri^a;

d_eficaz=d/(1+(d/(35*exp(-0.015*Ri))));

Ar=Yr*d_eficaz

A_chuva=Ar*0.12*p^-(0.546+0.043*log10(p))

%Modelação dos Efeitos Refrativos na Atmosfera (Ponto 13)

K=1/(1+(r/no)*dndh)

re=(K)*r;
 
y=-((x-(d/2)).^2)/(2*re);

raio_d=sqrt((hr-he)^2+d^2);
raio_menor=sqrt(comp_onda*d/4);

%1Âª Elipse de Fresnel
a=1/2*sqrt((d-0)^2+(hr-he)^2);
b=sqrt(comp_onda*a)/2;

y0=-((0-(d/2)).^2)/(2*re)+he;
y1=-((d-(d/2)).^2)/(2*re)+hr;

t=linspace(0,2*pi,300);
X=a*cos(t);
Y=b*sin(t);

ang=atan((((y1)-(y0)))/d);

m=(0+d)/2+X*cos(ang)-Y*sin(ang);
n=(y0+y1)/2+X*sin(ang)+Y*cos(ang);

syms a;
S=double(real(solve(a.^3-1.5*d*a.^2+(0.5.*d^2-r*(he+hr))*a+he*d*r==0,a,'MaxDegree',3)))

for i=1:length(S)
    if S(i)>0 && S(i)<d
        P_especular=S(i)
    end
end

y1=-((x-(d/2)).^2)/(2*re);

figure()
plot(x,y1);
xlim([-5*10^3 50*10^3]);
hold on
plot([0 0],[-((0-(d/2)).^2)/(2*re) -((0-(d/2)).^2)/(2*re)+he]);
plot([d d],[-((d-(d/2)).^2)/(2*re) -((d-(d/2)).^2)/(2*re)+hr]);       
plot([0 d],[-((0-(d/2)).^2)/(2*re)+he -((d-(d/2)).^2)/(2*re)+hr]);
plot(m,n,'-k');
plot(P_especular,-((P_especular-(d/2)).^2)/(2*re),'.k','MarkerSize',15);
grid on;
hold off


% Assuming Betha =1
X=d*(pi/(comp_onda*r^2))^(1/3);
Y1=2*he*(pi^2/(comp_onda^2*r))^(1/3);
Y2=2*hr*(pi^2/(comp_onda^2*r))^(1/3);

F=11+10*log10(X)-17.6*X;
GY1=17.6*sqrt(Y1-1.1)-5*log10(Y1-1.1)-8
GY2=17.6*sqrt(Y2-1.1)-5*log10(Y2-1.1)-8

Adif=-(F+GY1+GY2)

%Parameters from ITU-R
M=33.20;                    
Y_1=0.27;                   

drhe=sqrt(2*re*he)
drhr=sqrt(2*re*hr)

dm= drhe+drhr

x1=drhe-0;
x2=drhr-d;

y1=(-((drhe-(d/2)).^2)/(2*re))-(re+he);
y2=(-((drhr-(d/2)).^2)/(2*re))-(re+hr);

ang_e=(d*10^3/(K*re));
ang_t=acos(((y1*y2)+(x1*x2))/(x1*x2));
ang_r=ang_t;

ang_retas=real(ang_t+ang_r+ang_e);

H=(ang_retas*d)/(4*10^3)
h=(ang_retas^2*re)/(8*10^6)

NHh=20*log10(5+Y_1*H)+4.34*Y_1*h;

Lc=0.07*exp(0.055*(Ge+Gr));
Gp=Ge+Gr-Lc

A=30*log10(fmhz)+10*log10(dkm)+30*log10(ang_retas)+M+NHh+Gp
