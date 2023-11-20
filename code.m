clear al; close all; clc;


Pe=10;                                  % Emisson Power
he=40;                                  % Height of emitter
hr=100;                                 % Height of receptor
f=10*10^9;                              % Frequence in Hertz
d=45000;                                % Distance in meters
comp_onda=3*10^8/f;                     % Wavelength
l=3;                                    % Diameter of the antenna's dish
rend=0.5;                               % Efficiency
Ag=44.3;                                % Atenuation dos guias por km
Dg=60;                                  % Length dos guias em m
r = 6370*10^3;                          % Earth's ray
%AtenuaÃ§Ã£o dos guias
ar=10^-((Ag*Dg/1000)/10)               
ae=ar

Pre=1013;                               %Atmospheric pressure in hPa
Temp=25;                                %Temperature in Celsius Degree

%AtenuaÃ§Ã£o do espaÃ§o livre
A0=(comp_onda/(4*pi*d))^2   

%Cálculo do ângulo de fogo
ang = atan((hr-he)/d)

%Ganhos das antenas
Gr=20*log10(pi*l*f*1000/300)+10*log10(rend) 
Ge=Gr

x=0:45000;
y=-((x-(d/2)).^2)/(2*r);

raio_d=sqrt((hr-he)^2+d^2);
raio_menor=sqrt(comp_onda*d/4);

%1Âª Elipse de Fresnel
a=1/2*sqrt((d-0)^2+(hr-he)^2);
b=sqrt(comp_onda*a)/2;

t=linspace(0,2*pi,300);
X=a*cos(t);
Y=b*sin(t);
ang=atan((he+hr-80)/d);

m=(0+d)/2+X*cos(ang)-Y*sin(ang);
n=(he+hr-80)/2+X*sin(ang)+Y*cos(ang);

%tangente do emissor
tanE = d*x/(2*r) + he - d^2/(8*r);

%tangente do Recetor
tanR = -d*x/(2*r) + (3/8)* d^2/r + hr;

%ponto especular terra esferica
syms a;
S=double(real(solve(a.^3-1.5*d*a.^2+(0.5.*d^2-r*(he+hr))*a+he*d*r==0,a,'MaxDegree',3)))

for i=1:length(S)
    if S(i)>0 && S(i)<d
        P_especular=S(i)
    end
end

%Distâncias d1 e d2
d1=P_especular;
d2=d-d1;

%Alturas equivalentes
he_equi=he-d1^2/(2*r)
hr_equi=hr-d2^2/(2*r)

%ângulo de incidência
ang_incidencia=atan(he_equi/d1);

%fator de divergência
D=1/(sqrt(1+(2*d1*d2)/(r*d*sin(ang_incidencia))))

%figuras
figure()
plot(x,y);
xlim([-5*10^3 50*10^3]);
hold on
plot(tanE);
plot(tanR);
plot([0 0],[-40 -40+he]);
plot([d d],[-40 -40+hr]);       
plot([0 d],[-40+he -40+hr]);
plot(m,n,'-k');
plot(P_especular,-((P_especular-(d/2)).^2)/(2*r),'.k','MarkerSize',15);
grid on;
hold off

%Fresnel para altura equivalente

a2=1/2*sqrt((d-0)^2+(hr_equi-he_equi)^2);
b2=sqrt(comp_onda*a2)/2;

t1=linspace(0,2*pi,300);
X1=a2*cos(t1);
Y1=b2*sin(t1);
ang_i=atan((he_equi+hr_equi-40)/d);

m_equi=(d)/2+X*cos(ang_i)-Y1*sin(ang_i);
n_equi=(he_equi+hr_equi)/2+X1*sin(ang_i)+Y1*cos(ang_i);

%figura do equivalente em terra plana
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

%Coeficiente de Fresnel
Rr=4; %?? página 12 reflexão no solo

Coef_Fres=Rr*exp(atan((he+hr)/d))

%Comparação de  fase da Terra Esférica - Terra Plana
dif_fase=2*pi/comp_onda*(d^4/16)*(1/(r*d))*(he/(d^2/4))

%Atenuação vapor de água e do oxigénio
Rp=Pre/1013;
Rt=288/(273+Temp);

Ao57=((7.27*Rt)/((57)^2+0.351*Rp^2)+7.5/((57-57)^2+2.44*Rp^2*Rt^5))*57^2*Rp^2*Rt^2*10^-3;
Ao63=((7.27*Rt)/((63)^2+0.351*Rp^2)+7.5/((63-57)^2+2.44*Rp^2*Rt^5))*63^2*Rp^2*Rt^2*10^-3;

if f<=57
    Ao=((7.27*Rt)/((f)^2+0.351*Rp^2)+7.5/((f-57)^2+2.44*Rp^2*Rt^5))*f^2*Rp^2*Rt^2*10^-3;
else
    if 57<f && f<=63
        Ao=((f-60)*(f-63)/18)*Ao57-1.66*Rp^2*Rt^8.5*(f-57)*(f-63)+((f-57)*(f-60)/18)*Ao63;
    else
       Ao=((2*10^-4*Rt^1.5*(1-1.2*10^-5*f^1.5))+(4/((f-63)^2+1.5*Rp^2*Rt^5))+(0.28*Rt^2)/((f-118.75)^2)+2.84*Rp^2*Rt^2)*f^2*Rp^2*Rt^2*10^-3;
    end
end

Aw=((3.24*10^-2*Rt+1.67^-3*Pre*Rt^7/Rp+7.7*10^-4*f^0.5)+(3.79/((f-22.235)^2+9.81*Rp^2*Rt))+((11.73*Rt)/((f-183.31)^2+11.85*Rp^2*Rt))+...
    ((4.01*Rt)/((f-325.153)^2+10.44*Rp^2*Rt)))*f^2*Pre*Rt*Rp*10^-4;

At=(Ao+Aw)*d;
