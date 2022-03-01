clc
clear

Lsum=1.5;

rpm=3600;

beta2=70*pi()/180;

L1=0.5;
    
L0=Lsum-L1;

Wa=2*pi()*rpm*L0/(60*tan(beta2));

V1=Wa;

segments=3;

unit=L1/(segments-1);

w=503707.267019818;

mdot=5000000/w;

alpha1=0;

for Lkar=1:segments

y(Lkar)=unit*(Lkar-1);

U(Lkar)=rpm*2*pi()*(y(Lkar)+L0)/60;

beta3(Lkar)=-(-1*(20*pi()/(L1*180))*y(Lkar)+beta2);

Wu2(Lkar)=tan(beta2)*Wa;
W2(Lkar)=sqrt(Wu2(Lkar)^2+Wa^2);

alpha2(Lkar)=atan((Wu2(Lkar)+rpm*2*pi()*(y(Lkar)+L0)/60)/Wa);
V2(Lkar)=Wa/cos(alpha2(Lkar));

W3(Lkar)=Wa/(cos(beta3(Lkar)));
Wu3(Lkar)=Wa*tan(beta3(Lkar));

alpha3(Lkar)=-atan(tan(beta3(Lkar))+(2*pi()*3600/(Wa*60))*(y(Lkar)+L0));

V3(Lkar)=Wa/cos(alpha3(Lkar));

reactions (Lkar)=-(W3(Lkar)^2-W2(Lkar)^2)/(W3(Lkar)^2-W2(Lkar)^2+V2(Lkar)^2-V1^2);

Vu3=(Wu3+U);

Vu2=Wu2+U;

end

beta2=zeros(1, segments);

beta2(1,:)=70*pi()/180;

%now for finding the losses

T01=zeros(1, segments);
T03=zeros(1, segments);

T01(1,:)=1600;
T02=T01;
T03(1,:)=1400;

h01=zeros(1, segments);
h03=zeros(1, segments);

h01(1,:)=1757.57*1000;
h03(1,:)=1400*1000;

%the cp that will be used is the mean of both temperatures it must be in
%J/(kg.K)
cp=(h01./T01+h03./T03)./2;

%usually R stays constant, therefore gamma changes
R=287;
gamma=cp/(cp-R);

%now to find T1, and T3
T1=T01-(V1.^2)./(2*cp);
T2=T02-(V2.^2)./(2*cp);
T3=T03-(V3.^2)./(2*cp);

%finsing ksi and psi
psi=Wa./U;
ksi=2*(1-reactions-psi.*tan(alpha3));

for Lkar=1:1:segments-2
    perlocation=Lkar*1/(segments-1);
location(Lkar)=round(perlocation*100,1)+"% of the blade";
end

%defining the average temperature and velocity

beta2prime=70*pi()/180;
segmentsprime=10000;
unitprime=L1/(segmentsprime-1);
for Lkar=1:segmentsprime
yprime(Lkar)=unitprime*(Lkar-1);
Uprime(Lkar)=rpm*2*pi()*(yprime(Lkar)+L0)/60;
beta3prime(Lkar)=-(-1*(20*pi()/(L1*180))*yprime(Lkar)+beta2prime);
Wu2prime(Lkar)=tan(beta2prime)*Wa;
alpha2prime(Lkar)=atan((Wu2prime(Lkar)+rpm*2*pi()*(yprime(Lkar)+L0)/60)/Wa);
V2prime(Lkar)=Wa/cos(alpha2prime(Lkar));
end

T01prime=zeros(1, segmentsprime);
T03prime=zeros(1, segmentsprime);

T01prime(1,:)=1600;
T02prime=T01prime;
T03prime(1,:)=541.6;

h01prime=zeros(1, segmentsprime);
h03prime=zeros(1, segmentsprime);

h01prime(1,:)=1757.57*1000;
h03prime(1,:)=1400*1000;

%the cp that will be used is the mean of both temperatures it must be in
%J/(kg.K)
cpprime=(h01prime./T01prime+h03prime./T03prime)./2;
V2avg=mean(V2prime);
T2avg=mean(T02prime-(V2prime.^2)./(2*cpprime));

V1avg=mean(V1);
V3avg=mean(V3);

%defining the design and shape of the blade
L1array=zeros(1,segments);
L1array(1,:)=L1;
cx=L1array/3;

%for c from an angle of 58 for xi
c=cx/cos(58*pi()/180);

%calculating the density from the thermodynamic table
%rho=0.282628;
rho1=mdot/(V1avg*(Lsum^2-L0^2)*pi());
rho2=mdot/(V2avg*(Lsum^2-L0^2)*pi());
rho3=mdot/(V3avg*(Lsum^2-L0^2)*pi());

%now to calculate p and subsequently mu with P=rho*R*T
P1=rho1*R*T1/100;
P2=rho2*R*T2/100;
P3=rho3*R*T3/100;

%at a temperature of around 1250 K or 977 cx, the viscosity by interpolation in the thermodynamic table, the viscosity is 4.77402*10^(-5) kg/(m*s)
mu=4.77402*10^(-5);

%"s" will be defined from the solidity using zweifel's corolation
alpha3avg=mean(alpha3);
alpha2avg=mean(alpha2);
beta3avg=mean(beta3);
V3avg=mean(V3);
ksiavg=mean(ksi);
cxavg=mean(cx);

%to find the stagnation pressure, we already have P2, but for P02 we have
%to add to it 0.5rhoV2square
P01=P1+(V1.^2)*0.5*rho1;
P02=P2+(V2.^2)*0.5*rho2;
P03=P3+(V3.^2)*0.5*rho3;

%finding "s"
s=mean(ksi.*cx.*(V3.^2).*(0.5)./((Wa.^2).*(tan(alpha2)-tan(alpha3))));
s2=mean(ksi.*cx.*(0.5)./((cos(alpha3).^2).*(tan(alpha2)-tan(alpha3))));
pitch=zeros(1, segments);

pitch(1,:)=s;

%now to calculate the reynold's number
Dhstator=2*L1*s*cos(alpha2avg)/(s*cos(alpha2avg)+L1);
Dhrotor=2*L1*s*cos(beta3avg)/(s*cos(beta3avg)+L1);

Restator=rho2*V2avg*Dhstator/mu;
Rerotor=rho2*V2avg*Dhrotor/mu;

%now calculating zeta
zetarotor=(0.04+0.06*((beta2-beta3)*180/(2*pi())/100).^2)*(((10^5)/Rerotor)^(1/4));
zetastator=(0.04+0.06*((alpha2-alpha1)*180/(2*pi())/100).^2)*(((10^5)/Restator)^(1/4));

M2=V2./(sqrt(R*gamma*T2));
M3=V3./(sqrt(R*gamma*T3));

%pressure losses
PLOS=(gamma./2).*P02.*zetastator.*(M2.^2);
PLOR=(gamma./2).*P03.*zetarotor.*(M3.^2);

Ypa=(-0.627*(alpha2*180/(pi())/100).^2+0.821*(alpha2*180/(pi())/100)-0.129).*(s./c).^2+(1.489*(alpha2*180/(pi())/100).^2-1.676*(alpha2*180/(pi())/100)+0.242).*(s./c)+(-0.356*(alpha2*180/(pi())/100).^2+0.399*(alpha2*180/(pi())/100)+0.0077);
Ype=(-1.56*(alpha2*180/(pi())/100).^2+1.55*(alpha2*180/(pi())/100)-0.064).*(s./c).^2+(3.73*(alpha2*180/(pi())/100).^2-3.43*(alpha2*180/(pi())/100)+0.290).*(s./c)+(-0.83*(alpha2*180/(pi())/100).^2+0.78*(alpha2*180/(pi())/100)+0.078);

alpham=atan(0.5*(tan(alpha2)+tan(alpha1)));
betam=atan(0.5*(tan(beta2)+tan(beta3)));

C_Lstator=2*(s./c).*(tan(alpha2)-tan(alpha1)).*tan(alpham);
C_Lrotor=2*(s./c).*(tan(beta2)-tan(beta3)).*tan(betam);

%shrouded blades
B=0.37;

Y_s_Y_k_rotor=(c./L1).*(0.0334*cos(beta3)./cos(beta2) +B.*(0.02).^0.78 ).*(C_Lrotor./(s./c)).*(cos(beta3).^2)./(cos(betam).^3); 
Y_stator=(c./L1).*(0.0334*cos(alpha2)./cos(alpha1)).*(C_Lstator./(s./c)).*(cos(alpha2).^2)./(cos(betam).^3); 

%finding Fu or driving force
Fu=mdot*ksi.*U;


mdotprime=zeros(1, segments);
mdotprime(1,:)=mdot;

P01prime=P01;
P01=zeros(1, segments);
P01(1,:)=P01prime;

P03prime=P03;
P03=zeros(1, segments);
P03(1,:)=P03prime;

pressureratio=P01./P03;
polytropiceffeciency=(log((0.9*(1-(1/mean(pressureratio))^((gamma-1)/gamma))-1)*-1)/log(1/mean(pressureratio)))*gamma/(gamma-1);

polytropiceffeciencyprime=polytropiceffeciency;
polytropiceffeciency=zeros(1, segments);
polytropiceffeciency(1,:)=polytropiceffeciencyprime;

all_components=["properties", "hub", location, "casing";
   "y (m)", y;
   "mass flow rate (kg/s)", mdotprime;
   "beta3 (rad)", beta3; 
   "beta2 (rad)", beta2;
   "U (m/s)", U; 
   "Wu2 (m/s)", Wu2; 
   "W2 (m/s)", W2; 
   "alpha2 (degree)", alpha2; 
   "V2 (m/s)", V2; 
   "W3 (m/s)", W3;
   "Wu3 (m/s)", Wu3;
   "alpha3 (degree)", alpha3;
   "V3 (m/s)", V3;
   "Vu3 (m/s)", Vu3;
   "Vu2 (m/s)", Vu2;
   "R", reactions;
   "T1 (K)", T1;
   "T2 (K)", T2;
   "T3 (K)", T3;
   "P1 (kPa)", P1;
   "P2 (kPa)", P2;
   "P3 (kPa)", P03;
   "T01 (K)", T01;
   "T02 (K)", T02;
   "T03 (K)", T03;
   "P01 (kPa)", P01;
   "P02 (kPa)", P02;
   "P03 (kPa)", P03;
   "ksi", ksi;
   "psi", psi;
   "Fu (N)", Fu;
   "blade length (m)", L1array;
   "blade chord (m)", cx;
   "pitch", pitch;
   "zeta stator", zetastator;
   "zeta rotor", zetarotor;
   "mach number at 2", M2;
   "mach number at 3", M3;
   "pressure loss at the stator (kPa)", PLOS;
   "pressure loss at the rotor (kPa)", PLOR;
   "profile losses stator", Ypa;
   "profile losses rotor", Ype;
   "secondary flow and clearance tip losses rotor",Y_s_Y_k_rotor;
   "secondary flow losses stator",Y_stator;
   "total loss coefficient stator",Ypa+Y_stator ;
   "total loss coefficient stator",Ype+Y_s_Y_k_rotor;
   "pressure ratio", pressureratio;
   "polytropic effeciency", polytropiceffeciency;];