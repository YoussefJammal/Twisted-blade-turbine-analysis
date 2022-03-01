clc
clear

Lsum=1.5;

rpm=3600;

beta2=70*pi()/180;

L1=0.5;
    
L0=Lsum-L1;

Wa=2*pi()*rpm*L0/(60*tan(beta2));

V1=Wa;

segments=10000;

unit=L1/(segments-1);

w=503707.267019818;

for Lkar=1:segments

y=unit*(Lkar-1);

U(Lkar)=rpm*2*pi()*(y+L0)/60;

beta3(Lkar)=-(-1*(20*pi()/(L1*180))*y+beta2);

Wu2(Lkar)=tan(beta2)*Wa;
W2(Lkar)=sqrt(Wu2(Lkar)^2+Wa^2);

alpha2(Lkar)=atan((Wu2(Lkar)+rpm*2*pi()*(y+L0)/60)/Wa);
V2(Lkar)=Wa/cos(alpha2(Lkar));

W3(Lkar)=Wa/(cos(beta3(Lkar)));
Wu3(Lkar)=Wa*tan(beta3(Lkar));

alpha3(Lkar)=atan(tan(beta3(Lkar))+(2*pi()*3600/(Wa*60))*(y+L0));

V3(Lkar)=Wa/cos(alpha3(Lkar));

reactions (Lkar)=(W3(Lkar)^2-W2(Lkar)^2)/(W3(Lkar)^2-W2(Lkar)^2+V2(Lkar)^2-V1^2);

Vu3=(Wu3+U);

Vu2=Wu2+U;

end

disp("done")

k=linspace(0,L1,segments);

figure (1)
plot(k,alpha2*180/pi())
title('alpha2')
xlabel('position (m)')
ylabel('angle (degree)')

figure (2)
plot(k,alpha3*180/pi())
title('alpha3')
xlabel('position (m)')
ylabel('angle (degree)')

figure (3)
plot(k,W2)
title('W2')
xlabel('position (m)')
ylabel('velocity (m/s)')

figure (4)
plot(k,V2)
title('V2')
xlabel('position (m)')
ylabel('velocity (m/s)')

figure (5)
plot(k,W3)
title('W3')
xlabel('position (m)')
ylabel('velocity (m/s)')

figure (6)
plot(k,V3)
title('V3')
xlabel('position (m)')
ylabel('velocity (m/s)')

figure (7)
plot(k,U)
title('U')
xlabel('position (m)')
ylabel('velocity (m/s)')

figure (8)
plot(k,abs(reactions))
title('reactions')
xlabel('position (m)')
ylabel('reactin (%)')

figure (9)
plot(k,Wu2)
title('Wu2')
xlabel('position (m)')
ylabel('velocity (m/s)')

figure (10)
plot(k,Wu3)
title('Wu3')
xlabel('position (m)')
ylabel('velocity (m/s)')

figure (11)
plot(k,beta3*180/pi())
title('beta3')
xlabel('position (m)')
ylabel('angle (degree)')

figure (12)
plot(k,Vu2)
title('Vu2')
xlabel('position (m)')
ylabel('velocity (m/s)')

figure (13)
plot(k,Vu3)
title('Vu3')
xlabel('position (m)')
ylabel('velocity (m/s)')