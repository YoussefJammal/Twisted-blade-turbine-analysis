clc
clear

Lsum=1.5;

rpm=3600;

beta2=70*pi()/180;

Wa=2*pi()*rpm*Lsum/(60*tan(pi()*70/180));


for Lkar=1:1:1
clc
clear

Lsum=1.5;

rpm=3600;

beta2=70*pi()/180;

Wa=2*pi()*rpm*Lsum/(60*tan(pi()*70/180));


for Lkar=1:1:1

h0=Lkar;    
    
L1=0.5;
    
L0=Lsum-L1;

f1=@(x)(x+L0).*(tan(beta2)+tan(-(20*pi()/(L1*180))*x+beta2));

q1=integral(f1,0,L1);

fun1=q1*rpm*Wa*2*pi()/(60);

f2=@(x)cos(-(20*pi()/(L1*180))*x+beta2);

q2=integral(f2,0,L1);

fun2=1757.57*1000-1515.42*1000;

test1(Lkar)=fun1;
test2(Lkar)=fun2;
test3(Lkar)=q2/L1;

eta(Lkar)=fun1/fun2;

end

disp("done")
h0=Lkar;    
    
L1=0.5;
    
L0=Lsum-L1;

f1=@(x)(x+L0).*(tan(beta2)+tan(-(20*pi()/(L1*180))*x+beta2));

q1=integral(f1,0,L1);

fun1=q1*rpm*Wa*2*pi()/(60);

f2=@(x)cos(-(20*pi()/(L1*180))*x+beta2);

q2=integral(f2,0,L1);

fun2=1757.57*1000-1515.42*1000;

test1(Lkar)=fun1;
test2(Lkar)=fun2;
test3(Lkar)=q2/L1;

eta(Lkar)=fun1/fun2;

end

disp("done")