function deri = fourtank1(t,h,v,w)
h1=h(1);
h2=h(2);
h3=h(3);
h4=h(4);
A1=28;
A3=A1;
A2=32;
A4=A2;
a1=0.071;
a3=a1;
a2=0.057;
a4=a2;
kc=0.5; g=981;
x1=0.7; x2= 0.6;
v1=v(1); v2=v(2);
k1=3.33; k2=3.35;

deri(1,1)= (-a1/A1)*(sqrt(2*g*h1))+ (a3/A1)*(sqrt(2*g*h3)) + (x1*k1*v1/A1)+w(1);
deri(2,1)= (-a2/A2)*(sqrt(2*g*h2))+ (a4/A2)*(sqrt(2*g*h4)) + (x2*k2*v2/A2)+w(2);
deri(3,1)= (-a3/A3)*(sqrt(2*g*h3))+ ((1-x2)*k2*v2/A3)+w(3);
deri(4,1)= (-a4/A4)*(sqrt(2*g*h4))+ ((1-x1)*k1*v1/A4)+w(4);