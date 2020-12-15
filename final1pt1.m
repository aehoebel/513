clear
n=1000;
w=0.005;
L=1;
C=1;
c=1/(j*w*C);
l=j*w*L;
V_0=1;
Z_L=300;
t=0;
%t=(2*pi/w)/4;

Z=zeros(1,n);
Z(1)=Z_L;

for i=1:n
%Z(i+1)=(j*w*L)+(Z(i))/(((j*w*C)*Z(i))+1);
Z(i+1)=l+(Z(i)*c)/(c+Z(i));
end

V=zeros(1,n);
I=zeros(1,n);
V(1)=V_0;
I(1)=V_0/Z(n);

for i=1:n-1
%V(i+1)=(V(i)*(1-w^2*C*L)-j*w*C*I(i))*exp(j*w*t);
V(i+1)=(V(i)-I(i)*l)*exp(j*w*t);
%V(i)=(V_0*(c+l+Z(i))/(l*(2*c+l+Z(i))+c*Z(i)))*exp(j*w*t);
%V(i+1)=V(i)-j*w*L*I(i);
I(i+1)=((I(i)-j*w*C*V(i+1))/(1-w^2*C*L));
end
%V
%I
x = linspace(0,n,length(V)); % x coordinate of spatial samples
figure(1);
plot(x,V,'k');
title('V_k(t) of Ladder Network with t=0');
xlabel('N');
ylabel('Impedance');
print('Final_Part_II_Figure_1','-dpdf','-fillpage')
%ylim([0.95 1.05]);
drawnow;
