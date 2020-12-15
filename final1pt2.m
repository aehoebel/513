clear
n=1000;
w=0.005;
L=1;
C=1;
V_0=1;
Z_L=300;
t=(2*pi/w)/4;

Z=zeros(1,n);
Z(1)=Z_L;

for i=1:n
Z(i+1)=(j*w*L)+(Z(i))/(((j*w*C)*Z(i))+1);
end

V=zeros(1,n);
I=zeros(1,n);
V(1)=V_0;
I(1)=V_0/Z(n);
for i=1:n-1
%I(i)=V(i)/Z(i);
V(i+1)=(V(i)*(1-w^2*C*L)-j*w*C*I(i))*exp(j*w*t);
%V(i+1)=V(i)-j*w*L*I(i);
I(i+1)=(I(i)-j*w*C*V(i+1))/(1-w^2*C*L);
end
%V
%I
x = linspace(0,n,length(V)); % x coordinate of spatial samples
figure(1);
plot(x,V,'k');
title('V_k(t) of Ladder Network with 2\Pi/\omega/4');
xlabel('N');
ylabel('Impedance');
ylim([-2 2]);
print('Final_Part_II_Figure_2','-dpdf','-fillpage')
drawnow;
