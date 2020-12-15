clear
n=3;
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

y=200;
Zb = [y, y/2, y/3, y/4];
Z2(length(Zb)) = Zb(end);
d = [1/4,1/8,1/16,1/20];

for i = [length(Zb):-1:2]
    r = Z2(i)/Zb(i-1);
    R = exp(-4*pi*1i*d(i-1))*(r - 1)/(r + 1);
    Z2(i-1) = Zb(i-1)*(1 + R)/(1 - R);
end

x = linspace(0,n,length(Z)); % x coordinate of spatial samples
figure(1);
plot(x,Z,'k');
hold on;
plot(x,Z2,'b');
title('Discrete Approximation vs Exact Solution of Continuous Transmission Line');
legend('Discrete approximation with impedance Z_L \neq \surd(L/C)','Exact solution with impedances of [200, 200/2, 200/3, 200/4]');
xlabel('N');
ylabel('Impedance');
print('Final_Part_II_Figure_4','-dpdf','-fillpage')
drawnow;

