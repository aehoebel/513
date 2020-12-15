clear
n=4000;
w=0.005;
L=1;
C=1;
V_0=1;
Z_L=sqrt(L/C);
t=(2*pi/w)/4;

Z=zeros(1,n);
Z(1)=Z_L;

for i=1:n
Z(i+1)=(j*w*L)+(Z(i))/(((j*w*C)*Z(i))+1);
end

Zb = ones(1,n+1);
Z2(length(Zb)) = Zb(end);
for i = [1:i]
d = ones(1,n+1)+i*(1/4);
end

for i = [length(Zb):-1:2]
    r = Z2(i)/Zb(i-1);
    R = exp(-4*pi*1i*d(i-1))*(r - 1)/(r + 1);
    Z2(i-1) = Zb(i-1)*(1 + R)/(1 - R);
end

x = linspace(0,n,length(Z)); % x coordinate of spatial samples
figure(1);
plot(x,Z,'b');
hold on;
plot(x,Z2,'k');
title('Discrete Approximation vs Exact Solution of Continuous Transmission Line');
legend('Discrete approximation with impedance Z_L = \surd(L/C)','Exact solution');
xlabel('N');
ylabel('Impedance');
ylim([0.95 1.05]);
print('Final_Part_II_Figure_3','-dpdf','-fillpage')
drawnow;

