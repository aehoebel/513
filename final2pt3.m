clear
n=4000;
w=[0.005,0.0028,0.0005];
L=1;
C=1;
%Z_L=300;
V_0=1;
%Z_L=300;
Z_L=sqrt(L/C);

Z1=zeros(1,n);
Z1(1)=Z_L;
Z2=zeros(1,n);
Z2(1)=Z_L;
Z3=zeros(1,n);
Z3(1)=Z_L;

for i=1:n
Z1(i+1)=(j*w(1)*L)+(Z1(i))/(((j*w(1)*C)*Z1(i))+1);
Z2(i+1)=(j*w(2)*L)+(Z2(i))/(((j*w(2)*C)*Z2(i))+1);
Z3(i+1)=(j*w(3)*L)+(Z3(i))/(((j*w(3)*C)*Z3(i))+1);
end

Z_4 = sqrt((L/C)-((w(1).^2)*(L.^2)/4));
Z4 = zeros(1,n+1) + Z_4;

x = linspace(0,n,length(Z1)); % x coordinate of spatial samples
figure(1);
plot(x,Z1,'g');
hold on;
plot(x,Z2,'b');
hold on;
plot(x,Z3,'r');
hold on;
plot(x,Z4,'k');
title('Discrete Approximation vs Exact Solution of Continuous Transmission Line');
legend('\omega = 0.0050','\omega = 0.0028','\omega = 0.0005','Exact solution');
xlabel('N');
ylabel('Impedance');
ylim([0.99 1.01]);
print('Final_Part_II_Figure_5','-dpdf','-fillpage')
drawnow;
