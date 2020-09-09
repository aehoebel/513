
V = [80,60,100,20];

x = [1/3,2/3,1/3,2/3];

y = [2/3,2/3,1/3,1/3];

v = [0,0,0,0];

v(1)=(V(3)+V(1)+((V(1)+V(2)+V(3)+V(4))/4)+((V(1)+V(2)+V(3)+V(4))/4))/4;
v(2)=(V(3)+V(4)+v(1)+((V(1)+V(2)+V(3)+V(4))/4))/4;
v(3)=(v(1)+V(1)+((V(1)+V(2)+V(3)+V(4))/4)+V(2))/4;
v(4)=(V(4)+V(2)+v(2)+v(3))/4;

fprintf("\nThe potentials from columb 1 of the table:\n");
v

v = [0,0,0,0];
for i=1:1:4
vl_s = 0;
vb_s = 0;
vt_s = 0;
vr_s = 0;
    for n=1:2:11
    vl=(4*V(1)/(pi)).*sinh(n*pi*(1-x(i))).*sin(n*pi*y(i))./(n*sinh(n*pi));
    vb=(4*V(2)/(pi)).*sinh(n*pi*(1-y(i))).*sin(n*pi*x(i))./(n*sinh(n*pi));
    vt=(4*V(3)/(pi)).*sinh(n*pi*y(i)).*sin(n*pi*x(i))./(n*sinh(n*pi));
    vr=(4*V(4)/(pi)).*sinh(n*pi*x(i)).*sin(n*pi*y(i))./(n*sinh(n*pi));
    vl_s=vl_s+vl; %summation of each value at n
    vb_s=vb_s+vb;
    vt_s=vt_s+vt;
    vr_s=vr_s+vr;
    end
v(i) = vl_s + vb_s + vt_s + vr_s;
end

fprintf("My attempt at the correct potentials using the solutions from 2.1.1:\n");
v