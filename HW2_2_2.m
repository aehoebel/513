
V = [80,60,100,20];

x = [1/3,2/3,1/3,2/3];

y = [2/3,2/3,1/3,1/3];

v = [0,0,0,0];

for i=1:1:4
vl=(4*V(1)/(pi))*sinh(pi*(1-x(i)))*sin(pi*y(i))./sinh(pi);
vb=(4*V(2)/(pi))*sinh(pi*(1-y(i)))*sin(pi*x(i))./sinh(pi);
vt=(4*V(3)/(pi))*sinh(pi*y(i))*sin(pi*x(i))./sinh(pi);
vr=(4*V(4)/(pi))*sinh(pi*x(i))*sin(pi*y(i))./sinh(pi);
v(i) = vl + vb + vt + vr;
end

v