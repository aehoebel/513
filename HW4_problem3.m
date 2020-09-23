clc;
clear all;
close all;

n=1000;
r=1;
XYZ=zeros(n,3);
for i=1:n
  XYZ(i,1)=cos(2 * pi * i / n);
  XYZ(i,2)=r * sin(2 * pi * i / n);
  XYZ(i,3)=0;
end

p = [0, 0, -2];
I = 1;
B = HW4_BiotSavart(p, XYZ, I);
f1 = B*pi*(r^2);

p = [0, 0, 0];
I = 1;
B = HW4_BiotSavart(p, XYZ, I);
f2 = B*pi*(r^2);

fprintf('\nCalculated flux of 2 rings at z=0 and z=2:\n');
F = f1+f2

r=1;
z=2;
u0=4*pi*10.^-7; %u0 value
B1 = u0*I*(r.^2)./(2*((r.^2+z.^2).^(3/2)));
f1 = B1*pi*(r^2);
z=0;
u0=4*pi*10.^-7; %u0 value
B2 = u0*I*(r.^2)./(2*((r.^2+z.^2).^(3/2)));
f2 = B2*pi*(r^2);

fprintf('\nFlux of 2 rings at z=0 and z=2 using ring equation:\n');
F = f1+f2

function[BT] = HW4_BiotSavart(p, XYZ, I)
N=size(XYZ,1); %how many line segments are there?
u0=4*pi*10.^-7; %u0 value
B0=[0,0,0]; %to calculate magnetic field vector
B=0;
r=1;

if N==2 %Calculates magnetic value if there is only 1 line segment
    L=[abs(XYZ(2,1)-XYZ(1,1)),abs(XYZ(2,2)-XYZ(1,2)),abs(XYZ(2,3)-XYZ(1,3))]; %to calculate length of line segments
    r0=[((XYZ(2,1)+XYZ(1,1))-p(1))./2,((XYZ(2,2)+XYZ(1,2))-p(2))./2,((XYZ(2,3)+XYZ(1,3))-p(3))./2];
    r=sqrt(sum(r0.^2));
    for i=1:3
        B0(i)=(u0/(4*pi*r))*I.*(L(i)./((r.^2 + (L(i)./2).^2).^(1/2)));
    end
    B=B+B0;
else    
    for n=1:N
        if n~=N
            L=[abs(XYZ(n+1,1)-XYZ(n,1)),abs(XYZ(n+1,2)-XYZ(n,2)),abs(XYZ(n+1,3)-XYZ(n,3))]; %to calculate length of line segments
            r0=[abs(XYZ(N,1)+XYZ(1,1))./2,abs(XYZ(N,2)+XYZ(1,2))./2, p(3)];
        else
            L=[abs(XYZ(N,1)-XYZ(1,1)),abs(XYZ(N,2)-XYZ(1,2)),abs(XYZ(N,3)-XYZ(1,3))]; %to calculate length of line segment at end of array
            r0=[abs(XYZ(N,1)+XYZ(1,1))./2,abs(XYZ(N,2)+XYZ(1,2))./2, p(3)];
        end
        r=sqrt(sum((r0.^2)));
        for i=1:3
            B0(i)=(u0/(2*pi*r))*I./(((r./(L(i)/2)).^2 + 1).^(1/2));
        end
    B=B+B0;
    end
end
if N==2
    BT=(N-1)*sqrt(sum(B0.^2));
else
    BT=N*sqrt(sum(B0.^2));
end
end