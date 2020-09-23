clc;
clear all;
close all;

p = [0, 0, 0];
I = 1;

XYZ = [1000, 1, 0; -1000,1,0];
fprintf('\nCalculated B value of 1 line segment of very long length:\n');
B = HW4_BiotSavart(p, XYZ, I)

fprintf('\nActual B value of 1 line segment of very long length at r=1 and I=1:\n');
u0=4*pi*10.^-7; %u0 value
r=1;
B = u0*I/(2*pi*r)

n=4;
r=1;
XYZ=zeros(n,3);
for i=1:n
  XYZ(i,1)=cos(2 * pi * i / n);
  XYZ(i,2)=r * sin(2 * pi * i / n);
  XYZ(i,3)=0;
end
B1 = HW4_BiotSavart(p, XYZ, I);

n=6;
r=1;
XYZ=zeros(n,3);
for i=1:n
  XYZ(i,1)=cos(2 * pi * i / n);
  XYZ(i,2)=r * sin(2 * pi * i / n);
  XYZ(i,3)=0;
end
B2 = HW4_BiotSavart(p, XYZ, I);

n=8;
r=1;
XYZ=zeros(n,3);
for i=1:n
  XYZ(i,1)=cos(2 * pi * i / n);
  XYZ(i,2)=r * sin(2 * pi * i / n);
  XYZ(i,3)=0;
end
B3 = HW4_BiotSavart(p, XYZ, I);

n=10;
r=1;
XYZ=zeros(n,3);
for i=1:n
  XYZ(i,1)=cos(2 * pi * i / n);
  XYZ(i,2)=r * sin(2 * pi * i / n);
  XYZ(i,3)=0;
end
B4 = HW4_BiotSavart(p, XYZ, I);

n=100;
r=1;
XYZ=zeros(n,3);
for i=1:n
  XYZ(i,1)=cos(2 * pi * i / n);
  XYZ(i,2)=r * sin(2 * pi * i / n);
  XYZ(i,3)=0;
end
B5 = HW4_BiotSavart(p, XYZ, I);

n=1000;
r=1;
XYZ=zeros(n,3);
for i=1:n
  XYZ(i,1)=cos(2 * pi * i / n);
  XYZ(i,2)=r * sin(2 * pi * i / n);
  XYZ(i,3)=0;
end
B6 = HW4_BiotSavart(p, XYZ, I);

Sides=[4;6;8;10;100;1000];
B=[B1;B2;B3;B4;B5;B6];

fprintf('\nTable of B values for polygons of n number of sides\n');
T = table(Sides,B)

fprintf('Higher n sides produce a B value close to that of a circle.\n\nActual B value of a circle of r=1, I=1, z=0:\n');
B = u0*I/(2*r)


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
            r0=[XYZ(n,1), XYZ(n,2), p(3)];
        else
            L=[abs(XYZ(N,1)-XYZ(1,1)),abs(XYZ(N,2)-XYZ(1,2)),abs(XYZ(N,3)-XYZ(1,3))]; %to calculate length of line segment at end of array
            r0=[XYZ(1,1), XYZ(1,2), p(3)];
        end
        r=sqrt(sum(r0.^2));
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

