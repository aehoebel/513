function[B] = HW4_BiotSavart(p, XYZ, I)
N=size(XYZ,1); %how many line segments are there?
u0=4*pi*10^-7; %u0 value
B0=[0,0,0]; %to calculate magnetic field vector
B=0;


if N==2 %Calculates magnetic value if there is only 1 line segment- this is for HW4_Check
    L=[abs(XYZ(2,1)-XYZ(1,1)),abs(XYZ(2,2)-XYZ(1,2)),abs(XYZ(2,3)-XYZ(1,3))]; %to calculate length of line segments
    r0=[((XYZ(2,1)+XYZ(1,1))-p(1))./2,((XYZ(2,2)+XYZ(1,2))-p(2))./2,((XYZ(2,3)+XYZ(1,3))-p(3))./2];
    r=sqrt(sum(r0.^2));
    for i=1:3
        B0(i)=(u0/(2*pi*r))*I./(((r./(L(i)/2)).^2 + 1).^(1/2)); %formula for magnetic field
    end
    B=B+B0;
else    
    for n=1:N
        if n~=N
            L=[abs(XYZ(n+1,1)-XYZ(n,1)),abs(XYZ(n+1,2)-XYZ(n,2)),abs(XYZ(n+1,3)-XYZ(n,3))]; %to calculate length of line segments
            r0=[((XYZ(n+1,1)+XYZ(n,1)))./2,((XYZ(n+1,2)+XYZ(n,2)))./2,((XYZ(n+1,3)+XYZ(n,3)))./2];
        else
            L=[abs(XYZ(N,1)-XYZ(1,1)),abs(XYZ(N,2)-XYZ(1,2)),abs(XYZ(N,3)-XYZ(1,3))]; %to calculate length of line segment at end of array
            r0=[((XYZ(N,1)+XYZ(1,1)))./2,((XYZ(N,2)+XYZ(1,2)))./2,((XYZ(N,3)+XYZ(1,3)))./2]; 
        end
        r=sqrt(sum((r0.^2)));
        for i=1:3
            B0(i)=(u0/(2*pi*r))*I./(((r./(L(i)/2)).^2 + 1).^(1/2)); %formula for magnetic field
        end
    B=B+B0;
    end
end
end

