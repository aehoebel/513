clc;clear all;
k = 8.987E9; 
p = [1,0; -1,0]; % location of points
Q = [1; -1]; % 1 and -1 Coulomb charge
t=0.05;
[X,Y] = meshgrid(-2:t:2);

V = zeros(size(X)); 
for ii = 1:numel(Q) 
    V = V + k * Q(ii)./ (((p(ii,1)-X).^2 + (p(ii,2)-Y).^2).^(1/2)); %calculate electric potential at points
end

contourf(X,Y,V)
hColorbar = colorbar;
ylabel(hColorbar,'Electric potential (V)')
title('Field Lines of Point Charges 1C and -1C')

[Ex,Ey] = gradient(V);
C = all(isfinite(Ex) & isfinite(Ey));
field = streamslice(X(:,C),Y(:,C),-Ex(:,C),-Ey(:,C));
set(field,'Color','k');
