clc;clear all;
k = 8.987E9;
d = 1;
p = [d,0; -d,0]; % location of points
Q = [1; -1]; % 1 and -1 Coulomb charge
t=0.05;
[X,Y] = meshgrid(-2:t:2);

V = zeros(size(X)); 
for ii = 1:numel(Q) 
    V = V + k * Q(ii)./ (((p(ii,1)-X).^2 + (p(ii,2)-Y).^2).^(1/2)); %calculate electric potential at points
end

colormap(parula(15)) % It does not makes sense for colorbar to have so
% many colors when plot essentially shows 3. Should also plot lot(V), but
% negative values would need special treatment.
Vo = (k*Q(1)/d^2); % I'd use dimensionless variables
contourf(X,Y,V/Vo)
hColorbar = colorbar;
ylabel(hColorbar,'Electric potential V/V_o')
title('Field Lines of Point Charges 1C and -1C')

[Ex,Ey] = gradient(V);
C = all(isfinite(Ex) & isfinite(Ey));
% You found it!
field = streamslice(X(:,C),Y(:,C),-Ex(:,C),-Ey(:,C));
xlabel('x/d');
ylabel('y/d');
set(field,'Color','k');
axis square
print -dpdf -bestfit HW1_1.pdf