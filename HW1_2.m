clc;clear all;

%I calculated "Actual value" with this equation: 
%E_p = K*Q(L/2)/(z*sqrt(z^2 + (L/2)^2)) 
%from this source:
%http://dev.physicslab.org/Document.aspx?doctype=3&filename=Electrostatics_ContinuousChargedRod.xml

fprintf("Running... this takes about a minute to complete\n");
n=50; %number of max iterations for loop
j=2; %number of charges, starting with 2 (smallest number of charges you can have to approximate a line charge)
e=1; %error value
tab = zeros(j-1,2); %for table to be printed
k = 8.987E9; 
while ((e > 0.01) && j<n) %run loop until error reaches < or = 0.01, or for n iterations, to prevent endless loop
    p = ones(j,2); %for location of charges
    for i=1:1:j
       p(i,:) = [2*((i-1)/(j-1))-1,0]; %array of point locations
    end
    Q = ones(j,1); %array of point charges of 1 Coulomb charge
    j=j+1;
    t=0.01;
    [X,Y] = meshgrid(-2:t:2);

    V = zeros(size(X)); 
    for ii = 1:numel(Q) 
        V = V + k * Q(ii) *(p(ii,1)-X)./ (((p(ii,1)-X).^2 + (p(ii,2)-Y).^2).^(1/2)); %calculate electric potential at points
    end
    % This introduces another error - we want to know the error due
    % to discretization and this introduces an error due to the
    % approximated method used to compute the electric field. Should
    % use the exact formula for E for point charges instead. Did you 
    % check your first few calculations by hand? Even given this
    % approxmation, your N seems too high.
    v=-gradient(V); %calculate electric field
    c=(v((1/t)+1, ((4/t)/2)+1)); %electric field at point [0,L], with L=1
    e=(((1.2710E10)-c)/(1.27106E10)); %error
    tab(j-2,:)=[j,e]; %fill table
end
fprintf("Number of point charges: " + string(j)+'\n');
fprintf("Error: " + string(e)+'\n');
fprintf("Electric field value at [0,L] with 1C charge at each point: " + '%d\n',c);
fprintf("Distance between charges: " + 1/j + "L\n");
T = array2table(tab,'VariableNames',{'Num_of_Charges','Error'})
