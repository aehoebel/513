Nx = 550;
x = [1:Nx];
x2 = [Nx:-1:1];
lambda = 100;

figure(1);clf
for i = 1:500
    Ei = cos(2*pi*(x-i)/lambda);
    Ei2 = 0.5*cos(2*pi*(x2-i)/lambda);
    % Set values to right of wave front to NaN so they won't be plotted.
    Ei(i+1:end) = NaN; 
    Ei2(1:end-i) = NaN; 

    % Plot current time step as light grey.
    plot(Ei,'k','LineWidth',1,'Color',[1,1,1,0.4]/2);
    hold on;
    plot(Ei2,'k','LineWidth',1,'Color',[1,1,1,0.4]/2);
    % Keep past time steps
    hold on;
    grid on;
    if i > 1
        % Delete previous current time step thick black line
        delete(h)
        delete(h2)
    end
    % Plot current time step as thick black line
    h = plot(Ei,'k','LineWidth',2);
    h2 = plot(Ei2,'k','LineWidth',2);
    
    set(gca,'Ylim',[-2,2]);
    set(gca,'Xlim',[1,Nx]);
    %legend('V^{+}');
    % Uncomment the following to hide past time steps
    %hold off;
    if mod(i,100) == 0
        % Allow early termination of animation
        input('Continue?');
    end
    drawnow;
end