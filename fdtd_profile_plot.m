function fdtd_profile_plot(profile, N, b, varargin)

[epsr,mur,sigma] = fdtd_profile(profile, N, b, varargin);

ms = 1; % Marker (dot) size.
if N < 10
    ms = 20; % Show points if few grid points.
end

clf;
subplot(2,1,1);grid on;hold on;box on;
    stairs(epsr,'b-','Marker','.','MarkerSize',ms,'LineWidth',2);
    stairs(mur,'r-','Marker','.','MarkerSize',ms,'LineWidth',2);
    legend('\epsilon_r','\mu_r');
    set(get(gca,'YLabel'),'Rotation',90,'HorizontalAlignment','Right');
subplot(2,1,2);grid on;hold on;box on;
    stairs(sigma,'k-','Marker','.','MarkerSize',ms,'LineWidth',2);
    axis([1 N min(sigma)*0.9-0.001 max(sigma)*1.1+0.001]);
    ylabel('\sigma [\Omega]');
    set(get(gca,'YLabel'),'Rotation',0,'HorizontalAlignment','Right');
    xlabel('i_x');
drawnow
