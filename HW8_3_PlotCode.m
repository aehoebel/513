% The following is a modified version of fdtd_original.m obtained from
% https://www.mathworks.com/matlabcentral/fileexchange/7459-fdtd1d-m
% The additions include the ability to extract the field at certain
% locations and to define different runs.
% To speed up the run, comment out the plotting in the loop.
% The basics of the algorithm when sigma = 0 is described in
% https://my.ece.utah.edu/~ece6340/LECTURES/lecture%2014/FDTD.pdf
% See also https://eecs.wsu.edu/~schneidj/ufdtd/ufdtd.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scott Hudson, WSU Tri-Cities
%1D electromagnetic finite-difference time-domain (FDTD) program.
%Assumes Ey and Hz field components propagating in the x direction.
%Fields, permittivity, permeability, and conductivity
%are functions of x. Try changing the value of "profile".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(0,'DefaultFigureWindowStyle','docked');
set(0,'DefaultAxesFontName','Times');
set(0,'DefaultAxesTitleFontWeight','normal');
set(0,'DefaultTextFontName','Times');
set(0,'DefaultTextFontSize',16);
set(0,'DefaultAxesFontSize',16);

% https://www.mathworks.com/matlabcentral/answers/338733-how-to-stop-legend-from-adding-data1-data2-when-additional-data-is-plotted
set(0,'DefaultLegendAutoUpdate','off');

%close all;
clear all;

animate = 0;

eps0 = 8.854e-12; % permittivity of free space
mu0  = pi*4e-7;   % permeability of free space
run = 1;

if run == 0
    profile = 0; % eps = eps_o, mu = mu_o, sigma = 0.
    source = 2;  % Gaussian pulse at left boundary
    xg = Lx;
    Niter = 600; % # of iterations to perform
end

if run == 1
    profile = 1;
    source  = 1;
    Lx  = 5;       % Domain length in meters
    Nx  = 500;     % Spatial samples in domain
    ixb = Nx/2;
    fs = 300e6;   % Source frequency in Hz
    fstr = '300 MHz';
    Niter = 500;  % Number of iterations to perform
    ip = Nx/2 - 2; % Index of probe
    ylims = [-3, 3];
end

ds = Lx/Nx; % spatial step in meters
dt = ds/fs; % "magic time step"
% See https://my.ece.utah.edu/~simpson/ECE5340/Taflove%20Chpt.%202.pdf
% for definition of magic time step.

% Scale factors for E and H
ae = ones(Nx,1)*dt/(ds*eps0);
am = ones(Nx,1)*dt/(ds*mu0);
as = ones(Nx,1);

% Create grid of epsilon, mu, sigma.
[epsr,mur,sigma] = fdtd_profile(profile, Nx, ixb);

figure(1);
fdtd_profile_plot(profile, Nx, ixb);

ae = ae./epsr;
am = am./mur;
ae = ae./(1+dt*(sigma./epsr)/(2*eps0));
as = (1-dt*(sigma./epsr)/(2*eps0))./(1+dt*(sigma./epsr)/(2*eps0));

% Initialize fields to zero
Hz = zeros(Nx,1);
Hz1 = zeros(Nx,1);
Hz1p = zeros(Nx,1);
Hz1m = zeros(Nx,1);
Ey = zeros(Nx,1);
Ey1 = zeros(Nx,1);
Ey2 = zeros(Nx,1);
Ey1p = zeros(Nx,1);
Ey1m = zeros(Nx,1);

figure(2);clf
    set(gcf,'doublebuffer','on'); % For smoother graphics
    grid on;
    plot(Ey,'b','LineWidth',2);
    hold on;
    plot(377*Hz,'r','LineWidth',2);
    set(gca,'YLim',ylims); 

c = sqrt(1/eps0/mu0);
fprintf('-------------------------\n')
fprintf('Nx = %d\n',Nx);
fprintf('Lx = %.1f [m]\n',Lx);
fprintf('dx = Lx/Nx = %.1e [m]\n',ds);
fprintf('fs = %.1e [Hz]\n',fs);
fprintf('dt = ds/fs = %.1e [s]\n',dt);
fprintf('i_lamda = %.1f\n',Nx*c/fs/Lx);
fprintf('lamda   = %.2e [m]\n',c/fs/Lx);
fprintf('max(sigma)*2*pi*f/epsilon_o = %.1e\n',max(sigma)*2*pi*fs/eps0);
fprintf('-------------------------\n')

a=0;
V0=1;                          %Voltage
I0=0.05;                       %Current
Z0=42;                         %Characteristics Impedance
V01=(V0+I0*Z0);             %Positive travelling Voltage wave
V02=(V0-I0*Z0);             %Negative travelling Voltage wave

Es=0;
m=1;

for iter=1:Niter
    % Source
    if source == 1
        Ey(2) = sin(2*pi*fs*dt*iter);
        Ey1(2) = sin(2*pi*fs*dt*iter);
        Ey2(2) = sin(2*pi*fs*dt*iter);
   %     Ey1p(2) = (V01*exp(-a*iter).*cos(2*pi*fs*dt*iter));
        Ey1p(2) = sin(2*pi*fs*dt*iter);
        Ey1m(2) = (V02*exp(a*iter).*cos(2*pi*fs*dt*iter));
    end
    if source == 2
        % Gaussian pulse
        Ey(3) = exp(-((iter-10)/5)^2);
    end

  %  ae_2 = ones(length(ae)/2,1);
  %  i=length(ae)/2;
  
     ae_2=376.4777;
    
    % The next 10 or so lines of code are where we actually integrate Maxwell's
    % equations. All the rest of the program is basically bookkeeping and plotting.
    Hz(1) = Hz(2); % Absorbing boundary conditions for left-propagating waves
    for i=2:Nx-1 % Update H field
      Hz(i) = Hz(i)-am(i)*(Ey(i+1)-Ey(i));
      Hz1(i) = Hz1(i)-am(i)*(Ey1(i+1)-Ey1(i));
      Hz1p(i) = Hz1p(i)-am(i)*(Ey1p(i+1)-Ey1p(i));
      Hz1m(i) = Hz1m(i)-am(i)*(Ey1m(i+1)-Ey1m(i));
    end
    Ey(Nx) = Ey(Nx-1); % Absorbing boundary conditions for right-propagating waves
    Ey1(Nx) = Ey1(Nx-1);
    Ey2(Nx) = Ey2(Nx-1);
    Ey1p(Nx) = Ey1p(Nx-1);
    Ey1m(Nx) = Ey1m(Nx-1);
    for i=2:Nx-1 % Update E field
      Ey(i) = as(i)*Ey(i)-ae(i)*(Hz(i)-Hz(i-1));
      Ey1(i) = as(i)*Ey1(i)-ae(i)*(Hz1(i)-Hz1(i-1));
      Ey1p(i) = as(i)*Ey1p(i)-ae_2*(Hz1p(i)-Hz1p(i-1));
      Ey1m(i) = Ey1(i)-Ey1p(i);
    end
        
  %  Hz_ip(iter) = Hz(ip);
  %  Ey_ip(iter) = Ey(ip);
    
    Ey_1=zeros(Nx/2,1);
    Ey_1p=zeros(Nx/2,1);
    Ey_1m=zeros(Nx/2,1);
    i=1;
    for n=1:Nx/2
        Ey_1(n)=Ey1(i);
        Ey_1p(n)=Ey1p(i);
        Ey_1m(n)=Ey1m(i);
        i=i+1;
    end
    
    Ey_T=zeros(Nx/2,1);
    i=1;
    for n=1:Nx/2
        Ey_T(n)=Ey_1p(i)+Ey_1m(i);
        i=i+1;
    end
    x1=[1:Nx];
    lambda=100;
    w=2*pi*fs;
    
    if (animate || iter > 425)
            subplot(2,1,1);grid on;hold off;box on;
            plot(Ey_T,'b','LineWidth',2);
            hold on;
            plot(Ey_1p,'g','LineWidth',2)
            hold on;
            plot(Ey_1m,'r','LineWidth',2);
            title(sprintf('i_t = %03d; f = %s [Hz]; L_x = %.1f [m]',...
                iter,fstr,Lx));
            if 1
                legend('E_y_1','E^+_y_1','E^-_y_1','location','northwest');
                set(get(gca,'YLabel'),'Rotation',90,'HorizontalAlignment','Left'); 
                set(gca,'Ylim',[-2,2]);
            else
                legend('E_y [V/m]','377H_z [T]') % Slows down rendering
            end
         %   fdtd1d_annotate;
            print('HW8_2','-dpdf','-fillpage')
            subplot(2,1,2);grid on;hold on;box on;
            title(sprintf('Standing Wave Pattern: |E|_m_a_x/|E|_m_i_n = 1.5/0.5 = 3'));
            Es=2*cos(2*pi*x1/lambda)*cos(w*m*dt)-1.5*cos(2*pi*(x1+m)/lambda);
            plot(Es,'y','LineWidth',1);
            m=m+1;
            set(gca,'Ylim',[-2,2]);
            set(gca,'Xlim',[1,Nx/2]);
            xlabel('i_x');
            if (iter==Niter)
         %       print('HW8_3_2','-dpdf','-fillpage')
            end
            drawnow
            %pause(0);
    end

    if iter == 1
        fprintf('i_t (Time Step) = %04d',iter);
    end
    if iter > 1
        fprintf('\b\b\b\b');
        fprintf('%04d',iter);
    end

end
fprintf('\n');


%fdtd1d_plot;