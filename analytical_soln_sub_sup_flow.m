% Here we find the analytical solution of Subsonic-Supersonic flow in a CD nozzle
clc;
clear all;
close all;


%% plot the nozzle geometry A(x)= 1+ 2.2*(x-1.5)^2  for 0<x<3

L = 3;              % Length of the nozzle
n = 31;             % Number of points (odd so as to get midpoint throat index)
x = linspace(0, L, n); % Discretize the nozzle length
dx = L/(n-1) ;      % Step size  
x_dash = x / L;     % Normalized x-coordinate

A = 1 + 2.2 * (x - 1.5).^2; % Nozzle area function
% y = sqrt(A)./pi;
%  
% figure();
% hold on
% plot(x,0);
% plot(x,y);
% plot(x,-y);
% hold off

%% Analytical solution with A` i.e A/A0 value with x
gam = 1.4;                  % Specific heat ratio

% Find the throat index
idx_throat = ceil(n / 2);   
A_throat = A(idx_throat);    % Throat area

% Initialize M array to store Mach number solutions
M_dash = zeros(1, n);

% Initialize non-dimenstionlized P_dash array
P_dash = zeros(1, n);
pressure = @(M_i) ( 1 + 0.5*(gam-1)*M_i^2 )^( - gam / ( gam -1 )) ;

% Initialize non-dimenstionlized T_dash array
T_dash = zeros(1, n);
temperature= @(M_i) ( 1 + 0.5*(gam-1)*M_i^2 )^( - 1 ) ;

% Initialize non-dimenstionlized Rho_dash array
Rho_dash = zeros(1, n);
density = @(M_i) ( 1 + 0.5*(gam-1)*M_i^2 )^( - 1 / ( gam -1 )) ;

for i = 1:n
    func = @(M_i) A(i) / A_throat - (1 / M_i) * ((2 / (1 + gam)) * (1 + 0.5 * (gam - 1) * M_i^2))^((gam + 1) / (2 * (gam - 1)));
    
    % Use a different initial guess for M depending on the position in the nozzle
    if i < idx_throat
        M_i_guess = 0.1;  % Subsonic guess for converging section
    else
        M_i_guess = 2.0;  % Supersonic guess for diverging section
%         M_i_guess = 0.2;  % Subsonic guess for diverging section (adjust as needed)
    end
    
    % Solve for Mach number at each location
    M_dash(i) = fsolve(func, M_i_guess);

    % Solve for the primitive variables

    P_dash(i) =  pressure(M_dash(i)) ;
    T_dash(i) = temperature(M_dash(i)) ;
    Rho_dash(i) = density(M_dash(i)) ;
end

% % Display results
% figure();
% subplot(4,1,1);
% plot(x, M_dash ,"LineStyle","--");
% title('Mach number M'); 
% 
% subplot(4,1,2);
% plot(x, P_dash, "LineStyle","--");
% title('Non-dimenstionlized P');

% Initial conditions for Rho and T
T_dash_initial = zeros(1,n) ;
Rho_dash_initial = zeros(1,n) ;
for i= 1:n
    if i*dx < 0.5
        T_dash_initial(i) = 1.0 ;
        Rho_dash_initial(i) = 1.0 ;
    elseif i*dx>=0.5 && i*dx <=1.5
        T_dash_initial(i) = 1.0 - 0.167*(i*dx -0.5) ;
        Rho_dash_initial(i) = 1.0 - 0.336*(i*dx -0.5) ;
    else
        T_dash_initial(i) = 0.833 - 0.3507*(i*dx -1.5) ;
        Rho_dash_initial(i) = 0.634 - 0.3879*(i*dx -1.5) ;
    end
end

subplot(2,1,1);
hold on
plot(x, T_dash, "LineStyle","--");
plot(x, T_dash_initial);
ylabel("T'");
title('Non-dimenstionlized Temperature');

subplot(2,1,2);
hold on
plot(x, Rho_dash, "LineStyle","--");
plot(x, Rho_dash_initial);
title('Non-dimenstionlized Density');
ylabel("Rho'")
xlabel("x");



