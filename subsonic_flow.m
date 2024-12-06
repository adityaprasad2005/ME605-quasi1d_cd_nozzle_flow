% Simulation of fluid flow through a CD nozzle for Purely Subsonic Case
clear all;
clc;
close all;


% Parameters
n = 201;                % Odd Number of grid points 
L = 3;                  % Length of the nozzle
x = linspace(0, L, n); % Spatial domain (0 to 3)
dx = L/(n-1);           % Grid spacing
gam = 1.4;           % Specific heat ratio
CFL = 0.6;             % Courant number for stability
max_iter = 10000;              % Max_number of iterations
eps = 1e-6 ;            % convergence criteria

% Initialization of non-dimenstionalized nozzle flow variables and area profile
rho = zeros(1, n); % Density
V = zeros(1, n);   % Velocity
T = zeros(1, n);   % Temperature
P = zeros(1, n);   % Pressure
A = 1 + 2.2 * (x - 1.5).^2;  % Converging-diverging area profile


% Define the initial conditions along the nozzle
for i = 1:n
    if x(i) <= 0.5
        rho(i) = 1;       % Initial density before x = 0.5
        T(i) = 1;         % Initial temperature
    elseif (x(i) > 0.5 && x(i) <= 1.5)
        rho(i) = 1 - 0.363 * (x(i) - 0.5); % Linear decrease in density
        T(i) = 1 - 0.167 * (x(i) - 0.5);   % Linear decrease in temperature
    elseif (x(i) > 1.5 && x(i) <= 3)
        rho(i) = 0.637 +  0.363* (x(i) - 1.5); % Further decrease beyond x = 1.5
        T(i) = 0.835 - 0.165 * (x(i) - 1.5);   % Further decrease in temperature
    else
        rho(i) = 1.0 ;
        T(i) = 1.0 ;
      
    end
    V(i) = 0.59 / (rho(i) * A(i));      % Velocity
end

% Stagnation values
rho(1) = 1 ;
T(1) = 1; 

% initial pressure 
P = rho.*T;
P_out = 0.995;     % Outlet non-dimensional pressure
P(n) = P_out;


% Identifying the throat location (minimum area)
throat_idx = ceil(n/2);

% Defining the initial values of the Conservative variables
U1 = rho .* A;                            
U2 = rho .* A .* V;                        
U3 = rho .* A .* (T / (gam - 1) + 0.5*gam* V.^2);

% Defining the initial values of the Flux variables
F1 = U2;
F2 = ((gam-1)/gam) .*U3  +  ((3-gam)/2) * ((U2.^2)./U1);
F3 = gam*(U2.*U3./U1) + 0.5*(gam*(1-gam))*((U2.^3)./(U1.^2)); 


%% Starting the time-marching approach

% predictor derivative term initialisation to avoid repeted initialisations
dU1_dt_pred = zeros(1, n);
dU2_dt_pred = zeros(1, n);
dU3_dt_pred = zeros(1, n);

% predictor conservative term initialisation
U1_pred = zeros(1, n);
U2_pred = zeros(1, n);
U3_pred = zeros(1, n);


% corrector derivative term initialisation 
dU1_dt_corr = zeros(1, n);
dU2_dt_corr = zeros(1, n);
dU3_dt_corr = zeros(1, n);

time = 0 ;
for iter = 1:max_iter

    % Time step calculation using the grid values from the previous step
   if min(T) < 0
       fprintf("The simulation broke at iter %d due to T < 0", iter)
       break
   else
        dt = min( (CFL*dx)./( T.^0.5 + V ) );
   end
    
    % Predictor Step
    for i = 1:n-1
        dU1_dt_pred(i) = -( F1(i+1)-F1(i) )/dx;
        dU2_dt_pred(i) = -( F2(i+1)-F2(i) )/dx + (1/gam) * ( P(i))*((A(i+1) - A(i))/dx );
        dU3_dt_pred(i) = -( F3(i+1)-F3(i) )/dx;

        U1_pred(i) = U1(i) + dt*dU1_dt_pred(i);
        U2_pred(i) = U2(i) + dt*dU2_dt_pred(i);
        U3_pred(i) = U3(i) + dt*dU3_dt_pred(i);
    end
    
    % Predictor values calculation
    F1_pred = U2_pred;
    F2_pred = ((gam-1)/gam) .*U3_pred  +  ((3-gam)/2) * ((U2_pred.^2)./U1_pred);
    F3_pred = gam*(U2_pred.*U3_pred./U1_pred) + 0.5*(gam*(1-gam))*((U2_pred.^3)./(U1_pred.^2));  


       
    % Corrector Step
    for i = 2:n
        dU1_dt_corr(i) = -(F1_pred(i)-F1_pred(i-1))/dx;
        dU2_dt_corr(i) = -(F2_pred(i)-F2_pred(i-1))/dx + (1/gam)*(P(i))*((A(i) - A(i-1))/dx);
        dU3_dt_corr(i) = -(F3_pred(i)-F3_pred(i-1))/dx;

    end
    
    % Average derivative of conservative variables
    avg_dU1_dt = 0.5.*(dU1_dt_corr + dU1_dt_pred);
    avg_dU2_dt = 0.5.*(dU2_dt_corr + dU2_dt_pred);
    avg_dU3_dt = 0.5.*(dU3_dt_corr + dU3_dt_pred);

    % new conservative variables 
    U1_new = U1 + dt*avg_dU1_dt;
    U2_new = U2 + dt*avg_dU2_dt;
    U3_new = U3 + dt*avg_dU3_dt;

    % inlet boundary conditions such that rho,T,P are new at atleast at their first values
    U1_new(1) = rho(1) * A(1);
    U2_new(1) = 2*U2_new(2) - U2_new(3);
    V(1) = U2_new(1) ./ U1_new(1);
    U3_new(1) = U1_new(1) * ( T(1) / (gam - 1) + 0.5*gam* V(1)^2);  

    % outlet boundary conditions
    U1_new(n) = 2*U1_new(n-1) - U1_new(n-2);
    U2_new(n) = 2*U2_new(n-1) - U2_new(n-2);
    V(n) = U2_new(n)/U1_new(n);
    U3_new(n) = P(n)*A(n)/(gam-1) + 0.5*gam*U2_new(n)*V(n);

    % new Fluz terms calculation
    F1 = U2_new;
    F2 = (U2_new.^2)./(U1_new) + ((gam-1)/gam)*(U3_new - 0.5*gam*(U2_new.^2)./(U1_new));
    F3 = gam*(U2_new.*U3_new./U1_new) + 0.5*(gam*(1-gam))*((U2_new.^3)./(U1_new.^2));  

    % new primitives variables calculation
    rho = U1_new ./ A;
    V = U2_new ./ U1_new;
    T = (gam - 1) * (U3_new ./ U1_new - 0.5* gam * V.^2);
    P = rho .* T;
    M = V ./ sqrt(T); 

    % Updating U
    U1 = U1_new;
    U2 = U2_new;
    U3 = U3_new;

    % convergence check 
    check1 = abs(max(avg_dU1_dt)) < eps ;
    check2 = abs(max(avg_dU2_dt)) < eps ;
    check3 = abs(max(avg_dU3_dt)) < eps ;
    if check1 && check2 && check3
        fprintf("The simulation was successfully completed at iter : %d", iter);
        break
    end

    time = time + dt;
end


%% Plot results 

figure();
plot(x, P, 'LineWidth', 1.5);
title('Pressure Distribution along the Nozzle');
xlabel('Nozzle X-Distance');
ylabel('Non-dimensional Pressure');
grid on;

figure();
plot(x, rho, 'LineWidth', 1.5);
title('Density Distribution along the Nozzle');
xlabel('Nozzle X-Distance');
ylabel('Non-dimensional Density');
grid on;

figure();
plot(x, T, 'LineWidth', 1.5);
title('Temperature Distribution along the Nozzle');
xlabel('Nozzle X-Distance');
ylabel('Non-dimensional Temperature');
grid on;

figure();
plot(x, M, 'LineWidth', 1.5);
title('Mach Number along the Nozzle');
xlabel('Nozzle X-Distance');
ylabel('Mach Number');
grid on;


