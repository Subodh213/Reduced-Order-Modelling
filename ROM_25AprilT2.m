
clc;
clear;

load('C:\ifiss3.6\datafiles\square_unsteadycd.mat', 'U');
Usnap = U; % snapshot Matrix
%%
Umean = mean(Usnap, 2);
Usnap_centered = Usnap - Umean;
%%

%first column%
U_1 = Usnap(:, 1);
%%
%  SVD on the snapshot matrix U_snap
[Uvec, S, V] = svd(Usnap);  % 
singularValues = diag(S);
figure;
semilogy(singularValues, 'o-');
xlabel('Mode number');
ylabel('Singular value (log scale)');
title('Decay of Energy');
grid on;

%%
r = 50; 
V_red = Uvec(:, 1:r);  % Reduced basis matrix 
load('C:\ifiss3.6\datafiles\ref_grid.mat', 'y');
load('C:\ifiss3.6\datafiles\ref_grid.mat', 'x');
load('C:\ifiss3.6\datafiles\ref_grid.mat', 'xy');

% load('C:\ifiss3.6\datafiles\square_cd_nobc.mat', 'A');
load('C:\ifiss3.6\datafiles\square_cd_nobc.mat', 'N');
% load('C:\ifiss3.6\datafiles\square_cd_nobc.mat', 'Q');

load('C:\ifiss3.6\datafiles\square_cd.mat', 'Asupg');
load('C:\ifiss3.6\datafiles\square_cd.mat', 'fsupg');
load('C:\ifiss3.6\datafiles\square_cd.mat', 'Q');
% load('C:\ifiss3.6\datafiles\square_cd_nobc.mat', 'fsupg');
% Q=Mass Matrix
% N=Convection
% A= Convection-Diffusion
% Stiffness Matrix = A-N
 f = fsupg;
 S = Asupg;

 S_r = V_red' * S * V_red;
 C_r =  V_red' * N * V_red;
 M_r =  V_red' * Q * V_red;
 f_r =  V_red' * f;


%%

% disp(S_r);
% disp(C_r);
% disp(M_r);

% Time-stepping setup
delt = 0.1; % Time step size
t_f= 15; % Final time
t = 1:delt:t_f; % Time vector
M= t_f/delt;
num_steps = length(t) ;
% ur = zeros(r, M+1);
u_int = 0.000;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  ; 
ur(:,1) = u_int * ones(r, 1);


for i = 1:M
    
    %  the implicit step
    a = M_r + delt*S_r + delt*C_r; 
    b = (M_r * ur(:,i))+delt*f_r; 

    ur(:,i+1) = a\b;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    % Check for solution divergence
    if any(abs(ur(:,i+1)) > 1e10)
        warning('Solution diverg at timestep %d', i);
        break;
    end
       % Check for solution convergence
  
    %    % Check for solution convergence
    % if (abs(ur(:,i+1)-ur(:,i)) < 1e-5)
    %     warning('Solution converge at timestep %d', i);
    %     tF = i;
    %     break;
    % end
    
end
disp(ur);
% plot(ur);
tvec = 1:delt:t_f; % Time vector
%%

u_full = V_red * ur;

    % t = 100;

square_heatplot(u_full(:,2:end), tvec(:,2:end), xy,x, y, 0.5, 29);


%%
%Interpolation of FOM solution for uniform time steps
load('C:\ifiss3.6\datafiles\Tfine.mat', 'time');
load('C:\ifiss3.6\datafiles\Ufine.mat', 'U');
Ufine=U;
%%
timeFOM = time;
% Define the ROM time steps, which are uniform
t_final = t_f; % final time of 5 seconds
n_time_steps_ROM = M+1; % number of time steps in ROM
timeROM = linspace(0, t_final, n_time_steps_ROM);

% Preallocate the interpolated matrix
UsnapInterpolated = zeros(size(Ufine, 1), n_time_steps_ROM);

% Interpolate each row of Usnap onto the new time grid
for i = 1:size(Ufine, 1)
    UsnapInterpolated(i, :) = interp1(timeFOM, Ufine(i, :), timeROM);
end

% % Calculate the error between FOM and ROM solutions
error = u_full - UsnapInterpolated;

L2_norm = sqrt(trapz((error).^2, 1))/i;  % Summing along the columns


figure;

semilogy(timeROM, L2_norm, '-o','LineWidth',1.5); 
title('L^2 Norm of Error (log-log)');
xlabel('Time');
ylabel('L^2 Norm');
grid on;

% error = Usnap-u_full;
% L2_norm = sqrt(trapz((error).^2, 1));  % Summing along the columns
% semilogy(nonUniformTimeSteps, L2_norm); 
% title('L^2 Norm of Error (log-log)');
% xlabel('Time');
% ylabel('L^2 Norm');
% grid on;

% Calculate L-infinity norm of the error over time
Linfty_norm = max(abs(error), [], 1);  % Taking max along the columns

L1_L2_error = sum((L2_norm ) * delt);
