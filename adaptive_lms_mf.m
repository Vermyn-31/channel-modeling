%% ADAPTIVE ALGORITHM IMPLEMENTATION USING LEAST MEAN SQUARE (LSM)

% Rx and Tx relative path
script_folder = fileparts(mfilename('fullpath'));

rx_path = fullfile(script_folder, 'RX', 'R2','OFDM');
tx_path = fullfile(script_folder, 'TX');

% Load Tx data 
tx_struct = load(fullfile(tx_path, 'OFDM.mat'));
tx_data = tx_struct.tx_data_final;

% Modulation scheme: OFDM
modulation_scheme = 'OFDM';

% ARMAX parameters: na=2, nb=2, nc=2, nk=1
na_armax = 2;  
nb_armax = 2;  
nc_armax = 2;
nk_armax = 1;  

% Total number of parameters
M_armax = na_armax + nb_armax + nc_armax;

% LMS parameters
mu_armax = 0.01;      

% Initialize parameter vectors for real and imaginary parts
theta_armax_real = zeros(M_armax, 1); 
theta_armax_imag = zeros(M_armax, 1); 

% Store first cycle parameters
first_theta_real = [];
first_theta_imag = [];

fprintf('Modulation scheme: %s\n', modulation_scheme);
fprintf('Tx data length: %d samples\n\n', length(tx_data));

% Initialize storage variables
first_fit_real = 0;
first_fit_imag = 0;
last_fit_real = 0;
last_fit_imag = 0;

first_y_test_real = [];
first_y_pred_real = [];
first_y_test_imag = [];
first_y_pred_imag = [];
last_y_test_real = [];
last_y_pred_real = [];
last_y_test_imag = [];
last_y_pred_imag = [];

% Process multiple cycles
for cycle_num = 1:10
    % Format cycle filename
    cycle_file = sprintf('cycle%02d.mat',cycle_num);

    fprintf('=== Processing Cycle %d: %s ===\n', cycle_num, cycle_file);
    
    % Load Rx data for current cycle
    rx_struct = load(fullfile(rx_path, cycle_file));
    rx_data = rx_struct.rx_data_final;
    rx_fs = rx_struct.rx_fs;
    
    fprintf('Rx data length: %d samples\n', length(rx_data));
    fprintf('Sampling frequency: %.0f Hz\n', rx_fs);
    
    % Create iddata objects for real and imaginary parts
    z_real = iddata(real(rx_data), real(tx_data), 1/rx_fs);
    z_imag = iddata(imag(rx_data), imag(tx_data), 1/rx_fs);
    
    % Get total number of samples
    total_samples = length(z_real.OutputData);
    fprintf('Total samples: %d\n', total_samples);
    
    % 70/30 data split for both real and imaginary parts
    train_samples = round(0.7 * total_samples);
    
    % Real part
    z_train_real = z_real(1:train_samples);
    z_test_real = z_real(train_samples+1:end);
    u_train_real = z_train_real.InputData;  
    y_train_real = z_train_real.OutputData; 
    u_test_real = z_test_real.InputData;    
    y_test_real = z_test_real.OutputData; 
    
    % Imaginary part
    z_train_imag = z_imag(1:train_samples);
    z_test_imag = z_imag(train_samples+1:end);
    u_train_imag = z_train_imag.InputData;  
    y_train_imag = z_train_imag.OutputData; 
    u_test_imag = z_test_imag.InputData;    
    y_test_imag = z_test_imag.OutputData; 
    
    fprintf('Training samples: %d (70%%)\n', train_samples);
    fprintf('Testing samples: %d (30%%)\n\n', total_samples - train_samples);
    
    % Initialize vectors for past values (real part)
    y_past_armax_real = zeros(na_armax, 1);
    u_past_armax_real = zeros(nb_armax, 1);
    e_past_armax_real = zeros(nc_armax, 1);
    
    % Initialize vectors for past values (imaginary part)
    y_past_armax_imag = zeros(na_armax, 1);
    u_past_armax_imag = zeros(nb_armax, 1);
    e_past_armax_imag = zeros(nc_armax, 1);
    
    %% LMS Adaptive ARMAX Model for Real Part
    fprintf('=== LMS ADAPTIVE ARMAX MODEL (REAL PART) ===\n');
    
    % LMS adaptive filtering for ARMAX (real part)
    for n = 1:length(y_train_real)
        % Create regression vector: φ = [-y(n-1), -y(n-2), u(n-1), u(n-2), e(n-1), e(n-2)]
        phi_armax_real = [-y_past_armax_real; u_past_armax_real; e_past_armax_real];
        
        % Compute prediction
        y_hat_armax_real = phi_armax_real' * theta_armax_real;
        
        % Compute error
        e_armax_real = y_train_real(n) - y_hat_armax_real;
        
        % Update parameters using LMS
        theta_armax_real = theta_armax_real + mu_armax * e_armax_real * phi_armax_real;
        
        % Update past values
        y_past_armax_real = [y_train_real(n); y_past_armax_real(1:end-1)];
        u_past_armax_real = [u_train_real(n); u_past_armax_real(1:end-1)];
        e_past_armax_real = [e_armax_real; e_past_armax_real(1:end-1)];
    end
    
    % Test the adaptive ARMAX model (real part)
    y_past_test_armax_real = zeros(na_armax, 1);
    u_past_test_armax_real = zeros(nb_armax, 1);
    e_past_test_armax_real = zeros(nc_armax, 1);
    y_pred_adaptive_armax_real = zeros(size(y_test_real));
    
    for n = 1:length(y_test_real)
        % Create regression vector
        phi_test_armax_real = [-y_past_test_armax_real; u_past_test_armax_real; e_past_test_armax_real];
        
        % Compute prediction
        y_pred_adaptive_armax_real(n) = phi_test_armax_real' * theta_armax_real;
        
        % Compute error for MA part
        e_test_armax_real = y_test_real(n) - y_pred_adaptive_armax_real(n);
        
        % Update past values for test
        y_past_test_armax_real = [y_test_real(n); y_past_test_armax_real(1:end-1)];
        u_past_test_armax_real = [u_test_real(n); u_past_test_armax_real(1:end-1)];
        e_past_test_armax_real = [e_test_armax_real; e_past_test_armax_real(1:end-1)];
    end
    
    % Calculate fit for adaptive ARMAX (real part)
    fit_adaptive_armax_real = 100 * (1 - norm(y_test_real - y_pred_adaptive_armax_real) / norm(y_test_real - mean(y_test_real)));
    fprintf('Test Fit (Real Part): %.2f%%\n', fit_adaptive_armax_real);
    
    %% LMS Adaptive ARMAX Model for Imaginary Part
    fprintf('=== LMS ADAPTIVE ARMAX MODEL (IMAGINARY PART) ===\n');
    
    % LMS adaptive filtering for ARMAX (imaginary part)
    for n = 1:length(y_train_imag)
        % Create regression vector: φ = [-y(n-1), -y(n-2), u(n-1), u(n-2), e(n-1), e(n-2)]
        phi_armax_imag = [-y_past_armax_imag; u_past_armax_imag; e_past_armax_imag];
        
        % Compute prediction
        y_hat_armax_imag = phi_armax_imag' * theta_armax_imag;
        
        % Compute error
        e_armax_imag = y_train_imag(n) - y_hat_armax_imag;
        
        % Update parameters using LMS
        theta_armax_imag = theta_armax_imag + mu_armax * e_armax_imag * phi_armax_imag;
        
        % Update past values
        y_past_armax_imag = [y_train_imag(n); y_past_armax_imag(1:end-1)];
        u_past_armax_imag = [u_train_imag(n); u_past_armax_imag(1:end-1)];
        e_past_armax_imag = [e_armax_imag; e_past_armax_imag(1:end-1)];
    end
    
    % Test the adaptive ARMAX model (imaginary part)
    y_past_test_armax_imag = zeros(na_armax, 1);
    u_past_test_armax_imag = zeros(nb_armax, 1);
    e_past_test_armax_imag = zeros(nc_armax, 1);
    y_pred_adaptive_armax_imag = zeros(size(y_test_imag));
    
    for n = 1:length(y_test_imag)
        % Create regression vector
        phi_test_armax_imag = [-y_past_test_armax_imag; u_past_test_armax_imag; e_past_test_armax_imag];
        
        % Compute prediction
        y_pred_adaptive_armax_imag(n) = phi_test_armax_imag' * theta_armax_imag;
        
        % Compute error for MA part
        e_test_armax_imag = y_test_imag(n) - y_pred_adaptive_armax_imag(n);
        
        % Update past values for test
        y_past_test_armax_imag = [y_test_imag(n); y_past_test_armax_imag(1:end-1)];
        u_past_test_armax_imag = [u_test_imag(n); u_past_test_armax_imag(1:end-1)];
        e_past_test_armax_imag = [e_test_armax_imag; e_past_test_armax_imag(1:end-1)];
    end
    
    % Calculate fit for adaptive ARMAX (imaginary part)
    fit_adaptive_armax_imag = 100 * (1 - norm(y_test_imag - y_pred_adaptive_armax_imag) / norm(y_test_imag - mean(y_test_imag)));
    fprintf('Test Fit (Imaginary Part): %.2f%%\n', fit_adaptive_armax_imag);
    
    % Store first cycle parameters and data
    if cycle_num == 1
        first_theta_real = theta_armax_real;
        first_theta_imag = theta_armax_imag;
        first_fit_real = fit_adaptive_armax_real;
        first_fit_imag = fit_adaptive_armax_imag;
        first_y_test_real = y_test_real;
        first_y_pred_real = y_pred_adaptive_armax_real;
        first_y_test_imag = y_test_imag;
        first_y_pred_imag = y_pred_adaptive_armax_imag;
    end
    
    % Store last cycle data (always update, will end with cycle 10)
    last_fit_real = fit_adaptive_armax_real;
    last_fit_imag = fit_adaptive_armax_imag;
    last_y_test_real = y_test_real;
    last_y_pred_real = y_pred_adaptive_armax_real;
    last_y_test_imag = y_test_imag;
    last_y_pred_imag = y_pred_adaptive_armax_imag;
    
    fprintf('\n');
end

fprintf('=== Final Model Parameters ===\n');
fprintf('Real Part Parameters: a1=%.4f, a2=%.4f, b1=%.4f, b2=%.4f, c1=%.4f, c2=%.4f\n', ...
        theta_armax_real(1), theta_armax_real(2), theta_armax_real(3), theta_armax_real(4), theta_armax_real(5), theta_armax_real(6));
fprintf('Imag Part Parameters: a1=%.4f, a2=%.4f, b1=%.4f, b2=%.4f, c1=%.4f, c2=%.4f\n\n', ...
        theta_armax_imag(1), theta_armax_imag(2), theta_armax_imag(3), theta_armax_imag(4), theta_armax_imag(5), theta_armax_imag(6));

%% Create equation strings
% First cycle equations
first_eq_real = sprintf('y(t) = %.4f*y(t-1) + %.4f*y(t-2) + %.4f*u(t-1) + %.4f*u(t-2) + %.4f*e(t-1) + %.4f*e(t-2)', ...
    -first_theta_real(1), -first_theta_real(2), first_theta_real(3), first_theta_real(4), first_theta_real(5), first_theta_real(6));

first_eq_imag = sprintf('y(t) = %.4f*y(t-1) + %.4f*y(t-2) + %.4f*u(t-1) + %.4f*u(t-2) + %.4f*e(t-1) + %.4f*e(t-2)', ...
    -first_theta_imag(1), -first_theta_imag(2), first_theta_imag(3), first_theta_imag(4), first_theta_imag(5), first_theta_imag(6));

% Last cycle equations
last_eq_real = sprintf('y(t) = %.4f*y(t-1) + %.4f*y(t-2) + %.4f*e(t-1) + %.4f*e(t-2) + %.4f*u(t-1) + %.4f*u(t-2)', ...
    -theta_armax_real(1), -theta_armax_real(2), theta_armax_real(3), theta_armax_real(4), theta_armax_real(5), theta_armax_real(6));

last_eq_imag = sprintf('y(t) = %.4f*y(t-1) + %.4f*y(t-2) + %.4f*e(t-1) + %.4f*e(t-2) + %.4f*u(t-1) + %.4f*u(t-2)', ...
    -theta_armax_imag(1), -theta_armax_imag(2), theta_armax_imag(3), theta_armax_imag(4), theta_armax_imag(5), theta_armax_imag(6));

%% Plot Last Cycle Results
figure(1);

% Real Part (Last Cycle)
subplot(2,1,1);
plot(1:length(last_y_test_real), last_y_test_real, 'k-', 'LineWidth', 1.5);
hold on;
plot(1:length(last_y_pred_real), last_y_pred_real, 'r--', 'LineWidth', 1);
title_str = sprintf('Cycle 10: Adaptive ARMAX Model - Real Part [%s] (Fit: %.2f%%)', ...
    modulation_scheme, last_fit_real);
title({title_str, last_eq_real});
xlabel('Sample');
ylabel('Amplitude');
legend('Validation Data', 'Prediction', 'Location', 'best');
grid on;

% Imaginary Part (Last Cycle)
subplot(2,1,2);
plot(1:length(last_y_test_imag), last_y_test_imag, 'k-', 'LineWidth', 1.5);
hold on;
plot(1:length(last_y_pred_imag), last_y_pred_imag, 'b--', 'LineWidth', 1);
title_str = sprintf('Cycle 10: Adaptive ARMAX Model - Imaginary Part [%s] (Fit: %.2f%%)', ...
    modulation_scheme, last_fit_imag);
title({title_str, last_eq_imag});
xlabel('Sample');
ylabel('Amplitude');
legend('Validation Data', 'Prediction', 'Location', 'best');
grid on;