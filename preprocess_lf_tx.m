%% PRE-PROCESSING OF THE TX DATA

% Relative file path of the current script
code_folder = fileparts(mfilename('fullpath'));

% Load the Tx data from SFI Smart Ocean
tx_fpath = fileparts(code_folder);
tx_fpath = fullfile(tx_fpath, 'SFI Smart Ocean Dataset', 'TX', 'wav', 'LF');

% Save the data on this folder
save_fpath = fullfile(code_folder,'TX');

% Center Frequency, Bandwidth and Target Sampling Rate
fc = 6000;
BW = 4500;
target_fs = 96000;

% Load the tx_data 
files_load = {'BFGN-LF.wav', 'A-PRBS-LF-16s.wav', 'P-PRBS-LF.wav', 'LFM-225-LF.wav', 'LFM-3000-LF.wav', 'MULTITONE-LF.wav', 'OFDM-LF.wav', 'OSDM-LF.wav', 'MFSK-LF.wav', 'S2C-QPSK-LF.wav', 'BCSK-LF.wav', 'FH-BCSK-LF.wav'};

for i = 1:12
    filename = files_load{i};
    [tx_data, fs_tx] = audioread(fullfile(tx_fpath, filename));
    
    % Store original data for comparison
    tx_data_original = tx_data;
    fs_original = fs_tx;
    
    % Extract modulation name
    mod_name = strsplit(filename, '-LF');
    mod_name = mod_name{1};
    
    % Resample to match receiver sampling rate
    if fs_tx ~= target_fs
        tx_data = resample(tx_data, target_fs, fs_tx);
        fs_tx = target_fs;
    end
    
    time_length = length(tx_data) / fs_tx;
    fprintf('Time length of %s: %.15f seconds\n', mod_name, time_length);
    fprintf('Number of samples: %d\n', length(tx_data));
    
    t_tx = (0:length(tx_data)-1)' / fs_tx;
    
    % Fix syntax error in this line (added missing closing parenthesis)
    tx_data_downconverted = tx_data .* (sqrt(2)*exp((-1j*2*pi*fc*t_tx)/fs_tx));
    
    % Add zero at the end
    tx_data_downconverted = [tx_data_downconverted; 0];
    
    nyquist = fs_tx/2;
    normalized_cutoff = BW / nyquist;
    filter_order = 100;
    b = fir1(filter_order, normalized_cutoff, 'low');
    
    tx_data_final = filter(b, 1, tx_data_downconverted);

    % Normalize amplitude
    tx_data_final = tx_data_final / max(abs(tx_data_final));

    % Save the file
    save_filename = sprintf('%s.mat', mod_name);
    save(fullfile(save_fpath, save_filename), 'tx_data_final', 'fs_tx');
    
    % Graph original vs resampled (downsampled) signal
    figure;
    subplot(2,1,1);
    plot((0:length(tx_data_original)-1)/fs_original, real(tx_data_original));
    title(sprintf('Original Signal: %s (fs = %d Hz)', mod_name, fs_original));
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
    subplot(2,1,2);
    plot((0:length(tx_data)-1)/fs_tx, real(tx_data));
    title(sprintf('Resampled Signal: %s (fs = %d Hz)', mod_name, fs_tx));
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    
end
