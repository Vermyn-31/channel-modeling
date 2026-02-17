%% PRE-PROCESSING OF DATA IN LOW FREQUENCY BAND (Cycles 1 to 13) FOR R1

function process_rx_lf_r1()
% ========================================================================
%  MAIN SCRIPT 
% ========================================================================

% Relative file path of the current script
code_folder = fileparts(mfilename('fullpath'));

% Load the Rx data from SFI Smart Ocean
fpath = fileparts(code_folder);
fpath = fullfile(fpath, 'SFI Smart Ocean Dataset', 'RX-LF', 'R1');

% Center Frequency and Bandwidth
fc = 6000;
BW = 4500;

% Modulation types mapping
modulation_types = {
    'BFGN',       16.00,   0, 96000;
    'A-PRBS',     16.021999999999998,  0, 90000;
    'P-PRBS',     16.162666666666667, 0, 90000;
    'LFM-225',    15.975000000000000,  0, 90000;
    'LFM-3000',   15.000000000000000,   0, 90000;
    'MULTITONE',  60.000000000000000,   0, 96000;
    'OFDM',       7.591666666666667,  0, 96000;
    'OSDM',       3.749500000000000,  0, 96000;
    'MFSK',       4.00,    0, 96000;
    'S2C-QPSK',   2.971437500000000,  0, 96000;
    'BCSK',       2.964291666666667,  0, 96000;
    'FH-BCSK',    2.964291666666667,  0, 96000;
};

% Save path
save_fpath = fullfile(code_folder, 'RX','R1');

% Process cycles 01 to 13
cycle_numbers = [1, 3, 6, 7, 8, 9, 10, 11, 12, 13];
for idx = 1:length(cycle_numbers)
    % Load the rx_data
    cycle_num = cycle_numbers(idx)
    filename = sprintf('LF_R1_cycle%02d.wav', cycle_num);
    [rx_data, fs_rx] = audioread(fullfile(fpath, filename));

    % Low energy detection 
    low_energy_segments = detect_low_energy(rx_data, fs_rx, 0.01, 0.1);
        fprintf('Cycle %02d: Detected %d low-energy segments\n', cycle_num, size(low_energy_segments, 1));

    % Downconvert the signal    
    rx_data = rx_data/max(abs(rx_data));
    t_rx = (0:length(rx_data)-1)' / fs_rx;
    rx_data_downconverted = rx_data .* (sqrt(2) * exp(-1j * 2*pi*fc*t_rx));

    % Apply LPF
    nyquist = fs_rx/2;
    cutoff = BW / nyquist;
    b = fir1(100, cutoff, 'low');
    rx_data_downconverted = filter(b, 1, rx_data_downconverted);


    segment_mapping = [1,2,3,4,5,10,11,12,13,14,15,16];
    
    % Fill-in the start time
    for i = 1:length(segment_mapping)
        idx = segment_mapping(i);
        modulation_types{i, 3} = t_rx(low_energy_segments(idx, 2));
    end

    % Iterate to all modulation_types
    for i = 1:size(modulation_types, 1)

        mod_name = modulation_types{i, 1};
        duration = modulation_types{i, 2};
        start_time = modulation_types{i, 3};
        rx_fs = modulation_types{i, 4};

        % Extract from downconverted signal using actual file sampling rate
        start_idx = round(start_time * fs_rx);
        end_idx = round((start_time + duration) * fs_rx);
        rx_data_final = rx_data_downconverted(start_idx:end_idx);

        % Save the file
        mod_folder = fullfile(save_fpath, mod_name);
        filename = sprintf('cycle%02d.mat', cycle_num);
        save(fullfile(mod_folder, filename), 'rx_data_final', 'rx_fs');

    end
end

fprintf('Done\n');
end % END OF MAIN FUNCTION


% ========================================================================
%  Function: Detect Low-Energy Segments 
% ========================================================================
function low_energy_segments = detect_low_energy(rx_data, fs_rx, low_energy_threshold, min_duration)

    % Calculate energy from signal
    energy = abs(rx_data).^2;

    low_energy_segments = [];
    in_low_segment = false;
    segment_start = 0;

    for i = 1:length(energy)
        if energy(i) < low_energy_threshold
            if ~in_low_segment
                segment_start = i;
                in_low_segment = true;
            end
        else
            if in_low_segment
                segment_end = i - 1;
                segment_duration = (segment_end - segment_start) / fs_rx;

                if segment_duration >= min_duration
                    low_energy_segments = [low_energy_segments; segment_start, segment_end];
                end

                in_low_segment = false;
            end
        end
    end

    if in_low_segment
        segment_end = length(energy);
        segment_duration = (segment_end - segment_start) / fs_rx;

        if segment_duration >= min_duration
            low_energy_segments = [low_energy_segments; segment_start, segment_end];
        end
    end
end