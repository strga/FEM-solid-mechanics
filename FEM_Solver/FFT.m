clear all;
close all;
clc;
format longG;

fft_files = {
    'Beam_5_FFT.dat'
};

index = 14;                % spectrum truncation
chunkCols = 512;    
useSingle = false;         

combined_frequencies = cell(length(fft_files), 1);
combined_fft_sums    = cell(length(fft_files), 1);

for file_idx = 1:length(fft_files)
    data = load(fft_files{file_idx});
    time = data(:, 1);
    displacements = data(:, 2:end);
    if useSingle, displacements = single(displacements); end

    N  = length(time);
    dt = time(2) - time(1);
    Fs = 1 / dt;
    delta_F = Fs / N;

    max_index  = floor(N / index);
    freq_vec   = (0:N-1) * (Fs / N);
    frequencies = freq_vec(1:max_index);

    nDOF = size(displacements, 2);

    % Accumulate the global spectrum without storing the whole FFT
    if useSingle
        global_fft_sum = zeros(max_index, 1, 'single');
    else
        global_fft_sum = zeros(max_index, 1);
    end

    for cStart = 1:chunkCols:nDOF
        cEnd = min(cStart + chunkCols - 1, nDOF);
        X = fft(displacements(:, cStart:cEnd), [], 1);     % N x (#chunk)
        X = X / N;                                          % Normalize by number of time steps
        X(2:end-1, :) = 2 * X(2:end-1, :);                  % One-sided spectrum scaling
        mag = abs(X(1:max_index, :));                      % Amplitude spectrum (in meters)
        global_fft_sum = global_fft_sum + sum(mag, 2);     % Sum across DOFs

        clear X mag;
    end

    % Detect peaks on the accumulated spectrum
    [peaks, peak_freqs] = findpeaks(double(global_fft_sum), double(frequencies), ...
                                    'MinPeakDistance', delta_F * 5);

    % -------- PLOTS --------
    figure;
    plot(frequencies, double(global_fft_sum), 'LineWidth', 1.5);
    hold on; plot(peak_freqs, peaks, 'ro');
    xlabel('Frequency (Hz)'); ylabel('Summed FFT Magnitude');
    title(['Global Frequency Spectrum - File ', num2str(file_idx)]);
    grid on; legend('FFT Magnitude', 'Detected Peaks');

    figure;
    semilogy(frequencies, double(global_fft_sum), 'LineWidth', 1.5);
    hold on; semilogy(peak_freqs, peaks, 'ro');
    xlabel('Frequency (Hz)'); ylabel('Summed FFT Magnitude (Log Scale)');
    title(['Logarithmic FFT Spectrum - File ', num2str(file_idx)]);
    grid on; legend('FFT Magnitude', 'Detected Peaks');

    fprintf('\n--- FFT File %d: %s ---\n', file_idx, fft_files{file_idx});
    fprintf('Detected dominant frequencies (Hz):\n');
    disp(peak_freqs');

    combined_frequencies{file_idx} = frequencies;
    combined_fft_sums{file_idx}    = double(global_fft_sum);   % store as double for plotting
end

% ==========================
% Combined Semilog Plot
% ==========================
figure; hold on;
colors = ['b', 'r', 'g', 'k', 'm', 'c'];
for file_idx = 1:length(fft_files)
    semilogy(combined_frequencies{file_idx}, combined_fft_sums{file_idx}, ...
        'Color', colors(mod(file_idx-1, numel(colors))+1), 'LineWidth', 1.5);
end
xlabel('Frequency (Hz)');
ylabel('Summed FFT Magnitude (Log Scale)');
title('Combined Logarithmic FFT Spectrum from Multiple Files');
legend(fft_files, 'Interpreter', 'none');
grid on;
