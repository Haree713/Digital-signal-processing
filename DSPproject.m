% Filter Specifications
normalized_cutoff = 3/48;  % Normalized cutoff frequency νc
suppression_freq = 6/48;   % Frequency for suppression
filter_order = 54;         % Maximum filter order
target_suppression = 40;   % Desired suppression in dB

% Design the filter with the fir1 function
fir_order = min(filter_order, round(3.3 / (normalized_cutoff))); % Choosing filter order
fir_coeffs = fir1(fir_order, normalized_cutoff, 'low', hamming(fir_order + 1));

% Frequency Response
[H, omega] = freqz(fir_coeffs, 1);

% Plot frequency response
figure;
plot(omega/pi, 20*log10(abs(H)));
title('Magnitude Frequency Response of Type I FIR Filter');
xlabel('Normalized Frequency (\nu)');
ylabel('Magnitude (dB)');
axis([0 0.5 -80 5]); % Set limits for better visualization
grid on;

% Impulse Response
figure;
stem(0:length(fir_coeffs)-1, fir_coeffs);
title('Impulse Response of Type I FIR Filter');
xlabel('Sample');
ylabel('Amplitude');
grid on;


% Given low-pass filter coefficients fir_coeffs from previous code

% Calculate the delay 'm' for constructing the high-pass filter
m = round((length(fir_coeffs) - 1) / 2);

% Construct the high-pass filter's impulse response hH[n]
hH = -fir_coeffs;                      % Invert the low-pass filter's response
hH(m + 1) = hH(m + 1) + 2 * normalized_cutoff;  % Adjust center amplitude to 1/8

% Normalize the high-pass filter coefficients
max_coeff_hH = max(abs(hH));
hH_normalized = hH / max_coeff_hH;

% Frequency Response of the high-pass filter
[H_HighPass, omega_HighPass] = freqz(hH_normalized, 1);

% Plot the magnitude frequency response of the high-pass filter
figure;
plot(omega_HighPass/pi, 20*log10(abs(H_HighPass)));
title('Magnitude Frequency Response of High-Pass FIR Filter');
xlabel('Normalized Frequency (\nu)');
ylabel('Magnitude (dB)');
axis([0 0.5 -80 5]); % Set limits for better visualization
grid on;

% Given low-pass filter coefficients hL[n] (obtained from earlier design)
% ... (Your code to generate the hL[n] coefficients)

% Filter Specifications
normalized_cutoff = 3/48;  % Normalized cutoff frequency νc
filter_order = 54;         % Maximum filter order
F = 10;                    % Chosen F value (adjustable, but not exceeding 16)

% Quantize the filter coefficients based on the specified F value
scaling_factor = 2^F;  % Scaling factor for quantization
hQ = round(fir_coeffs * scaling_factor) / scaling_factor;  % Quantize coefficients

% Frequency Response of the fixed-point filter
[H_FixedPoint, omega_FixedPoint] = freqz(hQ, 1);

% Plot the magnitude frequency response of the fixed-point filter
figure;
plot(omega_FixedPoint/pi, 20*log10(abs(H_FixedPoint)));
title('Magnitude Frequency Response of Fixed-Point Filter');
xlabel('Normalized Frequency (\nu)');
ylabel('Magnitude (dB)');
axis([0 0.5 -80 5]); % Set limits for better visualization
grid on;



% Your code to obtain x_L[n] using the designed low-pass filter

% Assuming x_L[n] is stored in the variable 'x_L'

% Quantize the signal x_L[n] to obtain xq_L[n] using fixed-point implementation
% Your code for quantization (xq_L) based on the previous fixed-point implementation

% Calculate the power of the signal x_L[n]
power_x_L = mean(fir_coeffs.^2);

% Calculate the quantization noise
quantization_error = fir_coeffs - hQ;
power_quantization_noise = mean(quantization_error.^2);

% Calculate SQNR in linear scale
SQNR_linear = power_x_L / power_quantization_noise;

% Convert SQNR to dB
SQNR_dB = 10 * log10(SQNR_linear);

% Display the SQNR in dB
disp(['SQNR (dB): ', num2str(SQNR_dB)]);




