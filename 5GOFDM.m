%% ========================================================================
%  5G NR-style OFDM Link Simulation
%  Compatible with MATLAB R2014 and later (uses base MATLAB + Communications
%  Toolbox functions only: qammod, qamdemod, awgn, fft, ifft)
%
%  This script builds an OFDM waveform from first principles (manual
%  IFFT/FFT, manual cyclic prefix insertion) rather than using higher-level
%  5G Toolbox objects, so every step is transparent and reproducible.
%
%  Link budget context (computed by hand beforehand, 3.5 GHz / n78 band,
%  3GPP UMa scenario, 100 MHz bandwidth, NF = 7 dB):
%     LOS  : P_RX = -48.27 dBm -> SNR ~ 38.7 dB
%     NLOS : P_RX = -79.91 dBm -> SNR ~  7.0 dB
%  These two SNR values are used below as realistic reference points.
% ========================================================================

clear; clc; close all;

%% ---------------------- 1. SYSTEM PARAMETERS ---------------------------
N_fft        = 256;          % FFT size (scaled down from a real 5G NR
                              % ~2048/4096-point FFT for fast simulation;
                              % the underlying principle is identical)
N_usedSC     = 200;           % Number of active (data-carrying) subcarriers
                              % out of N_fft (rest are guard bands, as in
                              % real OFDM systems to avoid spectral leakage)
CP_ratio     = 1/4;           % Cyclic prefix length as a fraction of N_fft
                              % (typical for 5G NR "normal CP")
N_cp         = round(N_fft * CP_ratio);
M            = 16;            % QAM modulation order (16-QAM, common for
                              % mid-range SNR in real 5G deployments)
bitsPerSym   = log2(M);
N_symbols    = 500;           % Number of OFDM symbols to simulate

% Reference SNR values derived from the manual link budget calculation
SNR_LOS_dB   = 38.7;
SNR_NLOS_dB  = 7.0;

fprintf('=== 5G NR-style OFDM Simulation ===\n');
fprintf('FFT size: %d | Active subcarriers: %d | CP length: %d | Modulation: %d-QAM\n', ...
    N_fft, N_usedSC, N_cp, M);

%% ---------------------- 2. TRANSMITTER ---------------------------------
% Generate random bit stream for all OFDM symbols
totalBits = N_usedSC * bitsPerSym * N_symbols;
txBits    = randi([0 1], totalBits, 1);

% QAM modulation: map bits -> complex symbols
% qammod expects integer symbols, so bits are grouped and converted first
bitsMatrix   = reshape(txBits, bitsPerSym, []).';      % each row = bitsPerSym bits
symInts      = bi2de(bitsMatrix, 'left-msb');          % convert to integers 0..M-1
qamRaw       = qammod(symInts, M);                     % raw (unnormalized) QAM symbols

% Manual power normalization (R2014-compatible: no 'UnitAveragePower' option
% existed yet in that release, so we normalize by hand).
% For a standard square M-QAM constellation (levels +-1,+-3,...), the
% average symbol energy is E_avg = 2*(M-1)/3. Dividing by sqrt(E_avg)
% gives unit average power, which is required for SNR to be well-defined.
E_avg        = 2*(M-1)/3;
qamSymbols   = qamRaw / sqrt(E_avg);

% Reshape QAM symbols into subcarriers x OFDM symbols
qamGrid      = reshape(qamSymbols, N_usedSC, N_symbols);

% Map active subcarriers onto the full FFT grid (zero-pad guard bands)
% Active subcarriers are centered around DC, mimicking real resource grids
freqGrid = zeros(N_fft, N_symbols);
firstSC  = floor((N_fft - N_usedSC)/2) + 1;
freqGrid(firstSC:firstSC+N_usedSC-1, :) = qamGrid;

% IFFT: convert each OFDM symbol from frequency domain to time domain
% This is the core OFDM operation -- it's what makes many parallel
% narrowband subcarriers combine into one time-domain waveform
timeGrid = ifft(freqGrid, N_fft, 1) * sqrt(N_fft); % sqrt(N_fft) preserves energy

% Add cyclic prefix: copy the last N_cp samples of each symbol to its front
% This combats inter-symbol interference (ISI) from multipath delay spread
txWaveform = [timeGrid(end-N_cp+1:end, :); timeGrid];
txSerial   = txWaveform(:);   % serialize all OFDM symbols into one stream

fprintf('Transmit waveform generated: %d time-domain samples\n', length(txSerial));

%% ---------------------- 3. CHANNEL + BER SWEEP -------------------------
% We sweep SNR across a range AND mark our two computed link-budget points
SNR_range_dB = 0:2:40;
BER_sim      = zeros(size(SNR_range_dB));

for k = 1:length(SNR_range_dB)
    snr_dB = SNR_range_dB(k);

    % Pass through AWGN channel at this SNR
    rxSerial = awgn(txSerial, snr_dB, 'measured');

    % ---------------- RECEIVER ----------------
    rxWaveform = reshape(rxSerial, N_fft + N_cp, N_symbols);
    rxNoCp     = rxWaveform(N_cp+1:end, :);          % remove cyclic prefix
    rxFreqGrid = fft(rxNoCp, N_fft, 1) / sqrt(N_fft); % back to frequency domain

    rxQamGrid  = rxFreqGrid(firstSC:firstSC+N_usedSC-1, :); % extract active SCs
    rxSymbols  = rxQamGrid(:);

    rxSymbolsRaw = rxSymbols * sqrt(E_avg);  % undo normalization before demod
    rxInts     = qamdemod(rxSymbolsRaw, M);
    rxBitsMat  = de2bi(rxInts, bitsPerSym, 'left-msb');
    rxBits     = reshape(rxBitsMat.', [], 1);

    % Bit error rate for this SNR point
    [~, ber]      = biterr(txBits, rxBits);
    BER_sim(k)    = ber;
end

%% ---------------------- 4. THEORETICAL BER (reference) ------------------
% Theoretical BER for M-QAM over AWGN (Gray-coded), for comparison against
% our simulated curve -- validates that the simulation is implemented correctly
SNR_lin   = 10.^(SNR_range_dB/10);
k_bits    = log2(M);
BER_theory = (4/k_bits) * (1 - 1/sqrt(M)) * qfunc(sqrt(3*k_bits*SNR_lin/(M-1)));

%% ---------------------- 5. RESULTS --------------------------------------
figure('Name', '5G NR-style OFDM Link Performance');
semilogy(SNR_range_dB, BER_sim, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
semilogy(SNR_range_dB, BER_theory, 'r--', 'LineWidth', 1.5);

% Mark our two link-budget-derived operating points
yLimits = ylim;
plot([SNR_NLOS_dB SNR_NLOS_dB], yLimits, 'k:', 'LineWidth', 1.2);
plot([SNR_LOS_dB SNR_LOS_dB], yLimits, 'k:', 'LineWidth', 1.2);
text(SNR_NLOS_dB, yLimits(2)*0.5, ' NLOS point (~7 dB)', 'Rotation', 90);
text(SNR_LOS_dB, yLimits(2)*0.5, ' LOS point (~38.7 dB)', 'Rotation', 90);

grid on; xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title(sprintf('%d-QAM OFDM over AWGN — Simulated vs Theoretical BER', M));
legend('Simulated BER', 'Theoretical BER (AWGN)', 'Location', 'southwest');
hold off;

fprintf('\n--- BER at link-budget reference points (nearest simulated SNR) ---\n');
[~, idxLOS]  = min(abs(SNR_range_dB - SNR_LOS_dB));
[~, idxNLOS] = min(abs(SNR_range_dB - SNR_NLOS_dB));
fprintf('Near LOS  (%.1f dB): BER ~ %.2e\n', SNR_range_dB(idxLOS), BER_sim(idxLOS));
fprintf('Near NLOS (%.1f dB): BER ~ %.2e\n', SNR_range_dB(idxNLOS), BER_sim(idxNLOS));