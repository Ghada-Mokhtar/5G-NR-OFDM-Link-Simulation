%% ========================================================================
%  5G NR OFDM Link Simulation — 5G Toolbox version (MATLAB R2024b)
%
%  This is an upgrade of ofdm_5g_simulation.m: the manual IFFT/CP-based
%  OFDM engine is replaced with genuine 3GPP-standard 5G Toolbox functions
%  (nrCarrierConfig, nrPDSCHConfig, nrPDSCH, nrOFDMModulate/Demodulate).
%  The channel model (flat Rician/Rayleigh fading + perfect equalization)
%  and the theoretical AWGN BER reference are kept unchanged from the
%  validated version, since they do not depend on the OFDM engine itself.
%
%  Link budget context (computed by hand, 3.5 GHz / n78 band, 3GPP UMa,
%  100 MHz bandwidth, NF = 7 dB):
%     LOS  : P_RX = -48.27 dBm -> SNR ~ 38.7 dB
%     NLOS : P_RX = -79.91 dBm -> SNR ~  7.0 dB
% ========================================================================

clear; clc; close all;

%% ---------------------- 1. CARRIER & PDSCH CONFIGURATION ----------------
% n78 band, 100 MHz channel bandwidth, 30 kHz subcarrier spacing (SCS) --
% a common FR1 mid-band numerology choice. NSizeGrid = 273 resource blocks
% is the standard 3GPP value for 100 MHz @ 30 kHz SCS (TS 38.101-1,
% Table 5.3.2-1) -- this is a specified constant, not an assumption.
carrier = nrCarrierConfig;
carrier.SubcarrierSpacing = 30;
carrier.NSizeGrid         = 273;
carrier.CyclicPrefix      = 'normal';

pdsch = nrPDSCHConfig;
pdsch.Modulation = '16QAM';           % same modulation order as before
pdsch.NumLayers  = 1;                 % single layer -- no MIMO spatial
                                       % multiplexing (a real Ericsson gNB
                                       % would typically use multiple
                                       % layers/antennas; omitted here to
                                       % keep the link-level BER analysis
                                       % focused on modulation + channel)
pdsch.PRBSet     = 0:carrier.NSizeGrid-1;   % allocate the full bandwidth

% Resource element indices for the PDSCH allocation, and info about how
% many payload bits (info.G) are needed to fill it for one slot
[pdschIndices, pdschInfo] = nrPDSCHIndices(carrier, pdsch);

% Reference SNR values derived from the manual link budget calculation
SNR_LOS_dB   = 38.7;
SNR_NLOS_dB  = 7.0;

fprintf('=== 5G NR OFDM Simulation (5G Toolbox) ===\n');
fprintf('NSizeGrid: %d RB | SCS: %d kHz | Modulation: %s | Bits/slot: %d\n', ...
    carrier.NSizeGrid, carrier.SubcarrierSpacing, pdsch.Modulation, pdschInfo.G);

%% ---------------------- 2. CHANNEL + BER SWEEP -------------------------
% Same three channel conditions as the validated version:
%   1) AWGN only            -- ideal baseline, matches theory curve
%   2) Rician fading (LOS)  -- dominant direct path + scattered paths
%   3) Rayleigh fading(NLOS)-- no dominant path, purely scattered energy
%
% SIMPLIFYING ASSUMPTIONS (carried over, explicitly flagged):
% - Fading is FLAT across all subcarriers within one OFDM symbol (single
%   complex gain per OFDM symbol, not per-subcarrier frequency-selective
%   fading). A full nrTDLChannel model would be more realistic but adds
%   substantial complexity (delay taps, DMRS-based channel estimation)
%   not attempted in this version.
% - Equalization is PERFECT (receiver divides by the exact known fading
%   gain) -- i.e. genie-aided / perfect channel estimation, not estimated
%   from reference signals (DMRS) as a real receiver would do.
SNR_range_dB   = 0:2:40;
K_factor_dB    = 10;                       % assumed typical LOS Rician K-factor
K_lin          = 10^(K_factor_dB/10);

channelNames = {'AWGN', 'Rician (LOS)', 'Rayleigh (NLOS)'};
BER_sim_all  = zeros(3, length(SNR_range_dB));

% Number of OFDM symbols per slot (grid columns) -- used to generate one
% fading gain per OFDM symbol
L = carrier.SymbolsPerSlot;

for chType = 1:3
    for snrIdx = 1:length(SNR_range_dB)
        snr_dB = SNR_range_dB(snrIdx);

        % ---------------- TRANSMITTER ----------------
        txBits = randi([0 1], pdschInfo.G, 1);
        pdschSymbols = nrPDSCH(carrier, pdsch, txBits);  % modulated, NR-normalized symbols

        txGrid = nrResourceGrid(carrier, pdsch.NumLayers);
        txGrid(pdschIndices) = pdschSymbols;

        % Per-OFDM-symbol flat fading gain
        switch chType
            case 1  % AWGN only
                h = ones(1, L);
            case 2  % Rician (LOS)
                h = sqrt(K_lin/(K_lin+1)) + ...
                    sqrt(1/(2*(K_lin+1))) * (randn(1,L) + 1i*randn(1,L));
            case 3  % Rayleigh (NLOS)
                h = sqrt(1/2) * (randn(1,L) + 1i*randn(1,L));
        end

        fadedGrid = txGrid .* reshape(h, 1, L, 1);   % broadcast across subcarriers & ports

        [txWaveform, ofdmInfo] = nrOFDMModulate(carrier, fadedGrid);

        % ---------------- CHANNEL (AWGN) ----------------
        rxWaveform = awgn(txWaveform, snr_dB, 'measured');

        % ---------------- RECEIVER ----------------
        rxGrid = nrOFDMDemodulate(carrier, rxWaveform);

        % Perfect equalization: divide out the known fading gain per symbol.
        eqGrid = rxGrid ./ reshape(h, 1, L, 1);

        rxPDSCHSymbols = eqGrid(pdschIndices);

        % Noise variance estimate for the decoder (nominal, based on
        % target SNR and unit-average symbol power -- NR modulation
        % mapping functions normalize to unit average power internally)
        noiseVar = 1 / (10^(snr_dB/10));

        % nrPDSCHDecode is the correct inverse counterpart of nrPDSCH: it
        % performs layer de-mapping, LLR (soft-bit) computation, AND
        % descrambling together (nrSymbolDemodulate alone does NOT
        % reverse the scrambling nrPDSCH applies internally -- see
        % Section 6.3 of the project report for the full diagnosis).
        cwLLRs = nrPDSCHDecode(carrier, pdsch, rxPDSCHSymbols, noiseVar);
        rxBits = double(cwLLRs{1} < 0);

        ber = mean(rxBits ~= txBits);
        BER_sim_all(chType, snrIdx) = ber;
    end
end

%% ---------------------- 3. THEORETICAL BER (reference) ------------------
M = 16; k_bits = log2(M);
SNR_lin    = 10.^(SNR_range_dB/10);
BER_theory = (4/k_bits) * (1 - 1/sqrt(M)) * qfunc(sqrt(3*SNR_lin/(M-1)));

%% ---------------------- 4. RESULTS --------------------------------------
figure('Name', '5G NR OFDM Link Performance (5G Toolbox)');
semilogy(SNR_range_dB, BER_sim_all(1,:), 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
hold on;
semilogy(SNR_range_dB, BER_sim_all(2,:), 'gs-', 'LineWidth', 1.5, 'MarkerFaceColor', 'g');
semilogy(SNR_range_dB, BER_sim_all(3,:), 'md-', 'LineWidth', 1.5, 'MarkerFaceColor', 'm');
semilogy(SNR_range_dB, BER_theory, 'r--', 'LineWidth', 1.5);

ylim([1e-6 1]);
yLimits = ylim;
plot([SNR_NLOS_dB SNR_NLOS_dB], yLimits, 'k:', 'LineWidth', 1.2);
plot([SNR_LOS_dB SNR_LOS_dB], yLimits, 'k:', 'LineWidth', 1.2);
text(SNR_NLOS_dB, yLimits(2)*0.5, ' NLOS point (~7 dB)', 'Rotation', 90);
text(SNR_LOS_dB, yLimits(2)*0.5, ' LOS point (~38.7 dB)', 'Rotation', 90);

grid on; xlabel('SNR (dB)'); ylabel('Bit Error Rate (BER)');
title(sprintf('%s, %d RB, %d kHz SCS — AWGN vs Fading (5G Toolbox)', ...
    pdsch.Modulation, carrier.NSizeGrid, carrier.SubcarrierSpacing));
legend('Simulated: AWGN only', 'Simulated: Rician (LOS)', 'Simulated: Rayleigh (NLOS)', ...
    'Theoretical BER (AWGN)', 'Location', 'southwest');
hold off;

fprintf('\n--- BER at link-budget reference points (nearest simulated SNR) ---\n');
[~, idxLOS]  = min(abs(SNR_range_dB - SNR_LOS_dB));
[~, idxNLOS] = min(abs(SNR_range_dB - SNR_NLOS_dB));
fprintf('Near LOS  (%.1f dB) -- AWGN: %.2e | Rician: %.2e\n', ...
    SNR_range_dB(idxLOS), BER_sim_all(1,idxLOS), BER_sim_all(2,idxLOS));
fprintf('Near NLOS (%.1f dB) -- AWGN: %.2e | Rayleigh: %.2e\n', ...
    SNR_range_dB(idxNLOS), BER_sim_all(1,idxNLOS), BER_sim_all(3,idxNLOS));
