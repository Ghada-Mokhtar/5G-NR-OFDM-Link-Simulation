# 5G NR OFDM Link-Level Simulation

Link-level simulation of a 5G NR downlink (n78 band, 3.5 GHz, 100 MHz bandwidth), connecting a manual RF link budget to an OFDM waveform simulation and a bit-error-rate (BER) performance analysis under AWGN, Rician (LOS), and Rayleigh (NLOS) fading.

Built as part of a CV project portfolio targeting network/telecom engineering roles (e.g. Ericsson, Huawei).

## What this project does

1. **Link budget** — derives received signal power and SNR for LOS and NLOS conditions by hand, using the Friis equation and the 3GPP TR 38.901 Urban Macro (UMa) propagation model.
2. **OFDM waveform simulation** — two implementations:
   - `ofdm_5g_simulation.m` — built from first principles (manual IFFT/FFT, cyclic prefix, 16-QAM), MATLAB R2014-compatible.
   - `ofdm_5g_simulation_toolbox.m` — rebuilt using MATLAB's 5G Toolbox (`nrCarrierConfig`, `nrPDSCH`, `nrOFDMModulate`), matching real 3GPP numerology (273 resource blocks, 30 kHz SCS).
3. **Channel modeling** — AWGN, Rician (LOS), and Rayleigh (NLOS) flat-fading channels, with BER evaluated across a 0-40 dB SNR sweep.

## Key result

![BER results](results/ber_results.png)

The simulated AWGN curve matches the theoretical Gray-coded 16-QAM reference closely. Rician (LOS) fading degrades performance moderately; Rayleigh (NLOS) fading degrades it substantially and decays much more slowly with SNR — consistent with the ~32 dB path-loss penalty found in the link budget between LOS and NLOS conditions.

## Debugging notes (the most useful part of this repo)

Three real implementation bugs were found and fixed during verification — documented in full in [`report/5G_NR_OFDM_Simulation_Report.pdf`](report/5G_NR_OFDM_Simulation_Report.pdf), Section 6:

1. **Misleading plot scale** — y-axis auto-scaling to the theoretical curve's numerically meaningless tail (~1e-300) compressed the real, measurable BER range into an unreadable sliver.
2. **Es/N0 vs. Eb/N0 mismatch** — the theoretical BER formula expected energy-per-bit, but was fed energy-per-symbol directly, leaving an extra factor of *k* = 4 inside the Q-function and making the reference curve artificially optimistic.
3. **Missing PDSCH descrambling** — `nrPDSCH` applies 3GPP-standard bit scrambling internally; `nrSymbolDemodulate` only reverses modulation, not scrambling, which produced BER ≈ 0.5 at every SNR until `nrPDSCHDecode` (the correct inverse of `nrPDSCH`) was used instead.

## Tools used

MATLAB R2024b, 5G Toolbox, Communications Toolbox.

## Repository structure

```
├── ofdm_5g_simulation.m          # manual OFDM implementation (first principles)
├── ofdm_5g_simulation_toolbox.m  # 5G Toolbox implementation
├── results/
│   └── ber_results.png
└── report/
    └── 5G_NR_OFDM_Simulation_Report.pdf
```

## Author

Ghada Mokhtar — Telecommunications & Network Engineering student, ENSTA (École Nationale Supérieure des Technologies Avancées)
