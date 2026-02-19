# Active CRLH Transmission Line Modeler

A MATLAB simulation of an **active Composite Right/Left-Handed (CRLH) transmission line** with a frequency-dependent negative resistance element. The structure achieves controlled gain across a user-defined passband by embedding an active shunt element (modeled as a negative resistance `Rn`) within a periodic CRLH unit cell topology.
This is supplemental code for the paper "Negative Resistance Enabled Amplifying CRLH Transmission Lines With Uniform Insertion Gain" DOI: 10.1109/LMWT.2026.3655243 (https://ieeexplore.ieee.org/document/11366944)
---

## Background

CRLH transmission lines support both left-handed (LH) and right-handed (RH) wave propagation, enabling precise dispersion control. By introducing a frequency-shaped negative conductance into the shunt branch, this design compensates for losses and achieves net gain across the passband.

The unit cell uses a **symmetric T-network topology** with:
- **Series branch**: Right-handed inductance `LR` + left-handed capacitance `CL`
- **Shunt branch**: Right-handed capacitance `CR` + left-handed inductance `LL` + active negative resistance `Rn(f)`

---

## Features

- Calculates CRLH element values (`LR`, `CR`, `LL`, `CL`) analytically from two frequency/phase design targets
- Models frequency-dependent negative resistance as `Rn(f) = -A·exp(α·f)`, fit to desired shunt conductance at two frequencies
- Cascades `n` unit cells in a **passive–active–passive** symmetric arrangement
- Converts ABCD matrices to S-parameters
- Plots:
  - Cascaded S-parameters (`S11`, `S21`)
  - Dispersion diagram (`βp` vs. frequency)
  - Shunt conductance vs. frequency
  - Required `Rn(f)` vs. exact analytical solutions
  - Single unit cell S-parameters with approximate insertion loss overlay
  - Bloch impedance (real and imaginary)
- Outputs (in command window):
  - Constituent CRLH parameters
  - $R_n(\omega)$ exponential fit parameters
  - $S_{21}$ at desired frequencies ($\omega_1$ and $\omega_2$)
  - Maximum gain

---

## Dependencies

- [SPARAMS](https://github.com/njchorda/MATLAB-Touchstone-Reader) — provides `SPARAMS.abcd2s()` for ABCD-to-S-parameter conversion. Add as a submodule or clone separately and add to your MATLAB path.
- `closestIdx` — a utility function for finding the nearest index in a vector. Included in the bottom of the main file.

After cloning, add dependencies to your MATLAB path
---

## Usage

1. Open `Active_NRCRLH_Calcs.m` in MATLAB
2. Set your design parameters in the **Input parameters** section:

| Parameter | Description |
|-----------|-------------|
| `f1`, `f2` | Lower and upper frequency bounds of the passband (Hz) |
| `theta1`, `theta2` | Desired phase shifts at `f1` and `f2` (rad) |
| `n` | Number of sub-unit cells (must be **odd** for passive–active–passive symmetry) |
| `Z0` | Reference impedance (default: 50 Ω) |
| `G_desired` | Target shunt conductance (set as a fraction of `G_max`) |

3. Run the script — six figures will be generated automatically. You may need to tune the figX.Position parameter to fit them on your screen

---

## Design Equations

CRLH element values are solved analytically from the two phase/frequency constraints in Reference [1].

Then, the required $R_n(\omega)$ is calculated as:

$$R_n(\omega) = -R_0e^{\alpha\omega}$$

with

$$\alpha = \frac{\ln{\left(\frac{R_n(\omega_1)}{R_n(\omega_2)}\right)}}{\omega_1 - \omega_2}$$,

$$R_0 = -R_n(\omega_1)e^{-\alpha \omega_1}$$.


This will approximately fit the exact solution of:

$$R_n(\omega) = \frac{1 \pm \sqrt{1 - 4(G'_{sh} \omega L_L)^2}}{2G'_{sh}}$$

The required negative resistance at each frequency satisfies:

$$G_{sh} = \frac{R_n}{R_n^2 + (\omega L_L)^2}$$

which is solved exactly at `f1` and `f2`, then fit with an exponential model across the full band.

---

## Output Figures

| Figure | Contents |
|--------|----------|
| 1 | Cascaded `S11` and `S21` (dB) |
| 2 | Dispersion: `βp` (deg) vs. frequency |
| 3 | Shunt conductance `G(f)` vs. frequency |
| 4 | `Rn(f)`: exponential fit vs. exact analytical solutions |
| 5 | Single unit cell S-parameters + approximate `S21` |
| 6 | Bloch impedance `ZB` (real and imaginary) |

---

## Notes

- `n` must be **odd** to maintain the passive–active–passive cascade symmetry
- Loss-balanced condition (`Gsh = Rse/Z0²`) is available but commented out by default
- The Rollett (K-Δ) stability analysis block is present in the code but commented out — uncomment to assess stability across frequency
- Maximum achievable gain is estimated and printed to the console on each run

---

## License

MIT License — see `LICENSE` for details.
