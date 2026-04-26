# SSFP Chemical Exchange Analysis Package

A MATLAB-based computational suite for extracting multi-site chemical exchange kinetics and transverse relaxation rates from Steady-State Free Precession (SSFP) NMR data. This package utilizes a Liouville-space Bloch-McConnell formalism to perform **Global Joint Fitting** across multiple resonance sites.

## Overview

This package is designed to analyze NMR signal attenuation as a function of repetition time ($T_R$) and off-resonance conditions. It resolves the mathematical degeneracy between chemical exchange ($k_{ex}$) and transverse relaxation ($R_2$) by simultaneously optimizing multiple observation sites and utilizing rigorous non-parametric uncertainty quantification.

### Key Features
- **Liouville-Space Simulation**: Exact numerical solution of the Bloch-McConnell equations using matrix exponential propagators.
- **Global Joint Fitting**: Simultaneous optimization of Peak A and Peak B using a shared parameter set to constrain fractional populations.
- **Analytical Scale Projection**: Automatic hardware gain compensation via orthogonal projection of the simulated model onto experimental data.
- **Global Search Strategy**: Combines **Latin Hypercube Sampling (LHS)** with a **MultiStart** heuristic to navigate non-convex error surfaces.
- **Non-Parametric Bootstrap**: Dataset-segregated resampling to provide stable error estimates where Jacobian-based methods fail.

## Package Structure

### Core Processing Scripts
- **`raw_data_process.m`**: The primary entry point for raw Bruker data. It handles directory traversal, applies FFT, frequency axis calibration (ppm), and extracts peak intensities.
- **`ssfp_exchange_jointfit.m`**: The optimization engine. It executes the MultiStart global search, performs the joint fit, and runs the bootstrap error analysis.
- **`chem_exchange_sim.m`**: The physical forward model. It constructs the 10-dimensional Liouvillian matrix (relaxation, precession, and kinetics) and solves for the periodic steady state.

### Utility Scripts (Bruker I/O)
- **`read_ssfp_acqus.m`**: Parses Bruker `acqus` files for parameters such as Carrier Frequency, Transmitter Offset, and $T_R$.
- **`read_ssfp_procs.m`**: Parses Bruker `procs` files for processing parameters like `SI` and `NC_proc`.
- **`read_bruker_data.m`**: Binary reader for `fid` or `ser` files, compatible with TopSpin 3 (int32) and TopSpin 4 (double) formats.

## Mathematical Formalism

### 1. Objective Function & Joint Fitting
To prevent the solver from artificially scaling individual sites to mask exchange-induced attenuation, a global scale factor $c(\theta)$ is computed:
$$c(\theta) = \frac{\mathbf{M}_{sim}(\theta)^T \mathbf{M}_{exp}}{\mathbf{M}_{sim}(\theta)^T \mathbf{M}_{sim}(\theta)}$$
The objective function returns a dimensionless residual vector $\mathbf{r}(\theta)$, normalized by the mean absolute intensity of the global experimental signal:
$$\mathbf{r}(\theta) = \frac{c(\theta)\mathbf{M}_{sim}(\theta) - \mathbf{M}_{exp}}{\frac{1}{N}\sum_{i=1}^{N}|M_{exp,i}|}$$

### 2. Uncertainty Quantification
Because the Jacobian is often ill-conditioned, we utilize a **Non-Parametric Bootstrap**. The residuals are isolated per peak ($\mathbf{e}_A$, $\mathbf{e}_B$) and resampled with replacement to generate $B$ synthetic datasets:
$$\mathbf{M}_{syn,A}^{(b)} = \hat{M}_{A} + \mathbf{e}_A^{*(b)}, \quad \mathbf{M}_{syn,B}^{(b)} = \hat{M}_{B} + \mathbf{e}_B^{*(b)}$$
The standard deviation of the resulting distribution provides the reported parameter uncertainties ($\sigma_{\theta}$).

## Usage

1. **Extraction**: Run `raw_data_process.m` to extract intensity-vs-TR curves from your Bruker experiment directories.
2. **Optimization**: Run `ssfp_exchange_jointfit.m`. Ensure the `load` paths point to your generated `.mat` files.

1. **Clone the Repo**: 
   ```bash
   git clone [https://github.com/your-username/SSFP-Exchange-Analysis.git](https://github.com/your-username/SSFP-Exchange-Analysis.git)
