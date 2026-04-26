# SSFP Chemical Exchange Analysis Suite

A MATLAB-based computational framework for analyzing multi-site chemical exchange in Steady-State Free Precession (SSFP) NMR experiments. This package addresses the mathematical degeneracy between exchange rates ($k_{ex}$) and transverse relaxation ($R_2$) through a **Global Joint Fitting** architecture.

## 🔬 Scientific Framework

This software implements a Liouville-space Bloch-McConnell model to simulate the periodic steady state of magnetization. Key algorithmic features include:

### 1. Global Joint Fitting
To resolve parameter correlations, experimental data from multiple sites (e.g., Peak A and Peak B) are concatenated and fitted simultaneously. A single **analytical scale factor** $c(\theta)$ is computed via orthogonal projection at each iteration:
$$c(\theta) = \frac{\mathbf{M}_{sim}(\theta)^T \mathbf{M}_{exp}}{\mathbf{M}_{sim}(\theta)^T \mathbf{M}_{sim}(\theta)}$$
This ensures that fractional populations are constrained by the relative intensities of the peaks, preventing unphysical scaling.

### 2. MultiStart & Latin Hypercube Sampling (LHS)
The $\chi^2$ error surface is often non-convex. We utilize LHS to provide a space-filling distribution of $N_{start}$ initial guesses, which are then refined using a trust-region-reflective local solver to guarantee identification of the global minimum.

### 3. Segregated Residual Bootstrap
Uncertainty is quantified using a non-parametric bootstrap protocol. To preserve site-specific noise characteristics, residuals are isolated per peak and resampled independently to generate synthetic datasets:
$$\mathbf{M}_{syn,i}^{(b)} = c(\hat{\theta})\mathbf{M}_{sim,i}(\hat{\theta}) + \mathbf{e}_i^{*(b)}$$

## 📂 Script Descriptions

| Script | Purpose |
| :--- | :--- |
| `raw_data_process.m` | Extracts peak intensities from raw Bruker data and handles frequency calibration. |
| `ssfp_exchange_jointfit.m` | **Main Engine**: Executes global joint optimization and bootstrap error analysis. |
| `chem_exchange_sim.m` | **Forward Model**: Matrix exponential-based steady-state simulator for $N$ sites. |
| `read_ssfp_acqus.m` | Parses Bruker `acqus` files for $T_R$, flip angle, and carrier frequency. |
| `read_ssfp_procs.m` | Extracts processing parameters (scaling factors, SI). |
| `read_bruker_data.m` | Low-level binary reader for Bruker `fid` and `ser` files. |

## 🚀 Getting Started

1. **Clone the Repo**: 
   ```bash
   git clone [https://github.com/your-username/SSFP-Exchange-Analysis.git](https://github.com/your-username/SSFP-Exchange-Analysis.git)
