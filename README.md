# SSFP Chemical Exchange Analysis Suite

A MATLAB-based computational framework for extracting multi-site chemical exchange kinetics and transverse relaxation rates from Steady-State Free Precession (SSFP) NMR data.

## 🔬 Scientific Framework

This software implements a Liouville-space Bloch-McConnell model to simulate the periodic steady state of magnetization. 

### 1. Objective Function Formulation & Global Joint Fitting
Parameter extraction is treated as a bounded, non-linear least squares optimization problem.

Experimental datasets from all observable sites (e.g., Peak A and Peak B) are concatenated into a single global vector $\mathbf{M}\_{exp}$, which is evaluated against a global simulated signal $\mathbf{M}\_{sim}(\theta)$. 

To account for arbitrary global scaling differences, an analytical scale factor $c(\theta)$ is computed via orthogonal projection at each iteration:

$$c(\theta) = \frac{\mathbf{M}_{sim}(\theta)^T \mathbf{M}_{exp}}{\mathbf{M}_{sim}(\theta)^T \mathbf{M}_{sim}(\theta)}$$

The objective function returns a dimensionless residual vector $\mathbf{r}(\theta)$, normalized by the mean absolute intensity of the global experimental signal:

$$\mathbf{r}(\theta) = \frac{c(\theta)\mathbf{M}_{sim}(\theta) - \mathbf{M}_{exp}}{\frac{1}{N}\sum_{i=1}^{N}|M_{exp,i}|}$$

### 2. Parameter Optimization via LHS and MultiStart
The $\chi^{2}$ error surface for multi-site exchange is often non-convex. To identify the true global minimum, the processing pipeline utilizes a **MultiStart** global search heuristic combined with **Latin Hypercube Sampling (LHS)**. 

* **LHS**: Stratifies each parameter dimension into equally probable intervals to ensure a highly uniform, space-filling distribution of $N_{start}$ initial guesses.
* **MultiStart**: Deploys a local Trust-Region-Reflective solver from each starting coordinate, reliably identifying the global optimum $\hat{\theta}$ even in the presence of rugged error valleys.

### 3. Error Estimation via Non-Parametric Bootstrap 
Standard error estimation derived from the Jacobian matrix was found to be unreliable due to severe ill-conditioning. To rigorously quantify parameter uncertainty, a **Non-Parametric Bootstrap** protocol is utilized:

1. **Dataset-Segregated Resampling**: After identifying the global optimum $\hat{\theta}$, the global residual vector is partitioned back into site-specific error arrays $\mathbf{e}\_{A}$ and $\mathbf{e}\_{B}$. 
2. **Synthetic Joint Data Generation**: New noise vectors $\mathbf{e}\_{A}^{(b)}$ and $\mathbf{e}\_{B}^{(b)}$ are constructed by drawing from their respective site-specific residual pools with replacement. This preserves unequal noise variances native to individual peaks:

$$\mathbf{M}_{syn,A}^{(b)} = \mathbf{M}_{sim,A}(\hat{\theta}) + \mathbf{e}_{A}^{*(b)}$$

$$\mathbf{M}_{syn,B}^{(b)} = \mathbf{M}_{sim,B}(\hat{\theta}) + \mathbf{e}_{B}^{*(b)}$$

3. **Joint Stochastic Probing**: The synthetic datasets are concatenated, and the joint optimization is re-executed for $B$ iterations. The parameter uncertainties are defined as the standard deviation of the resulting bootstrapped distribution.

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
2.
3.
4.
