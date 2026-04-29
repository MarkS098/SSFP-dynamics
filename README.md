# SSFP Chemical Exchange Analysis Suite

A MATLAB-based computational framework for extracting multi-site chemical exchange kinetics and transverse relaxation rates from Steady-State Free Precession (SSFP) NMR data.

## ⚙️ Requirements

- MATLAB (R2021a or later recommended)
- Optimization Toolbox (for `lsqnonlin`, `MultiStart`)

## Framework

This software implements a Liouville-space Bloch-McConnell model to simulate and fit the periodic steady state magnetization signal. 

### 1. Objective Function Formulation & Global Joint Fitting
Parameter extraction is treated as a bounded, non-linear least squares optimization problem.

Experimental datasets from all observable sites (e.g., Peak A and Peak B) are concatenated into a single global vector $\mathbf{M}\_{exp}$, which is compared to a global simulated signal $\mathbf{M}\_{sim}(\theta)$. 

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

### 4. Visualize the Fit and Confidence Intervals
The following figures for the fit and $\chi^{2}$ maps will be output after running the fit:
<p float="left">
  <img src="https://github.com/user-attachments/assets/6e3e9e69-8c82-4583-be21-52b4c4c6ac32", width="45%"/>
  <img src="https://github.com/user-attachments/assets/5b6d2cbc-e6df-4fac-991e-7d2a359f93cc", width="45%"/> 
</p>
<p float="left">
   <img src="https://github.com/user-attachments/assets/106e4a98-0521-4155-af07-50e9eeb6bbc9", width="90%"/>
</p>
<p float="left">
   <img src="https://github.com/user-attachments/assets/42f44326-2aca-4854-8b14-4daedcefc5a3", width="43.7%"/>
</p>

## 📂 Script Descriptions

| Script | Purpose |
| :--- | :--- |
| `raw_data_process.m` | Processes raw Bruker data and outputs a .mat file containing the peak intensities, $TR$ values and flip angle. |
| `ssfp_exchange_jointfit.m` | **Main Engine**: Executes global optimization and bootstrap error analysis. |
| `chem_exchange_sim.m` | **Forward Model**: steady-state simulator for $3$ sites. |
| `read_ssfp_acqus.m` | Parses Bruker `acqus` files for $TR$, flip angle, and other experimental parameters. |
| `read_ssfp_procs.m` | Extracts processing parameters (scaling factor, SI etc...). |
| `read_bruker_data.m` | Low-level binary reader for Bruker `fid` and `ser` files. |

## 🚀 Getting Started

1. **Clone the Repo**: 
   ```bash
   git clone [https://github.com/your-username/SSFP-Exchange-Analysis.git](https://github.com/your-username/SSFP-Exchange-Analysis.git)
2. **Process Raw Data**: Open `raw_data_process.m`, point  `data_dir` to your Bruker experiments, select appropriate boundaries for the region of interest, noise region and minimum peak threshold and run it.
3. **Run The Fit**: Open `ssfp_exchange_jointfit.m`, load your processed .mat files, select either a 2 site or 3 site process, appropriate boundaries and execute the solver.
4. **Output**: Logarithmic $\chi^2$ maps for $k_{ex} vs\\,{R_2}$, $k_{ex} vs\\,{\nu_i}$, fitted curves for the SSFP profile and a table of the optimized parameters.
