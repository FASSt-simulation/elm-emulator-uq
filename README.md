[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18645784.svg)](https://doi.org/10.5281/zenodo.18645784)

# UQ Framework for E3SM Land Model (ELM) — Code to Reproduce Figures

This repository contains scripts and notebooks used to train Gaussian-process (GP) emulators, perform Bayesian calibration, and reproduce figures for the manuscript:

**“A framework for parametric and predictive uncertainty quantification in the E3SM Land Model: Assessing site and observable generalizability”**  
(Journal of Advances in Modeling Earth Systems)

## Repository Structure

- `*.ipynb`: Jupyter notebooks for emulator training, calibration, post-processing, and figure generation
- `emulation_functions.jl`: Julia helper functions used by selected workflows
- `fig2_ENF_5sites.ncl`: NCL script to generate the site-location/topography map

> Note: Large intermediate artifacts (e.g., trained emulator binaries `*.jld2` and large ensemble outputs `*.nc`) are not tracked in GitHub.
---

## Figure-to-Script Map (Main Text & Supporting Information)


### Maps / site overview figure
- `fig2_ENF_5sites.ncl`  
  → Site locations + topography map; reproduce **Fig. 2**


### Emulator training & evaluation
- `ESM_emulator_train_evaluation_GPP_26param_fig3.ipynb`  
  → Train and evaluate emulator (26-parameter case); reproduce **Fig. 3**

- `ESM_emulator_train_GPP_4par_US_Me2.ipynb`  
  → Example workflow for training a 4-parameter emulator (illustrated for US-Me2 and GPP); the same script can be adapted to other sites and observables.


### Global sensitivity analysis (GSA)
- `sobol_GSA_heatmap_fig4_5.ipynb`  
  → Sobol-based GSA post-processing; reproduce **Figs. 4–5**


### Bayesian calibration / parameter estimation (PE)
- `PE_3sites_4par_be_fig6_7_S2_S5.ipynb`  
  → Posterior estimation using **synthetic (best-estimate) data** and **FLUXNET data** across sites; reproduce **Figs. 6–7** and **SI Figs. S2 & S5**

- `PE_3sites_4par_fluxnet_fig9_S7.ipynb`  
  → Probabilistic prediction using calibrated parameters across sites; reproduce **Fig. 9** and **SI Fig. S7**

- `PE_4par_be_4qoi_Me2_1be_fig8_S3_S6.ipynb`  
  → Posterior estimation using **synthetic (best-estimate) data** and **FLUXNET data** across observables for specific site; reproduce **Fig. 8** and **SI Figs. S3 & S6**

- `PE_4par_be_4qoi_Me2_figS1.ipynb`  
  → Posterior estimation using **synthetic (best-estimate) data** across observables for specific site; reproduce **SI Fig. S1**

- `PE_4par_fluxnet_4qoi_Me2_fig10_S4_S8.ipynb`  
  → Probabilistic prediction using calibrated parameters across observables for specific site; reproduce **Fig. 10** and **SI Figs. S4 & S8**

---


