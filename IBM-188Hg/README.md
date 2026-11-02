# IBM-188Hg

Verification and fitting codes for the paper

> M. J. S. Al Musawi, S. M. Mutsher, A. A. Abd-Aljbar, A. A. Elbndag, M. A. Al-Jubbori, H. A. H. Al-Naffakh, F. I. Sharrad,
> *Shape Coexistence and Collective Excitations in ¹⁸⁸Hg within the IBM-1, IBM-2 and Configuration-Mixing IBM Frameworks*.

The repository contains:

- an independent m-scheme diagonaliser for IBM-1, IBM-2 and the configuration-mixing IBM (IBM-CM, [N]⊕[N+2]);
- the fit data set, with the uncertainties actually used;
- the scripts that produce every fit, robustness test and derived number in Tables 2–10 and Figs. 3–5 of the paper;
- the numerical results (JSON).

## Requirements

Python ≥ 3.9 with `numpy`, `scipy`, `sympy` and `matplotlib`:

```
pip install numpy scipy sympy matplotlib
```

Run all scripts from inside the `code/` directory. Several scripts read earlier results (e.g. `fit_kpboth.json`, `cm_fit14_main.json`); copy the files of `results/` into `code/` first (`cp ../results/* .`). The IBM-CM space for ¹⁸⁸Hg has 9191 states; one diagonalisation takes about a second, and a full fit takes a few minutes on a laptop.

## Conventions

- Bosons: s and d. The quadrupole operator is Q(χ) = [d†s̃ + s†d̃]⁽²⁾ + χ[d†d̃]⁽²⁾, and L̂·L̂ = L(L+1).
- IBM-CM (García-Ramos & Heyde, PRC 89, 014306 (2014)):
  - H = H_ecqf^N + (H_ecqf^{N+2} + Δ) + V_mix;
  - V_mix = w₀(s†s† + ss) + w₂((d†d†)⁽⁰⁾ + h.c.);
  - T(E2) = Σᵢ eᵢ Q(χᵢ), with the eᵢ given in √W.u.
- Energies are in keV and B(E2) values in W.u. (1 W.u. = 6.397×10⁻³ e²b² for A = 188).
- Each experimental level is matched to the calculated state of the same configuration character: regular weight > 0.5 (R) or < 0.5 (I).

## Files

### `code/`

| File | Purpose |
|---|---|
| `bos.py` | s–d boson algebra: creation operators, tensor couplings, Q(χ), L, n_d, Clebsch–Gordan coefficients |
| `ibm1.py`, `run2.py` | IBM-1 (multipole form, PHINT conventions); reproduces the O(6) energies of Table 7 |
| `ibm2.py`, `q2.py`, `t4.py`, `t5.py` | IBM-2 (N_π = 1, N_ν = 9): energies, B(E2), F-spin content |
| `ecqf.py` | least-squares ECQF IBM-1 fit to band A (Tables 2 and 4) |
| `ibmcm.py` | IBM-CM Hamiltonian, E2 operator, configuration weights |
| `c14.py` | original restricted IBM-CM fit and its curvature errors (Table 3b, restricted row) |
| `cmfit.py` | general IBM-CM fit driver used for the robustness tests of Table 6 and Fig. 5 (see below) |
| `grid_polish.py` | second pass on the (Δ, w) grid from an alternative starting basin |
| `w2scan.py` | scan of w₂ with w₀ re-optimised (shows the w₀ + 0.445 w₂ degeneracy) |
| `vmix.py` | mixing matrix elements between unperturbed regular and intruder states (Table 9) |
| `ratios.py` | mixing-sensitive branching ratios (code validation against Ref. [16]; unfitted tests in Table 10) |
| `grh2.py`, `grh3.py` | IBM-CM with the published parameters of Ref. [16], and sensitivity of the ratios to rounding and to the V_mix normalisation |
| `post.py` | χ² decomposition, mixing contents, regular-2⁺ prediction and E0 estimate for every fit |
| `su3chk.py` | numerical check of the SU(3) Q·Q eigenvalues (coefficient 3/8 of L(L+1)) |
| `c7c.py` | auxiliary analysis of the main IBM-CM solution |
| `extended.py`, `ext2.py` | full spectrum, B(E2) values and configuration weights of a fit (Tables 7, 9, 10) |
| `pes_ext.py`, `c7c_coh.py` | IBM-CM potential-energy surfaces from a fit (Fig. 4b) |
| `figs3.py`, `fig5.py` | figures 1–5 |

Configurations of `cmfit.py`:

```
python3 cmfit.py kpboth         # FINAL fit: kappa'_N and kappa'_{N+2} free (Tables 3b, 5)
python3 cmfit.py kpboth_w20     # final Hamiltonian, w fixed at 20 keV
python3 cmfit.py kpboth_prolate # final Hamiltonian, chi_{N+2} = -0.17
python3 cmfit.py main_check     # pulls of the restricted fit (Table 5)
python3 cmfit.py prolate        # restricted, chi_{N+2} = -0.17
python3 cmfit.py w20            # restricted, w fixed at 20 keV
python3 cmfit.py w0w2           # w0 and w2 free
python3 cmfit.py kpI            # kappa'_{N+2} free
python3 cmfit.py kpI_prolate2   # chi_{N+2} = -0.17 fixed, kappa'_{N+2} free
python3 cmfit.py relE           # relative energy uncertainties
python3 cmfit.py nolimits       # without the two one-sided limits
python3 cmfit.py grid_0         # (Delta, w) profile, Delta = 3000, 4000 keV
python3 cmfit.py grid_1         # (Delta, w) profile, Delta = 3480, 5000 keV
python3 grid_polish.py 0; python3 grid_polish.py 1   # second pass on the grid
python3 post.py                 # summary of all fits -> post.json
python3 extended.py fit_kpboth.json                  # final-fit spectrum and B(E2)
python3 pes_ext.py fit_kpboth.json                   # final-fit energy surface
```

### `data/`

`fit_dataset_188Hg.csv` contains all fitted observables, with their values, the uncertainties used, their weights and their sources.

### `results/`

JSON outputs of every fit (`fit_*.json`; `fit_kpboth.json` is the final fit, `fit_main_check.json` the restricted fit) and of the polished (Δ, w) grid (`grid_0p.json`, `grid_1p.json`), the post-processing summary (`post.json`), the mixing matrix elements (`vmix.json`) and the w₂ scan (`w2scan.json`). `cm_fit14_main.json` and `cm_results14.json` hold the main fit and its full spectrum.

## Experimental data

- ENSDF: F. G. Kondev, S. Juutinen, D. J. Hartley, Nucl. Data Sheets 150, 1 (2018).
- Lifetimes: M. Siciliano et al., Phys. Rev. C 102, 014318 (2020).
- Coulomb excitation: K. Wrzosek-Lipska et al., Eur. Phys. J. A 55, 130 (2019).

## Licence

MIT. If you use these codes, please cite the paper.

## Additional robustness runs (response to referees, round 2)

```
bash runspecs.sh        # 20 fits listed in specs.txt (cmfit.py spec mode):
                        #  F_eI (eps_{N+2} free), F_D / R_D (Delta free), F_cN0 (chi_N = 0),
                        #  F_w*/R_w* (w profile at Delta = 3480 keV), F_D*/R_D* (Delta profile)
python3 chain.py R      # Delta profile, restricted Hamiltonian, warm-started chain (chain_R.json)
python3 chain.py F      # Delta profile, final Hamiltonian (chain_F.json)
python3 lowD.py         # final Hamiltonian at Delta = 2.5 and 3.0 MeV from the Delta-free solution
python3 w2scan_final.py # w2 scan with w0 re-optimised, final Hamiltonian
python3 fig5.py         # Fig. 5 (chi^2 profiles in w and Delta)
NB=81 BMAX=2 NG=25 OUT=pes_fine_kpboth.npz python3 pes_ext.py fit_kpboth.json   # fine-grid surface for Fig. 4b
```
`F_*` = final Hamiltonian (kappa'_N, kappa'_{N+2} free); `R_*` = restricted Hamiltonian. `dprofile.json` collects the converged Delta profile with the ground-state intruder admixture.
