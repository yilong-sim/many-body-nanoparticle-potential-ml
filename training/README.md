# training

Extracts two- and three-body PMFs from LAMMPS force output files and
prepares training data for PIP fitting.

## File

- `extract_pmf.py` — integrates constraint forces → W2(d), ΔW3(d12,d13,d23)

## Usage

```bash
python extract_pmf.py --input-dir ./NPforced_data/ --output ./pmf_results/
```

## Input

`NPforced{N}.dat` files from `simulation/2_pmf_calculation/`.
Each file contains ensemble-averaged constraint forces at one grid position
along the reaction coordinate.

## Output

```
pmf_results/
  W2.csv    # two-body PMF: d12 (σ), W2 (kBT)
  W2.png    # PMF plot with minimum labeled
```

## Method

Numerical integration of ensemble-averaged constraint forces:

    W(d) = -∫ <F(ξ)> dξ

The three-body contribution is isolated as:

    ΔW3 = W3'(d12, d13, d23) - W2(d13) - W2(d23)

Reference: Sprik & Ciccotti, J. Chem. Phys. 109, 7737 (1998)

## Next step

Feed W2.csv and ΔW3 data into `fitting/mbfit/NP_Dimer_fitting.ipynb`
and `NP_Trimer_fitting.ipynb` for PIP fitting.
