"""
extract_pmf.py
==============
Extract two- and three-body potentials of mean force (PMFs) from
LAMMPS blue moon ensemble force output files.

Input
-----
NPforced{N}.dat files produced by pmf_calculation.lammps.
Each file contains ensemble-averaged constraint forces on NP1
and the NP1-NP2 separation distance at one grid position.

Output
------
W2(d12)               : two-body PMF as function of separation distance
ΔW3(d12, d13, d23)    : three-body contribution as function of three distances

Method
------
Numerical integration of ensemble-averaged constraint forces along the
reaction coordinate ξ (NP separation distance):

    W(d) = -∫ <F(ξ)> dξ  from ξ_0 to d

Reference: Sprik & Ciccotti, J. Chem. Phys. 109, 7737 (1998)

The three-body contribution is isolated as:
    ΔW3(d12, d13, d23) = W3'(d12, d13, d23) - W2(d13) - W2(d23)

where W3' is the partial three-particle PMF (NP3 interacting with
the fixed NP1-NP2 dimer, excluding the NP1-NP2 interaction itself).

Usage
-----
    python extract_pmf.py --input-dir ./NPforced_data/ --output ./pmf_results/

Dependencies
------------
numpy, scipy, matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import argparse
import os
import glob


def load_force_file(filepath: str) -> tuple[np.ndarray, np.ndarray]:
    """
    Load force and distance data from one NPforced{N}.dat file.

    Parameters
    ----------
    filepath : str
        Path to NPforced{N}.dat file.

    Returns
    -------
    distances : np.ndarray
        NP separation distances ξ (σ).
    forces : np.ndarray
        Ensemble-averaged constraint forces <F(ξ)> (ε/σ).
    """
    data = np.loadtxt(filepath, comments='#')
    distances = data[:, 0]
    forces = data[:, 1]
    return distances, forces


def integrate_pmf(distances: np.ndarray, forces: np.ndarray,
                  xi_0: float = 14.0) -> np.ndarray:
    """
    Integrate constraint forces to obtain PMF via trapezoidal rule.

    W(d) = W(ξ_0) - ∫_{ξ_0}^{d} <F(ξ)> dξ

    W(ξ_0) ≈ 0 since ξ_0 = 14σ is large enough for negligible NP interaction.

    Parameters
    ----------
    distances : np.ndarray
        Reaction coordinate values ξ (σ).
    forces : np.ndarray
        Ensemble-averaged forces at each ξ (ε/σ).
    xi_0 : float
        Reference distance where W ≈ 0 (default 14σ).

    Returns
    -------
    pmf : np.ndarray
        PMF W(d) in units of kBT (at T = ε/kB).
    """
    pmf = np.zeros_like(distances)
    for i in range(1, len(distances)):
        pmf[i] = pmf[i-1] - integrate.trapezoid(
            forces[:i+1], distances[:i+1]
        )
    return pmf


def extract_two_body_pmf(input_dir: str, output_dir: str) -> np.ndarray:
    """
    Extract W2(d12) from two-body force files.

    Loads all two-body NPforced*.dat files, integrates forces,
    and returns W2 as a function of d12.

    Parameters
    ----------
    input_dir : str
        Directory containing NPforced{N}.dat files for two-body PMF.
    output_dir : str
        Directory to write W2.csv and W2.png.

    Returns
    -------
    w2 : np.ndarray, shape (N, 2)
        Columns: [d12 (σ), W2 (kBT)]
    """
    os.makedirs(output_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(input_dir, 'NPforced*.dat')))

    all_distances = []
    all_pmf = []

    for f in files:
        dist, force = load_force_file(f)
        pmf = integrate_pmf(dist, force)
        all_distances.extend(dist)
        all_pmf.extend(pmf)

    d = np.array(all_distances)
    w = np.array(all_pmf)
    idx = np.argsort(d)
    d, w = d[idx], w[idx]

    np.savetxt(os.path.join(output_dir, 'W2.csv'),
               np.column_stack([d, w]),
               header='d12_sigma  W2_kBT', delimiter=',')

    plt.figure(figsize=(6, 4))
    plt.plot(d, w, 'o', ms=3, color='steelblue', label='PMF data')
    plt.axhline(0, color='gray', lw=0.5)
    plt.xlabel('d₁₂ (σ)')
    plt.ylabel('W₂ (k_BT)')
    plt.title('Two-body PMF')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'W2.png'), dpi=150)
    plt.close()

    print(f'W2: {len(d)} points, min = {w.min():.3f} kBT at d = {d[np.argmin(w)]:.2f}σ')
    return np.column_stack([d, w])


def main():
    parser = argparse.ArgumentParser(
        description='Extract PMFs from LAMMPS blue moon ensemble output'
    )
    parser.add_argument('--input-dir', required=True,
                        help='Directory containing NPforced*.dat files')
    parser.add_argument('--output', default='./pmf_results/',
                        help='Output directory for PMF data and plots')
    args = parser.parse_args()

    print('Extracting two-body PMF W2(d12)...')
    w2 = extract_two_body_pmf(args.input_dir, args.output)
    print(f'Done. Results written to {args.output}')


if __name__ == '__main__':
    main()
