# mandacaru-paw

Projector augmented-wave (PAW) datasets for
[Mandacaru](https://github.com/seixas-research/mandacaru), one file per element
for every element with **Z ≤ 92** (H through U). The datasets are generated
from scratch by Mandacaru's own LDA radial atomic solver and its
`mandacaru.pseudopotentials.paw` module; nothing here is copied
from another PAW code.

They live in this repository, not in Mandacaru itself, because of their size:
about 190 MB for the 92 files, against the 100 MB limit of a PyPI release.
The Troullier–Martins (NCPP) library, 11 MB, still ships inside the package.

## Using the datasets

Point Mandacaru at a checkout of this repository once; it creates a symbolic
link `library/paw` inside the installed package, and the loaders take it from
there:

```bash
git clone git@github.com:seixas-research/mandacaru-paw.git
python -m mandacaru.pseudopotentials.link_library --paw mandacaru-paw
```

Use `--files` to link each dataset individually instead of the directory,
`--force` to replace an existing link, `--status` to see what each family
folder serves. Alternatively set `MANDACARU_PSEUDO_PATH` to a directory that
contains this checkout as its `paw/` subfolder.

Then, in a calculation:

```python
from ase.build import molecule
from mandacaru import Mandacaru

atoms = molecule("H2O")
atoms.center(vacuum=4.0)          # the cell is the real-space box
atoms.calc = Mandacaru(method="adapt-vqe",
                       basis="PAW",
                       h=0.25)
atoms.get_total_energy()          # eV, valence-only Hamiltonian
```

The family is selected **as a basis**: `basis="PAW"`, or
`basis={"name": "PAW", "size": "DZ"}` for a larger valence basis, exactly like
an all-electron family. Its guide is the *Pseudopotentials* page of the
Mandacaru manual (`docs/source/guide/pseudopotentials.md`).

## What is in a file

Each `<Symbol>.parquet` is a self-describing Mandacaru pseudopotential record
(format `mandacaru-pseudopotential`, version 2, `family = "paw"`), readable
with `mandacaru.pseudopotentials.paw.get_paw(symbol)` or the
generic `io.load_pseudopotential(path)`. The table holds the radial grid
(3000 points, 0.01 bohr spacing) and, per angular momentum `l`:

- two all-electron partial waves `ae_wave_l{l}_{0,1}` and their smooth
  counterparts `pseudo_wave_l{l}_{0,1}` at the reference energies ε₁ (the
  bound valence eigenvalue) and ε₂ = ε₁ + 1 Ha,
- the dual projectors `projector_l{l}_{0,1}` (⟨p̃ᵢ|φ̃ⱼ⟩ = δᵢⱼ) and the raw
  projectors they were built from,
- the 2×2 one-center matrices: overlap correction q_ij, kinetic-energy
  difference ΔT_ij, potential difference ΔV_ij and the coupling D_ij.

The record also carries the local potential (screened and unscreened), the
frozen core density and its smooth counterpart, the pseudo valence density,
the monopole compensation charge and radius, the frozen one-center energy,
the cutoff radii and the reference energies. All quantities are in atomic
units (bohr, hartree); Mandacaru converts to eV and Å at its user-facing layer.

## Construction, in brief

Following P. E. Blöchl, *Phys. Rev. B* **50**, 17953 (1994), in the
frozen-core, one-center-expansion formulation:

1. all-electron LDA atom; frozen core density and valence partial waves at
   two energies per `l`;
2. smooth partial waves as spherical-Bessel expansions inside `r_c`, matched
   at `r_c` but **not** norm-conserving (a scaled generalized-norm condition
   keeps the augmented overlap positive definite);
3. projectors dual to the smooth waves inside the augmentation sphere, by the
   Vanderbilt construction;
4. a smooth local potential and a monopole compensation charge;
5. one-center terms ΔT, q and the Hartree/exchange–correlation part of D,
   evaluated at the LDA reference atom and stored as fixed matrices.

In a molecule the nonlocal term is Σ|p̃ᵢ⟩D_ij⟨p̃ⱼ| and the overlap becomes
S̃ + Σ⟨φ|p̃ᵢ⟩q_ij⟨p̃ⱼ|φ⟩, the generalized eigenvalue problem of PAW.

Every dataset was checked on generation: projector duality to 1e-8, the bound
eigenvalue reproduced to better than 1e-6 Ha with no ghost state, the
all-electron wave reconstructed from the smooth one to 1e-4, and logarithmic
derivatives matched at both reference energies.

## Approximations

Relative to a full PAW implementation: the one-center Hartree and
exchange–correlation terms are frozen at the LDA reference (no
self-consistent D[ρ]); compensation charges are monopole only; the core is
frozen with no nonlinear core correction (the smooth core density is stored
but unused); LDA only; no relativistic terms; no projectors above the valence
`l`. Mandacaru then solves the valence problem at the Hartree–Fock / FCI level
on these LDA-generated datasets.

## Regenerating

From the Mandacaru repository, with the `mandacaru` environment:

```python
from mandacaru.pseudopotentials.paw import build_paw_library
from mandacaru.pseudopotentials.io import library_elements

build_paw_library(library_elements(92), directory="path/to/mandacaru-paw")
```

Generation takes 5 s for light elements and up to 150 s for the heaviest, about
90 minutes for the full set, with under 1 GB of memory. The 2026-09-14 set
was generated with zero failures.

## License

MIT, see `LICENSE`.
