# **Guide to DFT Computation with Gaussian 16 on Bluehive**

---

**Platform:** Bluehive (University of Rochester HPC Cluster)  
**Software:** Gaussian 16  
**Purpose:** This guide is a comprehensive, step-by-step reference for performing Density Functional Theory (DFT) calculations using Gaussian 16 on the Bluehive high-performance computing (HPC) cluster.

---

## **Table of Contents**

1. [Quantum Mechanical Background](#quantum-mechanical-background)
2. [The Hartree-Fock (HF) Method](#the-hartree-fock-hf-method)
3. [Density Functional Theory (DFT)](#density-functional-theory-dft)
   - Hohenberg-Kohn Theorems
   - Kohn-Sham Approach
   - Total Energy in DFT
4. [Jacob’s Ladder of DFT Approximations](#jacob’s-ladder-of-density-functional-approximations)
5. [Basis Sets](#basis-sets)
6. [Accessing Bluehive](#accessing-bluehive)
   - Logging In
   - Partitions: Standard vs Preempt
   - Storage Limits
7. [Running Jobs Using Slurm](#running-jobs-using-slurm)
   - Job Script Template
8. [Building Molecules with Avogadro](#building-molecules-with-avogadro)
9. [Gaussian 16 Capabilities](#gaussian-16-capabilities)
10. [Gaussian Input (.gjf) File Structure](#gaussian-input-gjf-file-structure)
    - Link 0 Commands
    - Route Section
    - Title and Molecule Specification
11. [Charge and Multiplicity](#charge-and-multiplicity)
12. [Gaussian Output (.log) File](#gaussian-output-log-file)
13. [Checkpoint (.chk) Files](#checkpoint-chk-files)
14. [Geometry Optimization](#geometry-optimization)
15. [Molecular Orbitals (MO)](#molecular-orbitals-mo)
16. [Time-Dependent DFT (TDDFT)](#time-dependent-dft-tddft)
17. [Frequency and Vibrational Analysis](#frequency-and-vibrational-analysis)
    - IR and Raman Intensities
    - Preresonance Raman
18. [Resonance Raman Spectroscopy](#resonance-raman-spectroscopy)
    - Vertical Gradient Method
    - Adiabatic Shift Method
19. [Post-Processing Tools](#post-processing-tools)
    - Extracting Raman Spectra
20. [References and Resources](#references-and-resources)

---

## **1. Quantum Mechanical Background**

### **Mean Field Approximation**

- Simplifies the N-electron problem into **N separate one-electron problems**.
- Each electron moves in an average potential field created by all other electrons.

### **Electron Exchange Energy**

- A quantum mechanical contribution to total energy due to the **indistinguishability** of electrons and the **Pauli exclusion principle**.
- Requirement: the wavefunction must be **antisymmetric** under exchange of any two electrons.
- Physical Intuition: Identical fermions cannot occupy the same quantum state → they "keep their distance", lowering system energy.

> 💡 **Note:** In DFT, exact exchange is not used; instead, approximate **exchange-correlation functionals** are employed.

---

## **2. The Hartree-Fock (HF) Method**

### **Core Idea**

- Models the many-electron wavefunction as a **single Slater determinant**.
- Electrons occupy individual molecular orbitals.

### **Mean Field Treatment**

- Each electron experiences a **static average electric field** (mean field) from the rest of the electrons.
- No instantaneous Coulomb correlation included.

### **Strengths & Limitations**

- ✅ Exact treatment of **electron exchange energy**.
- ❌ Does **not include electron correlation**.
- ❌ Computationally expensive for large systems.

### **Transformation**

- Converts the N-electron Schrödinger-like problem into **N single-electron equations** (Fock equations).

> ⚠️ HF gives a good description of exchange but misses correlation – key for accurate energetics in weakly interacting systems.

---

## **3. Density Functional Theory (DFT)**

### **Fundamental Idea**

- Uses **electron density** (ρ(**r**)) instead of the full wavefunction.
- Electron density is a function of only **three spatial coordinates (x, y, z)**, independent of the number of electrons.

> 🔁 Reduction: N single-electron problems → **1 3D electron density function**

---

### **Hohenberg-Kohn Theorems**

**First Theorem:**

- The **ground-state electron density** uniquely determines **all properties** of a system, including its total energy.

**Second Theorem:**

- The true ground-state density is the one that **minimizes the energy functional**:

  $$
  E[\rho] = T[\rho] + V_{\text{ne}}[\rho] + J[\rho] + E_{\text{xc}}[\rho]
  $$

Where:

- $T[\rho]$: Kinetic energy functional
- $V_{\text{ne}}[\rho]$: Nucleus-electron attraction
- $J[\rho]$: Classical Coulomb (electron-electron repulsion)
- $E_{\text{xc}}[\rho]$: Exchange-correlation energy (unknown, must be approximated)

---

### **Kohn-Sham (KS) Approach**

- Introduces a **fictitious non-interacting system** that reproduces the same electron density as the real interacting system.
- Solves effective single-particle Schrödinger-like equations.

**Equation:**

$$
\left[ -\frac{1}{2} \nabla^2 + v_{\text{eff}}(\mathbf{r}) \right] \psi_i(\mathbf{r}) = \epsilon_i \psi_i(\mathbf{r})
$$

Where:

- $v_{\text{eff}}(\mathbf{r}) = v_{\text{ext}}(\mathbf{r}) + \int \frac{\rho(\mathbf{r'})}{|\mathbf{r}-\mathbf{r'}|} d\mathbf{r'} + v_{\text{xc}}(\mathbf{r})$

**Key Challenge:**

- Exchange-correlation ($E_{\text{xc}}$) term is **unknown**, and **approximations are essential** for accurate results.

---

### **Total Energy in DFT (Kohn-Sham Decomposition)**

The total energy of the system is computed as:

| Term                                   | Description                                           |
| -------------------------------------- | ----------------------------------------------------- |
| **KS Kinetic Energy**                  | Kinetic energy of non-interacting electrons           |
| **External Nuclei-Electron Potential** | Attraction between electrons and nuclei               |
| **Electron Coulomb Energy**            | Classical electrostatic repulsion (Hartree term)      |
| **Exchange–Correlation Energy**        | Non-classical quantum effects: exchange + correlation |

> 🛠 **Goal of DFT methods**: Accurately model $E_{\text{xc}}$ through functional approximations.

---

## **4. Jacob’s Ladder of Density Functional Approximations**

Hierarchy of increasingly accurate (but more expensive) approximations for $E_{\text{xc}}$:

| Level                                         | Description                                                                                                                 | Examples          |
| --------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- | ----------------- |
| **LDA**<br>Local Density Approximation        | Assumes uniform electron gas at each point; depends only on density $\rho(\mathbf{r})$                                      | S,VWN             |
| **GGA**<br>Generalized Gradient Approximation | Includes gradient of density $\nabla\rho(\mathbf{r})$, better for non-uniform systems                                       | PBE, BLYP         |
| **meta-GGA**                                  | Adds dependence on kinetic energy density or second derivative of $\rho$                                                    | TPSS, M06L        |
| **Hybrid**                                    | Mixes DFT with a portion (e.g., 20%) of **exact HF exchange**                                                               | B3LYP, PBE0       |
| **Range-Separated Hybrid (RSH)**              | Separates exchange interaction by distance via parameter **ω**:<br>– Short-range: DFT exchange<br>– Long-range: HF exchange | CAM-B3LYP, ωB97XD |

> 🔬 **Range**: Refers to spatial separation between electrons. Non-Coulomb parts of exchange functionals decay too rapidly at **large distances** → RSH fixes this.

---

## **5. Basis Sets**

Used to represent the **Kohn-Sham orbitals** as linear combinations of atomic orbital-like functions.

### **Types of Basis Sets**

#### **Zeta Levels**

- Number of functions per valence atomic orbital.

| Type                 | Description                                  | Example Basis Sets         |
| -------------------- | -------------------------------------------- | -------------------------- |
| **Double-Zeta (DZ)** | Two basis functions for each valence orbital | 6-31G, cc-pVDZ, def2-SVP   |
| **Triple-Zeta (TZ)** | Three basis functions for better accuracy    | 6-311G, cc-pVTZ, def2-TZVP |

#### **Polarized Basis Sets**

- Add **higher angular momentum functions** (d, f, p) to allow orbital distortion during bonding.

| Family              | Convention                              | Example             |
| ------------------- | --------------------------------------- | ------------------- |
| Pople (e.g., 6-31G) | Add `*` for polarization on heavy atoms | 6-31G\*, 6-311G\*\* |
| Dunning (cc-pVxZ)   | Automatically includes polarization     | cc-pVDZ, cc-pVTZ    |
| Ahlrichs (def2)     | Already polarized in standard versions  | def2-TZVP           |

> 💡 Tip: Use **polarized basis sets** for accurate geometries and energy differences.

#### **Diffuse Functions**

- Include very low-exponent Gaussians to describe **electrons far from nucleus**.
- Important for anions, excited states, non-covalent interactions.

| Family  | Convention                    | Example                                                             |
| ------- | ----------------------------- | ------------------------------------------------------------------- |
| Pople   | `+` adds diffuse function     | 6-31+G, 6-311++G                                                    |
| Dunning | `aug-` adds diffuse functions | aug-cc-pVDZ                                                         |
| def2    | Manually use **def2-TZVPD**   | Import from [Basis Set Exchange](https://www.basissetexchange.org/) |

> 🛑 Never perform anion calculations without **diffuse functions**.

---

## **6. Accessing Bluehive**

### **Prerequisites**

- You need a **UR account (URAD credentials)**.
- Request Bluehive access via:  
  📎 [https://info.circ.rochester.edu/](https://info.circ.rochester.edu/#BlueHive)

---

### **Logging In via SSH**

From Unix/Linux/Mac terminal:

```bash
ssh YourURADUsername@bluehive.circ.rochester.edu
```

You will be prompted for your password.

> 📍 Windows users can use **PuTTY** or **Windows Subsystem for Linux (WSL)**.

---

### **Compute Partitions: Standard vs Preempt**

| Partition    | Max Walltime        | Resources                     | Job Limits             | Notes                                                             |
| ------------ | ------------------- | ----------------------------- | ---------------------- | ----------------------------------------------------------------- |
| **standard** | 5 days (5-00:00:00) | CPU=120, Mem=3100G            | 2000 running/submitted | Longer time limit, but longer wait time                           |
| **preempt**  | 2 days (2-00:00:00) | CPU=120 (no mem limit stated) | 2000 (same)            | Jobs may be **preempted and restarted** (dangerous for Gaussian!) |

- You can use both partitions together.
- **Recommended Allocations:**
  - **Standard**: `36 CPUs, 100GB RAM` (short wait), or `56 CPUs, 190GB RAM` (expect hours wait)
  - **Preempt**: `48 CPUs, 160GB RAM` or `40 CPUs, 140GB RAM`

> 🔄 **Preempt jobs** can be **killed and rerun** – Gaussian may restart from scratch, risking wasted time.

---

### **Storage Limits (as of Oct 2025)**

| Directory           | Limit  | Grace Period           | Write Block                                    |
| ------------------- | ------ | ---------------------- | ---------------------------------------------- |
| `/home/username`    | 20 GB  | N/A                    | Immediate at 20 GB                             |
| `/scratch/username` | 200 GB | 7 days after exceeding | **Blocked at 1000 GB** – immediate job failure |

> ⚠️ **Exceeding 1000 GB on /scratch kills all running jobs instantly.**

#### **How to Check Usage**

```bash
quota                            # Shows storage quotas
du -h -s ~/scratch               # Check actual usage
```

#### **Tips for Managing Storage**

- Remove old `.RWF`, `.chk`, `.log` files from **failed jobs**.
- Download important results to **lab server or local machine**.
- Clear space regularly to avoid disruptions.

📌 **Transfer Files Guide:**  
📎 [https://info.circ.rochester.edu/#BlueHive/Transferring_Files/Transferring_Files/](https://info.circ.rochester.edu/#BlueHive/Transferring_Files/Transferring_Files/)

---

## **7. Running Jobs Using Slurm**

### **Job Submission Script Template**

Save as `run.sh`. Modify highlighted lines.

#### ✅ **Standard Partition Script**

```bash
#!/bin/bash
#SBATCH -p standard
#SBATCH -J your_job_name
#SBATCH -o Error-%j.out
#SBATCH --error=error-%j.err
#SBATCH --mail-type=all
#SBATCH --mail-user=your.email@address
#SBATCH --cpus-per-task=36
#SBATCH --mem-per-cpu=3000Mb
#SBATCH -t 5-00:00:00
#SBATCH --no-requeue

scontrol show job $SLURM_JOB_ID
module load gaussian/16-c01/avx2
g16 your_gaussian_input.gjf
```

#### 🔁 **Preempt Partition Script**

```bash
#!/bin/bash
#SBATCH -p preempt
#SBATCH -J your_job_name
#SBATCH -o Error-%j.out
#SBATCH --error=error-%j.err
#SBATCH --mail-type=all
#SBATCH --mail-user=your.email@address
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3500MB
#SBATCH -t 2-00:00:00

scontrol show job $SLURM_JOB_ID
module load gaussian/16-c01/avx2
g16 your_gaussian_input.gjf
```

### **Explanation of Key Directives**

| Directive                   | Purpose                                        |
| --------------------------- | ---------------------------------------------- |
| `#SBATCH -J job_name`       | Name of your job                               |
| `#SBATCH -o output.out`     | Output .out file                               |
| `#SBATCH --error=error.err` | Output .err file                               |
| `--mail-type=all`           | Email notifications (start, end, fail)         |
| `--cpus-per-task=N`         | Number of CPU cores to use                     |
| `--mem-per-cpu=XXXXMb`      | Memory per core                                |
| `-t D-HH:MM:SS`             | Walltime limit                                 |
| `--no-requeue`              | Prevent resubmission if preempted (important!) |

> ✅ Always include `--no-requeue` with Gaussian jobs — **restarting mid-calculation may be unsafe**.

### **Submit the Job**

Place `run.sh` and `your_gaussian_input.gjf` in the same directory, then:

```bash
sbatch run.sh
```

> ✅ Monitor job status with `squeue -u $USER`

---

## **8. Building Molecules with Avogadro**

### **Step-by-Step**

1. Open **Avogadro**.
2. Draw your molecule.
3. Use **Extensions → Auto-Optimize** to refine geometry with force fields (e.g., UFF).
4. Go to **Extensions → Gaussian**.
5. Copy the coordinates and input setup.
6. Save as **`input.gjf`**.

> 💡 Avogadro generates template input files – perfect starting point.

---

## **9. Gaussian 16 Capabilities**

Gaussian 16 supports a wide range of computational methods:

| Method Type                  | Supported Methods                                                       |
| ---------------------------- | ----------------------------------------------------------------------- |
| **Molecular Mechanics**      | Amber, UFF, Dreiding                                                    |
| **Semi-Empirical**           | AM1, PM6, PM7, DFTB                                                     |
| **Hartree-Fock (HF)**        | Standard, including gradients (`G`) and frequencies (`F`)               |
| **Density Functional (DFT)** | All common functionals; long-range corrections; dispersion (e.g., D3BJ) |
| **Post-HF Methods**          | MP2, MP4, CCSD, CCSD(T), CASSCF, EOM-CCSD (excited states)              |
| **TD-DFT**                   | Excited state energies, gradients                                       |
| **High-Accuracy Models**     | G1-G4, CBS, W1 series                                                   |

> `E`: Energy  
> `G`: Analytic Gradients  
> `F`: Analytic Frequencies  
> `F†`: Reimplemented with analytic frequencies

---

## **10. Gaussian Input (.gjf) File Structure**

Example:

```gjf
%Chk=h2o.chk
%nprocshared=48
%mem=160GB
# B3LYP/6-31G(d) Opt

Water Optimization

0 1
O  -0.464   0.177   0.0
H  -0.464   1.137   0.0
H   0.441  -0.143   0.0

```

### **Sections Explained**

#### **Link 0 Commands (Start with %)**

| Command            | Meaning                                              |
| ------------------ | ---------------------------------------------------- |
| `%Chk=file.chk`    | Name of checkpoint file                              |
| `%OldChk=file.chk` | Read geometry/wavefunction from previous calculation |
| `%NProcShared=N`   | Number of CPU cores to use                           |
| `%Mem=sizeGB`      | Total memory allocation (e.g., 160GB)                |

> ❗ These must be first in the file.

---

#### **Route Section (`#`)**

Specifies the method, basis set, and job type.

Example:

```
#p B3LYP/6-31G(d) Opt Pop=Regular
```

- **`B3LYP`**: Functional
- **`/6-31G(d)`**: Basis set (equal to 6-31G\*)
- **`Opt`**: Geometry optimization
- **`Pop=Regular`**: Compute and print molecular orbitals
- **`p`**: Prints more output (useful for debugging)

> 📎 Full input syntax: [https://gaussian.com/input/](https://gaussian.com/input/)

---

#### **Title Section**

- One line describing the job.
- Can be anything (e.g., "Water Optimization").

---

#### **Molecule Specification**

- Charge and multiplicity on first line.
- Atom list: Element, x, y, z (in angstroms).
- Must end with a **blank line**.

> ❗ No trailing spaces; blank line required for parsing.

---

## **11. Charge and Multiplicity**

| Multiplicity     | Formula                              | Examples                            |
| ---------------- | ------------------------------------ | ----------------------------------- |
| **Singlet (S₀)** | 2S + 1 = 1 → S = 0 → 0 unpaired e⁻   | Closed shell neutral molecule       |
| **Triplet (T₁)** | 2S + 1 = 3 → S = 1 → 2 unpaired e⁻   | O₂, excited states                  |
| **Doublet**      | 2S + 1 = 2 → S = 0.5 → 1 unpaired e⁻ | Radicals, cations with odd e⁻ count |

### **Common Cases**

| Species              | Charge | Multiplicity |
| -------------------- | ------ | ------------ |
| Neutral closed-shell | 0      | 1            |
| Neutral radical      | 0      | 2            |
| Triplet state        | 0      | 3            |
| Anion (M⁻)           | -1     | 2            |
| Dianion (M²⁻)        | -2     | 1 or 3       |
| Cation (M⁺)          | +1     | 2            |
| Dication (M²⁺)       | +2     | 1 or 3       |

> ⚠️ For **anions**, **always use a basis set with diffuse functions** (e.g., 6-31+G\*, aug-cc-pVDZ)

---

## **12. Gaussian Output (.log) File**

- Generated upon job execution.
- Contains **all results**: energies, coordinates, frequencies, orbitals, transitions.

### **Sign of Successful Completion**

```log
Normal termination of Gaussian 16 at Wed Nov  1 16:26:50 2023.
```

> ❌ If this does **not appear**, the job failed or is still running.

### **Visualization**

- Load `.log` or checkpoint file into **Avogadro**, **GaussView**, or **ChemCraft** to visualize:
  - Optimized geometry
  - Vibrational modes
  - Excited states
  - Molecular orbitals

---

## **13. Checkpoint (.chk) Files**

### **Purpose**

- Binary file containing **wavefunction, orbitals, density, geometry**, etc.
- Enables restarts and advanced analysis.

### **Uses**

- Initial guess for subsequent calculations (via `geom=check guess=read`)
- TDDFT, frequency, population analysis
- Visualization of **molecular orbitals** (requires `.chk` or `.fchk`)

### **Platform Compatibility**

- `.chk` files are **binary and Linux-compatible only**.
- For Windows/Mac: convert to **formatted checkpoint (.fchk)**.

### **Convert .chk to .fchk**

#### On Bluehive (Terminal)

```bash
export GAUSS_MEMDEF=2000MW
module load gaussian/16-c01/avx2
formchk -3 input.chk output.fchk
```

> ✅ `-3`: Creates double-precision `.fchk`

#### Using Slurm (Recommended for Large Files)

Save as `job.sh`:

```bash
#!/bin/bash
#SBATCH -p preempt
#SBATCH -J formchk
#SBATCH -o Error-%j.out
#SBATCH --error=error-%j.err
#SBATCH --mail-type=all
#SBATCH --mail-user=lcai5@ur.rochester.edu
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=3000Mb
#SBATCH -t 0-48:00:00

scontrol show job $SLURM_JOB_ID
export GAUSS_MEMDEF=2000MW
module load gaussian/16-c01/avx2
formchk -3 check.chk check.fchk
```

Submit: `sbatch job.sh`

---

## **14. Geometry Optimization**

Optimize molecular structure to minimize total energy.

### **Route Section Keyword**

```
Opt
```

### **Example Input**

```gjf
%Chk=h2o.chk
%nprocshared=48
%mem=160GB
# opt B3LYP/6-31G(d)

Water Optimization

0 1
O  -0.464   0.177   0.0
H  -0.464   1.137   0.0
H   0.441  -0.143   0.0

```

### **Starting from a Previous .chk File**

Useful to continue optimization or change method.

```gjf
%OldChk=old.chk
%NProcShared=48
%Mem=160GB
%Chk=new.chk
#p opt b3lyp/def2SVP geom=check guess=read

Restarted Optimization

0 1
```

> ✅ `geom=check`: read geometry from checkpoint  
> ✅ `guess=read`: read wavefunction for faster convergence

### **Solvent Effects**

PCM performs a reaction field calculation using the integral equation formalism model (IEFPCM). To include solvent effects using PCM model:

```
#p opt b3lyp/def2SVP SCRF=(Solvent=Ethanol) geom=check guess=read
```

Available solvents can be found in [https://gaussian.com/scrf/](https://gaussian.com/scrf/)

### **Convergence Criteria**

Successful optimization includes:

```log
Item                     Value     Threshold   Converged?
Maximum Force           0.000002   0.000015     YES
RMS Force               0.000000   0.000010     YES
Maximum Displacement    0.000008   0.000060     YES
RMS Displacement        0.000002   0.000040     YES
```

All lines should say **YES**.

---

## **15. Molecular Orbitals (MO)**

To analyze and visualize orbitals:

Add to route section:

```
pop=regular
```

or for more detail:

```
pop=(regular,mocoeffs)
```

### **Steps**

1. Run calculation with `pop=regular`
2. Generate `.fchk` from `.chk` (see below)
3. Open in **Avogadro** → Extensions → Surfaces → Show MOs

> 🔎 MOs help understand bonding, reactivity, and excitation origins.

---

## **16. Time-Dependent DFT (TDDFT)**

Calculates **excited state energies** and properties.

### **Route Section**

```
td=(nstates=N,root=M)
```

- `N`: Number of excited states to compute
- `M`: State for energy/gradient (e.g., S₁ = root=1)

### **Example**

```gjf
%OldChk=benzyl_d0.chk
%NProcShared=10
%Mem=10GB
%Chk=benzyl_d3_vg.chk
#p td=(nstates=6,root=3) b3lyp/cc-pvtz geom=check guess=read

Benzyl radical doublet: gradient calculation on third excited state at ground state geometry.

0 2
```

> 💡 `root=3` means TDDFT will perform a single-point gradient on the 3rd excited state.

---

### **Sample Output**

```log
Excited State   1:      Singlet-A      1.9615 eV  632.07 nm  f=0.8381  <S**2>=0.000
    478 -> 481       -0.12828
    479 -> 482        0.30228
    479 -> 483       -0.11884
    480 -> 481        0.57529
    480 -> 483        0.14482
This state for optimization and/or second-order correction.
Total Energy, E(TD-HF/TD-DFT) =  -8409.97137589
```

- **<S²>**: Spin contamination (should be near 0 for singlets)
- **f**: Oscillator strength (bright transitions > 0.1)

> 🔧 TDDFT can be combined with `Opt` to **optimize excited-state geometries** – computationally expensive (5–10× ground-state opt).

---

## **17. Frequency and Vibrational Analysis**

Computes:

- Harmonic vibrational frequencies
- IR/Raman intensities
- Zero-point energy (ZPE)
- Thermodynamic corrections

---

### **A. Frequencies and IR Intensity**

```gjf
%OldChk=opt.chk
%NProcShared=48
%Mem=160GB
%Chk=freq.chk
#p freq b3lyp/def2SVP geom=check guess=read

Frequencies and IR intensity

0 1
```

> ✅ Must use **optimized geometry** (from `Opt` step) as input.

---

### **B. Preresonance Raman**

Resonant enhancement without actual resonance (excitation near but not at absorption peak).

```gjf
#p freq=raman cphf=rdfreq b3lyp/def2SVP geom=check guess=read

Preresonance Raman intensity

0 1
```

> 📐 `CPHF=RdFreq`: Requests frequency-dependent polarizability derivatives.

---

### **C. Frequency Convergence Issues**

If frequency job fails to converge, add:

```
Int=Acc2E=13 CPHF(MaxInv=10000)
```

This increases integral accuracy and maximum iterations.

---

### **Scaling Frequencies**

Computed harmonic frequencies are typically **overestimated**. Apply scaling factors:

📊 Resource:  
📎 [https://cccbdb.nist.gov/vsfx.asp](https://cccbdb.nist.gov/vsfx.asp)

| Functional/Basis | Scaling Factor |
| ---------------- | -------------- |
| B3LYP/6-31G(d)   | 0.9603         |
| wB97XD/TZVP      | 0.9550         |

Multiply all frequencies by the appropriate factor.

---

## **18. Resonance Raman Spectroscopy**

Resonance Raman (RR) enhances Raman intensity when laser frequency matches an electronic transition.

Two proper methods:

---

### **Method 1: Vertical Gradient (VG)**

1. **Ground State Frequency** (Hessian at GS geometry)

```gjf
%OldChk=gsopt.chk
%NProcShared=48
%Mem=160GB
%Chk=gsfreq.chk
#p freq b3lyp/def2SVP geom=check guess=read

Simple

0 1
```

2. **Vertical Gradient Calculation** (Force at GS geometry on excited state)

```gjf
%OldChk=gsfreq.chk
%NProcShared=48
%Mem=160GB
%Chk=s1force.chk
#p force td=(nstates=6,root=3) b3lyp/def2SVP geom=check guess=read

S1 Optimization

0 1
```

3. **Resonance Raman with Vertical Gradient**

```gjf
%OldChk=gsfreq.chk
%NProcShared=48
%Mem=160GB
%Chk=rrvg.chk
#p freq=(FC,ReadFC,ReadFCHT) geom=check

Resonance Raman with vertical gradient method

0 1
! ReadFCHT input keywords go here
Method=VerticalGradient
Spectroscopy=ResonanceRaman
ResonanceRaman=(Omega=31250.0)
Spectrum=(Broadening=Lorentzian,HWHM=10.0)
! Excited state checkpoint file name read here
s1force.chk
```

> 🔧 `Omega=31250.0`: Laser frequency in cm⁻¹ (e.g., 31250 ≈ 320 nm)

---

### **Method 2: Adiabatic Shift (AS)**

1. **Ground State Frequency** (same as above)

2. **Excited State Geometry Optimization**

```gjf
%OldChk=gsfreq.chk
%NProcShared=48
%Mem=160GB
%Chk=s1opt.chk
#p opt td=(nstates=6,root=3) b3lyp/def2SVP geom=check guess=read

S1 Optimization

0 1
```

3. **Resonance Raman with Adiabatic Shift**

```gjf
%OldChk=gsfreq.chk
%NProcShared=48
%Mem=160GB
%Chk=rras.chk
#p freq=(FC,ReadFC,ReadFCHT) geom=check

Resonance Raman with adiabatic shift method

0 1
! ReadFCHT input keywords go here
Method=AdiabaticShift
Spectroscopy=ResonanceRaman
ResonanceRaman=(Omega=31250.0)
Spectrum=(Broadening=Lorentzian,HWHM=10.0)
! Excited state checkpoint file name read here
s1opt.chk
```

📘 Source: _Coordination Chemistry Reviews_, Volume 254, Issues 21–22, November 2010, Pages 2505–2518

---

## **19. Post-Processing Tools**

### **Extract Raman Spectra**

Use the script:
📎 [https://github.com/Alchemist-Aloha/raman_from_gaussian](https://github.com/Alchemist-Aloha/raman_from_gaussian)

- Parses `.log` file for Raman activities.
- Generates **plots** with Lorentzian/Gaussian broadening.
- Customizable FWHM (Full-width at half-maximum), linewidth.

> 🧩 Perfect for non-resonance and resonance Raman visualization.

---

## **20. References and Resources**

- **Bluehive Documentation:**  
  📎 [https://info.circ.rochester.edu/#BlueHive](https://info.circ.rochester.edu/#BlueHive)

- **Gaussian 16 Manual:**  
  📎 [https://gaussian.com/](https://gaussian.com/)

- **Basis Set Exchange:**  
  📎 [https://www.basissetexchange.org/](https://www.basissetexchange.org/)

- **NIST Frequency Scaling Factors:**  
  📎 [https://cccbdb.nist.gov/vsfx.asp](https://cccbdb.nist.gov/vsfx.asp)

- **Resonance Raman Reference Paper:**  
  _Coord. Chem. Rev._, **2010**, _254_(21–22), 2505–2518. DOI: [10.1016/j.ccr.2009.11.015](https://doi.org/10.1016/j.ccr.2009.11.015)

---

## **Final Tips**

✅ Always:

- Use `--no-requeue` for Gaussian jobs
- Monitor storage usage with `quota`
- Scale frequencies appropriately
- Use diffuse functions for anions/excited states
- Convert `.chk` to `.fchk` for visualization

🧠 Think carefully about:

- Functional choice (use RSH for charge transfer)
- Basis set completeness
- Excited state convergence
- Memory and CPU allocation

---

**End of Guide**  
Date: October 2025
