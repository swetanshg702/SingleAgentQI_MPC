# SingleAgent QI-MPC

A Julia implementation of a **single-agent Model Predictive Control (MPC)** system enhanced with **Quasi-Interpolation (QI)**-based approximation for real-time control synthesis in discrete-time linear systems.

---

## Overview

This project implements an MPC controller that steers a 2D discrete-time linear system from an initial state `X0 = [0.5, -2.0]` to a target state `XT = [2.0, -2.0]`. Instead of solving a full QP optimization at every time step, the controller uses **quasi-interpolation over a precomputed grid** of optimal cost and control values — falling back to a direct QP solver only when the QI approximation is unavailable or unreliable.

### Key Ideas

- **Block MPC**: The horizon is divided into blocks of `ν` micro-steps (where `ν` is the system's reachability index), reducing the number of QP variables.
- **Quasi-Interpolation (QI)**: A Gaussian kernel-based interpolation scheme approximates optimal controls from a precomputed grid, enabling fast online evaluation.
- **Fallback QP**: When QI fails (e.g., query point is far from the grid), the system falls back to solving the block QP exactly using the Clarabel solver.

---

## Project Structure

```
SingleAgent_QI_MPC/
├── run_qi_mpc.jl           # Main entry point — runs the full simulation
└── src/
    ├── Constants.jl         # System parameters (A, B, Q, R, horizons, filenames)
    ├── SystemDynamics.jl    # Dynamics utilities (block B matrix, reachability index, Uref)
    ├── BlockQP.jl           # Block QP formulation and solver (via JuMP + Clarabel)
    ├── QI.jl                # Quasi-interpolation kernel and grid lookup functions
    ├── QIGridPrecompute.jl  # Precomputes cost and control grids offline
    └── MPCPlotting.jl       # Saves simulation results as interactive HTML plots
```

---

## System Parameters

Defined in `src/Constants.jl`:

| Parameter         | Value                        | Description                        |
|-------------------|------------------------------|------------------------------------|
| `X0`              | `[0.5, -2.0]`                | Initial state                      |
| `XT`              | `[2.0, -2.0]`                | Target state                       |
| `N`               | `20`                         | MPC horizon (macro blocks)         |
| `max_micro_steps` | `200`                        | Maximum simulation steps           |
| `A`               | `[0.1 -0.5; 0.2 0.8]`       | System matrix                      |
| `B`               | `[10.0; 0.0]`                | Input matrix                       |
| `Q`               | `Diagonal([5.0, 5.0])`       | State cost weight                  |
| `R_scalar`        | `10.0`                       | Control cost weight                |

---

## Module Descriptions

### `Constants.jl`
Defines all system-wide constants used across modules. Acts as the single source of truth for tunable parameters.

### `SystemDynamics.jl`
Core dynamics utilities:
- `build_block_B`: Constructs the lifted input matrix for a block of `ν` steps.
- `compute_reachability_index`: Finds the minimum `ν` such that the system is reachable to `XT` in `ν` steps.
- `compute_Uref`: Computes the reference control sequence that steers the system to `XT` in one block.
- `apply_micro_step`: Applies one discrete-time step `x_{t+1} = Ax_t + Bu_t`.

### `BlockQP.jl`
Formulates and solves the block MPC optimization problem using **JuMP** and the **Clarabel** solver. Includes stage and terminal costs with a 10× terminal weight on `Q`.

### `QI.jl`
Implements 2D Gaussian quasi-interpolation:
- Grid defined over `[-3, 3] × [-3, 3]` with spacing `h = 0.2`.
- `quasi_interpolate_cost`: Interpolates the optimal cost at a query point.
- `quasi_interpolate_control_vector`: Interpolates the optimal control vector at a query point.
- Falls back to the nearest finite neighbor if kernel weights are too small.

### `QIGridPrecompute.jl`
Precomputes the cost and control grids by solving the block QP at every grid node offline. These grids are then used by `QI.jl` at runtime.

### `MPCPlotting.jl`
Saves simulation results (state trajectory, control inputs, tracking error) as an interactive HTML file using Plotly or similar.

---

## Getting Started

### Prerequisites

- Julia **1.8+**
- The following Julia packages:
  ```
  JuMP
  Clarabel
  LinearAlgebra
  Printf
  ```

Install dependencies via the Julia REPL:
```julia
using Pkg
Pkg.add(["JuMP", "Clarabel"])
```

### Running the Simulation

From the project root directory:
```bash
julia run_qi_mpc.jl
```

This will:
1. Compute the reachability index `ν` for the system.
2. Precompute the QI cost and control grids.
3. Run the MPC simulation loop for up to `max_micro_steps = 200` steps.
4. Print a summary of the final state, tracking error, and accumulated cost.
5. Save an interactive HTML plot to `subsampled_qi_mpc.html`.

### Output

```
Starting Simulation...
--- Simulation Summary ---
Final State: [...]
Final Error: ...
Accumulated Cost: ...
Plot saved to: .../subsampled_qi_mpc.html
```

Open `subsampled_qi_mpc.html` in any browser to view the trajectory and control plots.

---

## Algorithm Flow

```
1. Compute reachability index ν
2. Precompute QI grids (cost + control) over [-3,3]²
3. For each MPC step:
   a. Query QI grid at current state z
   b. If QI succeeds → use interpolated control block
   c. If QI fails   → solve Block QP exactly (fallback)
   d. Apply ν micro-steps using chosen control
   e. Record state, control, and cost
4. Save results and plot
```

---

## Notes

- The `B` matrix was calibrated to match multi-agent settings. Using `B = [10.0; 0.0]` (100× larger) caused QI-approximated controls to overshoot and produce a persistent period-2 orbit — this is documented in `Constants.jl`.
- The QI grid range `[-3, 3]` must cover the expected state trajectory. States outside this range will fall back to the QP solver.

---

## License

This project was developed as part of a B.Tech final year project (BTP). All rights reserved.
