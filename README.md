# Sitnikov Problem Simulator

## Overview
This program simulates and visualizes the Sitnikov problem - a special case of the three-body problem in celestial mechanics. It consists of two primary bodies moving in elliptical orbits around their center of mass, while a third body moves along a line perpendicular to the orbital plane of the primaries.

## Features
- Parallel computation using MPI (Message Passing Interface)
- Configurable orbital parameters via YAML configuration
- Phase portrait generation
- Multiple numerical integration methods support
- Progress tracking for parallel computations

## Requirements
- Python 3.x
- NumPy
- SciPy
- Matplotlib
- mpi4py
- PyYAML

All requirements are listed in `requirements.txt` file for Python environment.

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/sitnikov-problem.git
cd sitnikov-problem
```

2. Create and activate a virtual environment (optional but recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Unix/macOS
# or
.\venv\Scripts\activate  # On Windows
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Configuration
The simulation parameters can be configured in `config.yaml`:

- `orbit`: Primary bodies' orbital parameters
  - `a`: Semi-major axis
  - `e`: Eccentricity (0 ≤ e < 1)

- `grid`: Integration parameters
  - `nz`: Number of position grid points
  - `ndot`: Number of velocity grid points
  - `z0`, `z1`: Position range
  - `dotz0`, `dotz1`: Velocity range
  - `t0`, `tf`: Time range

- `integrator`:
  - `method`: Integration method. Supported methods from `scipy.integrate.solve_ivp`:
    - 'RK45' (default): Explicit Runge-Kutta method of order 5(4)
    - 'RK23': Explicit Runge-Kutta method of order 3(2)
    - 'DOP853': Explicit Runge-Kutta method of order 8
    - 'Radau': Implicit Runge-Kutta method of the Radau IIA family of order 5
    - 'BDF': Implicit multi-step variable-order method
    - 'LSODA': Adams/BDF method with automatic stiffness detection

    For detailed information about these methods, see [SciPy's solve_ivp documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)


- `plotter`: Visualization parameters
  - `zmin`, `zmax`: Position plot range
  - `dotzmin`, `dotzmax`: Velocity plot range
  - `dpi`: Output image resolution

- `parallel`:
  - `track`: Enable/disable progress tracking

## Usage

### Single-core execution:
```bash
python main.py
```

### Parallel execution:
```bash
mpirun -n <number_of_processes> python main.py
```

### Plotting results:
```bash
python plot.py
```

## Output
The program generates:
- Data files in the `result/` directory containing trajectories
- A phase portrait plot showing the system's dynamics

## Mathematical Background
The Sitnikov problem describes the motion of a massless particle moving perpendicular to the plane of two equal masses (primaries) in elliptical orbits. The equations of motion are:

$$\ddot{z} = -\frac{z}{(r(t, e)^{2} + z^{2})^{3/2}}$$

where:
- $z$ is the position of the third body
- $r$ is the distance of the primaries from the center of mass
- $t$ is time
- $e$ is orbital eccenrisity

