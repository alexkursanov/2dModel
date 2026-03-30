# Cardiac Model

2D cardiac tissue simulation with electromechanical coupling based on TNNP + Ekb models.

## Installation

```bash
pip install -e .
```

## Usage

```python
from src.cardiac_model import Tissue2D

tissue = Tissue2D(nx=50, ny=50)
tissue.apply_stimulus(region='center', radius=5)
results = tissue.simulate(t_max=500, dt=0.05)
```

## Project Structure

```
.
├── src/cardiac_model/     # Source code
│   ├── tissue2D.py         # 2D tissue simulation
│   ├── Ekb_mech.py         # Electromechanical cell model
│   ├── constants.py        # Model parameters
│   └── initial_conditions.py
├── tests/                  # Unit tests
├── examples/               # Example scripts
└── pyproject.toml
```

## Requirements

- numpy
- numba
- scipy
- matplotlib
