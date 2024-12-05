# Cell Motility Simulations

This repository contains the code and methodology for simulating the collective behavior of cells modeled as semiflexible ring polymers. The project investigates the impact of motility force, cell density, and confinement on cell dynamics and collective motion.

---

## Overview

This project models cells as two-dimensional ring polymers using molecular dynamics simulations. The approach incorporates:
- Shape anisotropy,
- Propulsion forces, and 
- Various interaction potentials.

These tools are used to study clustering and motion patterns under free and confined conditions.

---

## Features

### Cell Modeling  
- Cells are modeled as ring polymers, representing circular or elongated shapes based on bending potentials.

### Simulation Scenarios  
- Includes free migration and confinement inside a semiflexible boundary.

### Emergent Behavior  
- Observes transitions in motion patterns with varying motility forces, densities, and boundary stiffness.



---

## Model Description

- **Bonding Energy**: Ensures linear connections between adjacent monomers using a quadratic potential.  
- **Bending Energy**: Maintains cell shape through three-body interactions.  
- **Area Constraint**: Enforces consistent cell area using a harmonic potential.  
- **Propulsion Force**: Applies a motility force aligned with cell polarity.  

---

## Simulation Details

- **Tools**: LAMMPS Molecular Dynamics Simulator  
- **Parameters**:  
  - Circular cells:  
    - Preferred angle = 180°  
    - Area = 127r_b²  
  - Elongated cells:  
    - Preferred angle = 120°  
    - Area = 88.1r_b²  
  - System size: Lx = Ly = 400r_b  
  - Motility force: 0 to 48ε/r_b  
  - Temperature: k_BT = ε  

---

## Key Findings

### Free Migration  
- Elongated cells demonstrate stronger collective motion than circular cells due to shape anisotropy.  
- Increased motility force enhances cluster size and velocity alignment.

### Confinement  
- At low packing fractions, cells exhibit vortical motion.  
- High packing fractions result in ballistic, directed motion.

---
