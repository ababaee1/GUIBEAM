# GUIBEAM: FEM GUI Tool for Beam Static Analysis

GUIBEAM is a MATLAB GUI tool for analyzing beam structures using the Finite Element Method (FEM). Perfect for engineers, researchers, and students, GUIBEAM offers a user-friendly interface with robust tools for mesh generation, customizable loading, boundary conditions, and detailed post-processing of displacements, stresses, strains, and forces.

<img src="https://github.com/user-attachments/assets/f75267b1-10b8-471a-8c12-fcce7589a4b7" width="60%">

## 🎯 Features

### Mesh Generation
- Mesh Types: Two structured mesh options with customizable density for refined analysis.
- Plane Stress & Plane Strain: Choose the formulation based on the problem type for accurate simulations.

### Boundary & Loading Conditions
- Displacement Constraints: Apply fixed, pinned, or custom displacements to nodes.
- Force Boundary Conditions: Set concentrated and distributed forces in x/y directions.
- Body Forces: Include gravity or other uniform forces across the beam.

### 📈 Distributed Load
- Equation-Based Loading: Define distributed loads along the beam using custom functions—ideal for varying load distributions.

### Node & Element Numbering
- Automatic Labeling: Easily identify nodes and elements for assigning conditions and checking mesh accuracy.

### Post-Processing
- Displacement Visualization: View x/y displacements and structural deformations under load.
- Stress & Strain Analysis: Inspect stress (σ_x, σ_y, τ_xy) and strain (ε_x, ε_y, γ_xy) distributions.
- Reaction Forces: Display calculated reaction forces at boundaries.

### Output Options
- Graphical Plots: Visualize deformed shapes, stress, and strain contours for a clear view of results.
- Numerical Displays: Access detailed lists of displacements, stresses, and strains.

### Interactive GUI
- Dynamic Inputs: Set mesh parameters, materials, loads, and conditions directly in the GUI.
- Reset & Clear: Reset or clear inputs for new analyses without hassle.
