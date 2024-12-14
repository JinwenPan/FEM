# Finite Element Method

This project solves the 2D Helmholtz equation on a unit circle using the finite element method (FEM). The solution is specifically applied to describe a beam waveguide under homogeneous Neumann boundary conditions.

## Problem Description

The Helmholtz equation is given by:

`-∇²u - k²u = λu` (on the unit circle)

with homogeneous Neumann boundary conditions:

`∂u/∂n = 0` (on the boundary)

where:
- `u`: represents the eigenmode of the waveguide.
- `k²`: the wave number squared.
- `λ`: eigenvalue. 

## Features

- **Mesh Input**: The unstructured mesh is described in `unit_circle.txt`.
- **Matrix Assembly**: The stiffness and mass matrices are assembled using the external library Colsamm.
- **Output**: The numerical solution is saved in `eigenmode.txt`.

## How to Run

### Prerequisites
1. Install a C++ compiler supporting the C++17 standard or later.
2. Ensure all project files are in the same directory.

### Build
Run the following command to compile the project:
```bash
make
```

### Execute
To solve the problem, run the executable:
```bash
./waveguide delta epsilon
```
where `delta` is the wave number squared at the origin, and `epsilon` is the error tolerance of the solver for the eigenvalue problem. 