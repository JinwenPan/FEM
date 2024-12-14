# Solving the 2D Helmholtz Equation on a Unit Circle

This project solves the 2D Helmholtz equation on a unit circle using the finite element method (FEM). The solution is specifically applied to describe a beam waveguide under **homogeneous Neumann boundary conditions**.

---

## Problem Description

The Helmholtz equation is given by:

\[
-\nabla^2 u - k^2 u = \lambda \quad \text{on the unit circle,}
\]

with homogeneous Neumann boundary conditions:

\[
\frac{\partial u}{\partial n} = 0 \quad \text{on the boundary.}
\]

Here:
- \( u \): represents the eigenmode of the waveguide.
- \( k^2 \): the wave number squared.

---

## Features

- **Mesh Input**: The unstructured mesh is described in `unit_circle.txt`.
- **Matrix Assembly**: The stiffness and mass matrices are assembled using the external library [Colsamm](https://github.com/marhol/colsamm).
- **Output**: The numerical solution is saved in `eigenmode.txt`.

---

## How to Run

### Prerequisites
1. Install a C++ compiler supporting the C++11 standard or later.
2. Clone and install the [Colsamm library](https://github.com/marhol/colsamm).
3. Ensure all project files are in the same directory.

### Build
Run the following command to compile the project:
```bash
make
