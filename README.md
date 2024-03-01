# flow_dynamic

Contains serial implementations of following algorithms:
* 2D vorticity transfer equation solver
* 2D thermal conductivity equation solver
* 2D combined fluid and thermal solver

Note that present implementations are intended for academic purposes, as such they are not meant to be used in any sort of high-performance production code.

## Compilation

* Recommended compiler: Intel C++ Compiler
* Requires C++17 support

## Usage
Two dimensional formulation:<br>
$$ \vec{V} = \vec{V}(x, y, t) = (u, v, 0)^{T}, \quad \vec{W} = \nabla \times \vec{V} = (0, 0, \omega), \quad \omega = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}.$$

Vorticity transfer equation:<br>
$$\frac{\partial \omega}{\partial t} + (\vec{V} \cdot \nabla)\omega = \nu \Delta \omega, $$
$$ \omega = -\Delta \psi, \quad u = \frac{\partial psi}{\partial y}, \quad v = -\frac{\partial v}{\partial x},$$
where<br>
$\nu$ - kinematic viscosity, $\psi$ - fluid potential.<br>