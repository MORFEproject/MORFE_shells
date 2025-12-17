# MORFE Shells
**Direct Parametrisation of Invariant Manifolds for Nonlinear Shell Structures**

DPIM-Shell is a nonlinear model order reduction (MOR) library that combines the **Direct Parametrisation of Invariant Manifolds (DPIM)** with a **seven-parameter solid–shell finite element formulation**.  
The library provides implementations in both **Julia** and **MATLAB**, aiming at high-performance computation as well as clarity and ease of debugging.

---

## ✨ Key Features

- Direct parametrisation of invariant manifolds for large-scale nonlinear finite element systems  
- Seven-parameter solid–shell elements with **Enhanced Assumed Strain (EAS)** correction  
- Intrusive reduced-order modelling (ROM) for geometrically nonlinear shell structures  
- Treatment of **non-autonomous terms** directly within the DPIM framework  
- Efficient computation of nonlinear forced responses via **Harmonic Balance Method (HBM)**  

---

## 🧩 Julia Implementation

The Julia implementation focuses on **computational efficiency and scalability**.  
It follows recent developments by Vizzaccaro *et al.* on DPIM, where **non-autonomous terms are intrinsically included** in the parametrisation-based model reduction process.

In combination with the seven-parameter solid–shell finite element formulation, the EAS-corrected element operators are embedded directly into the ROM construction.  
The original benchmark example is adapted from Jain *et al.*, consisting of a **singly-curved shell model** that has been verified to exhibit a **1:2 internal resonance**.

To further increase the complexity of the problem, an additional configuration is considered in which the shell thickness varies **linearly along the direction of curvature**, allowing for a more challenging nonlinear dynamic response.

---

## 🧪 MATLAB Implementation

The MATLAB implementation mirrors the DPIM workflow used in the Julia code, with an emphasis on **readability, transparency, and ease of prototyping**.  
It is particularly suitable for understanding the DPIM algorithmic structure and for developing user-defined extensions.

Since DPIM involves highly nested and loop-intensive computations, the original `.m` files used for evaluating **high-order nonlinear internal forces** can be replaced by compiled **`mex` functions** for improved performance.  
For clarity, the corresponding `.m` source files are provided alongside the `mex` versions to document the finite element implementation in detail.

For the computation of nonlinear forced responses of the reduced-order models, an **in-house harmonic balance solver** is employed instead of standard collocation-based continuation methods.  
This choice is motivated by the poor convergence behaviour often encountered when using packages such as Matcont or COCO for complex resonance scenarios.

In addition, the open-source software NLvib developed by Krack *et al.* is integrated as an alternative solver, allowing users to select their preferred computational strategy.

---

## 📚 References

1. Xia, Z. X., Touzé, C., Cong, Y., *et al.* (2025).  
   *Reduced order modelling for shell finite element structures using the direct parametrisation of invariant manifolds:  
   Hardening/softening transition, resonant dynamics and mode selection.*  
   **Thin-Walled Structures**, 114329.

2. Vizzaccaro, A., Gobat, G., Frangi, A., *et al.* (2023).  
   *Direct parametrisation of invariant manifolds for generic non-autonomous systems including superharmonic resonances.*  
   **arXiv preprint**, arXiv:2306.09860.

3. Jain, S., & Haller, G. (2022).  
   *How to compute invariant manifolds and their reduced dynamics in high-dimensional finite element models.*  
   **Nonlinear Dynamics**, 107(2), 1417–1450.

4. Krack, M., & Gross, J. (2019).  
   *Harmonic Balance for Nonlinear Vibration Problems.*  
   Cham: **Springer International Publishing**.

5. Detroux, T., Renson, L., Masset, L., *et al.* (2015).  
   *The harmonic balance method for bifurcation analysis of large-scale nonlinear mechanical systems.*  
   **Computer Methods in Applied Mechanics and Engineering**, 296, 18–38.

6. Büchter, N., Ramm, E., & Roehl, D. (1994).  
   *Three-dimensional extension of nonlinear shell formulation based on the enhanced assumed strain concept.*  
   **International Journal for Numerical Methods in Engineering**, 37(15), 2551–2568.


