# Rendezvous control on Riemannian manifolds

This repository provides MATLAB simulation codes accompanying the paper:

> **Xiaoyu Li, Yuhu Wu, Lining Ru**, *Leader-following rendezvous control for generalized Cucker-Smale model on Riemannian manifolds*, SIAM Journal on Control and Optimization, 2024. DOI: [10.1137/23M1545811](https://doi.org/10.1137/23M1545811)

We consider a **double-integrator multi-agent system** moving on a Riemannian manifold.  
The generalized Cucker–Smale model describes the interaction among followers.  
A **feedback control law** is designed using **logarithmic maps** and **parallel transport** to achieve:

- **Position convergence**: all followers track a moving leader  
- **Velocity alignment**: follower velocities match the leader’s transported velocity  

The control is proven effective for:
- **Compact Riemannian manifolds** (e.g., sphere)
- **Flat Riemannian manifolds** (e.g., Euclidean space, circle, cylinder)

## Repository Structure

| File | Description |
|------|-------------|
| `Example4_5.m` | Unit sphere $\mathbb{S}^2$ |
| `Example4_8.m` | Euclidean plane $\mathbb{R}^2$ |
| `Example4_15.m` | Unit circle $\mathbb{S}^1$ |
| `Example4_25.m` | Infinite cylinder $\mathbb{S}^1\times\mathbb{R}$ |
