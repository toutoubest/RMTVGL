# RMTVGL  
Robust and Missing-Data-Aware Time-Varying Graphical Lasso

RMTVGL implements a robust extension of the Time-Varying Graphical Lasso (TVGL) tailored for high-dimensional temporal data with **missing values and heavy-tailed contamination**.  
The method integrates:

- **EM-based missing data imputation**,  
- **Huber loss for robustness against outliers**, and  
- **temporal smoothness regularization** to capture evolving network structures.

It is designed for applications such as dynamic gene regulatory network inference, where data are noisy, incomplete, and time-indexed.

This repository accompanies the methodology introduced in my paper:
> *“Robust and Missing-Data-Aware Time-Varying Graphical Lasso (RM-TVGL) via EM-Based Adaptive Estimation and Nonconvex Regularization” (2025).*


