This repository contains the code accompanying the paper:

**“Data-Adaptive Estimation in Time-to-Event Data with Competing Risks: A Super Learner Analysis of Alzheimer’s Disease Data from Clalit EHR.”**

---

## Overview

This project provides an implementation of **continuous-time** and **discrete-time Super Learner methods** for competing risks. It includes:

- The main workflow for running the Super Learner  
- Functions required for both continuous-time and discrete-time settings  
- Simulation code for reproducing results from the paper and supplementary materials  
- A case study using publicly available RHC data (since Clalit EHR data cannot be shared)

---

## Repository Structure

### `main.R`
This script demonstrates how to use the main Super Learner functions for both continuous-time and discrete-time competing risks analysis. It includes:

- Required R packages  
- Detailed explanation of input formats  
- Description of output objects  
- Example usage

---

### `function/`

This folder contains all functions needed for the project:

- Continuous-time Super Learner functions  
- Discrete-time Super Learner functions  
- Simulation utilities, including:
  - Data-generating mechanisms  
  - Plotting functions  

---

### `simulation.R`

This script reproduces the simulation results presented in the paper and supplementary materials.  
You can flexibly switch between different data-generating settings by modifying the annotated sections in the script. The comparison methods are also included.

---

### Case Study: RHC Data

Because the Clalit EHR dataset used in the paper is not publicly available, we include an alternative case study using the **public RHC dataset**.  
We also compare our method with the **causalCmprsk** package to illustrate performance differences.
Since the RHC dataset is much smaller than the Clalit EHR dataset, we did not use data splitting to construct confidence intervals.

---


