2D e-AWBLM
---------------------------------------------------------------------------------
# e-AWBLM

**Enhanced Atmospheric Wave Boundary Layer Model for physically based air–sea momentum exchange under strong winds and finite water depth**

## Overview

**e-AWBLM** is an open-source implementation of the **enhanced atmospheric wave boundary layer model** developed for estimating **wind stress over the ocean surface** under conditions that are especially important for **storm surge**, **storm wave**, and **compound flood** modeling.

This repository builds on the atmospheric wave boundary layer framework of **Xu and Yu (2021)** and further incorporates the **effect of wave breaking on air–sea momentum exchange**, following **Zhang and Yu (2024)**. The model is designed to improve wind-stress estimation not only in deep water, but also in **shallow and intermediate-depth coastal waters**, where conventional bulk drag formulations are often inadequate.

In our broader research, e-AWBLM serves as a key physical module for simulating **typhoon-induced storm surges, storm waves, and compound coastal–fluvial flooding**, particularly in large delta regions such as the **Pearl River Delta**.

---

## Why this repository matters

Accurate air–sea momentum exchange is one of the most fundamental requirements in coastal hazard modeling. During tropical cyclones and other extreme events, the ocean surface is strongly affected by:

- rapidly evolving wave states,
- shallow-water effects,
- wave breaking,
- nonlinear atmosphere–wave–ocean interactions.

These processes directly influence the **surface wind stress**, which is the main momentum input for both **circulation models** and **wave models**. If the wind stress is poorly represented, errors can propagate into:

- storm surge magnitude,
- wave setup and radiation stress effects,
- coastal water levels,
- inundation extent,
- and ultimately flood-risk assessment.

e-AWBLM is developed to address this problem using a **physics-based formulation** rather than relying solely on empirical drag parameterizations.

---

## Scientific motivation

Most traditional wind-drag formulations used in coastal and ocean models are simple and efficient, but they may become unreliable under **extreme wind conditions** or over **finite-depth waters**, where wave development and wave breaking substantially modify the atmospheric boundary layer.

The goal of e-AWBLM is to provide a more physically consistent estimate of wind stress by explicitly accounting for:

1. **wave-state dependence**,
2. **finite water depth effects**,
3. **breaking-induced stress**, and
4. the coupled behavior of the **atmospheric boundary layer** and the **wave field**.

This is particularly important for studies that aim to move beyond event hindcasts and toward **mechanism-based hazard assessment**, including the simulation of:

- historical typhoon storm surges,
- storm-wave–surge interactions,
- and compound flood processes in river deltas under changing tropical cyclone conditions.

---

## Connection to our research

This repository is not just a standalone drag-coefficient calculator. It is a core physical component in a broader modeling framework developed in our research on coastal and compound flooding.

### 1. Improved wind stress for storm surge and wave modeling
The enhanced model was developed to improve the estimation of air–sea momentum transfer under strong winds, with special attention to the role of **wave breaking** and **finite depth**.

### 2. Coupling with SWAN + ADCIRC
e-AWBLM can be coupled with **SWAN** and **ADCIRC** to improve simulations of **storm waves** and **storm surges** in coastal regions.

### 3. Support for compound flood modeling
In our coupled **land–river–ocean** modeling framework, e-AWBLM is used to provide more realistic ocean-side forcing for simulating **tropical-cyclone-induced compound floods**, where storm surges, waves, tides, river discharge, and rainfall-runoff interact.

### 4. Foundation for future hazard assessment
Because it is physics-based, e-AWBLM is especially useful for studies that need a robust representation of ocean surface forcing when linking **tropical cyclone changes** to **future coastal and compound flood risk**.

---

## What is included in this repository

The repository currently contains the core e-AWBLM implementation and example scripts for reproducing representative calculations.

### Main folders

- `v1.0/`  
  Original implementation used to reproduce the representative wind-stress-coefficient results.

- `v2.0/`  
  Updated and faster implementation with reduced computational cost.

- `function/`  
  Supporting MATLAB functions required by the model.

- `data/`  
  Precomputed or supporting data files used by the example scripts.

- `ADCIRC+SWAN/`  
  Materials related to model coupling and application within ocean circulation and wave simulations.

---

## Key features

- Physics-based estimation of **ocean-surface wind stress**
- Applicable to **finite-depth** as well as deep-water conditions
- Includes **wave-breaking effects** on air–sea momentum exchange
- Designed for use in **extreme-wind coastal environments**
- Suitable for coupling with **storm surge** and **wave** models
- Useful for both **process-based research** and **hazard modeling applications**

---

## Conceptual framework

The model is based on conservation of momentum and energy within the **atmospheric wave boundary layer** over the ocean surface. Compared with simpler drag laws, e-AWBLM explicitly links wind stress to the evolving wave field and its breaking behavior.

### Model concept
![Figure 1](https://github.com/anyifang/e-AWBLM/assets/89235013/c6d6a062-bffd-4721-adff-fbc6e1f06bde)

*Figure 1. Conceptual sketch of the enhanced atmospheric wave boundary layer model.*

### Representative wind stress coefficients
![Figure 3](https://github.com/anyifang/e-AWBLM/assets/89235013/b2a6a174-3be6-4fb1-a660-f93ee3f234f3)

*Figure 2. Representative wind stress coefficients under different water depths.*

### Coupled application with SWAN + ADCIRC
![Figure 8](https://github.com/user-attachments/assets/649e64ad-7cfa-41d7-96aa-a41662d5508e)

*Figure 3. Wind stress coefficients obtained from e-AWBLM coupled with SWAN + ADCIRC in shallow, intermediate, and deep water.*

---

## Typical applications

e-AWBLM is intended for applications such as:

- **storm surge simulation**
- **storm wave simulation**
- **typhoon / hurricane coastal hazard analysis**
- **air–sea interaction studies under extreme winds**
- **coupled wave–circulation modeling**
- **compound coastal–fluvial flood modeling**
- **future hazard projection using physically based coastal forcing**

---

## Quick start

### Requirements

- **MATLAB**
- Required input data placed in the repository `data/` folder
- Required function files placed in the repository `function/` folder

### Run the original implementation

Open MATLAB, change the working directory to:

Reference

Xu, Y., and X. Yu, 2021: Enhanced atmospheric wave boundary layer model for evaluation of wind stress over waters of finite depth. Prog. Oceanogr., 198, 102664, https://doi.org/10.1016/j.pocean.2021.102664.

Zhang, A. and Yu, X., 2024. A Major Improvement of Atmospheric Wave Boundary Layer Model for Storm Surge Modeling by Including Effect of Wave Breaking on Air–Sea Momentum Exchange. Journal of Physical Oceanography, 54(5), pp.1153-1168.

