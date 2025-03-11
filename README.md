# Characterization of Proton Beam Parameters using Monitor Cross Sections

<details>
  <summary><strong>Table of Contents</strong></summary>

- [Introduction](#introduction)
- [Experimental Setup](#experimental-setup)
- [Experimental Challenges](#experimental-challenges-and-considerations)
  - [Beam Alignment](#beam-alignment)
  - [Collimator](#collimator)
  - [Gamma Spectrometry](#gamma-spectrometry)
  - [Target Material](#target-materials)
  - [Data Analysis](#data-analysis)
  - [General Algorithm](#general-algorithm)
  - [Optimization Algoritms](#optimization-algorithms)

</details>

## Introduction

This repository conveys the struggle to characterize the energy and current of IP2 proton beam by using IAEA-recommended monitor cross sections. Monitor cross sections are powerful because it allows one to assess beam parameters, e.g., proton current and energy, by irradiating targets of specific materials without the control of all the irradiation parameters. The produced activity is determined with techniques, e.g., gamma or alpha spectrometry, liquid scintillation counting (LSC), etc. Considering the the role played by the cross section in determining the produced activity:

```math
\frac{dN}{dt}=\sigma\Phi N_0 - \lambda\cdot N
```

Where $N$ is the surface nuclide density ($cm^{-2}$), $\sigma$ the monitor cross section ($mbarn$), $\Phi=\frac{I_p}{A}$ the proton flux ($cm^{-2}s^{-1}$), $N_0=\frac{N_A}{MW}\cdot \frac{dm}{dA}$ the initial surface nuclide density, $\lambda$ the decay constant, $t$ the irradiation time. Solving the above equation:

```math
A(t)=\lambda\cdot N(t) = \sigma\Phi N_0 \cdot(1-e^{-\lambda \cdot t})
```

Assuming an homogeneous and uniform mass surface distribution for the irradiated targets, leading to $\frac{dm}{dA}\approx\frac{m}{A}$ with a mass thickness $\rho_m = \rho\cdot\Delta t$, the end-of-bombardment (EoB) equation becomes:

```math
A(t,E_p)=\frac{\sigma(E_p) I_p N_A \rho_m}{MW}(1-e^{-\lambda \cdot t})
```

To determine proton current and energy, we can use two approaches:

1. Use monitor cross sections of a single production reaction considering both parameters to be fitted at the same time, leading to a 2D fitting problem. For this, we could use $Ti(p,x)^{48}V$ and $Ni(p,x)^{57}Ni$ reactions to compare results obtained using two monitor cross sections. Moreover, the produced radionuclides possess few gamma rays with high intensities and energy, thereby facilitating measurements using gamma spectrometry.
1. Use the ratio of activities of two radionuclides co-produced in the same target material, e.g., $^{62}Zn$ and $^{63}Zn$ or $^{65}Zn$. In this case, we can get rid of the current dependency. To obtain low-uncertainty results, it would be advisable to exploit combinations of radionuclides which exhibit a cross section ratio function with steep gradients in the investigated energy range. Considering the half-lives of zinc radionuclides, only $^{63}Zn$ ($t_{1/2}=38\, min$) might constitute a limitation due to transfer time from the irradiation station to the HPGe for spectrometry measurements. Moreover, the ratio between $^{62}Zn$ and the other two useful zinc radioisotopes exhibit a low gradient, which might reflect in higher uncertainty in the fitted proton energy. Then, the following equation can be used:

```math
\frac{A_1}{A_2}=\frac{\sigma_1(E_p)}{\sigma_2(E_p)}\cdot\frac{(1-e^{-\lambda_1t})}{(1-e^{-\lambda_2t})}
```

## Experimental Setup

The intial idea consists in using the same aluminum capsule used for production tests at IP2, Paul Scherrer Institute, to host a stack of foils 6 mm in diameter, replicating the same irradiation conditions during production irradiations. Its pocket should be deep enough to contain the whole stack of (Ti)+(Ni)+(Energy degrader) to degrade the proton beam from 23 MeV to 13 MeV. It is good practice in the community using stacked foil technique to degrader the proton energy beam at maximum to half of the beam energy at the first foil. Aluminum foils should be suitable according to its stopping power and the energy range investigated, accounting for thickness and costs. The stack should be irradiated for short times and with low proton currents to avoid producing too high activities, resulting in long cooling times to prevent exposure to high radiation doses.

## Experimental Challenges and Considerations

### Beam alignment

If we proceed irradiating stacks of foils 6 mm in diameter, then we should ensure the beam is well aligned, meaning it is as much as possible perpendicular to the target surface plane. If this assumption is not satisfied, we risk to assume that all the foils in the stack are irradiated by the same proton flux. In this case two sources of error should be accounted for:
    - The beam has non uniform surface distribution. Therefore, if the beam direction is also tilted with respect to the normal of the target surface, the last foil in the stack will be irradiated by a different proton current.
    - The protons will travel a longer distance $l'=l/\cos(\theta)$, where $\theta$ is the angle formed between the normal to the target lid plane and the impingin direction of the proton beam. Therefore, the simulations for energy degradation will be erroneous.

### Collimator

To ensure no protons are scattered into the foil stack from the aluminum capsule, a niobium foil with a central hole of 6 mm in diameter concentric with respect to the foil center could be included in the target capsule. This prevents to have proton that traversed different materials other than those of the stack. Considering the target materials we want to utilize, they all exhibit higher stopping power than aluminum, which composes the target capsule. According to a SRIM simulation, a 200-micrometer thick niobium foild would be sufficient to completely stop a 23 $\pm$ 2 MeV proton beam.

### Gamma spectrometry

The assessment of the EoB activities produced should not represent a challenge if we are using titanium and nickel as target materials and we irradiate with a 23 MeV proton beam. Both $^{48}V$ and $^{57}Ni$ have high-energy gamma rays with high absolute intensities. Considering copper targets for zinc radioisotope production, the main challenge is represented by $^{63}Zn$ due to its short half-life, $t_{1/2}=38.5\,min$. However, if the irradiations and sample transfer from the target station to the location of the HPGe are properly planned, this should not constitute any problem. Preliminary evaluations of EoB activities for $^{63}Zn$ should be performed to assess whether it is possible to determine its activity with low uncertainty due to its half-life and low-intesity gamma rays: $E_{\gamma,1}=669.6\,keV$, $I_{\gamma,1}=8.3\,%$ and $E_{\gamma,2}=963.1\,keV$, $I_{\gamma,2}=6.5\,%$. Nonetheless, we could use $^{62}Zn$-$^{65}Zn$ pair to determine $E_p$ by means of Method 2. The energy of all the gamma rays of each possible radionuclide is >500 keV. Therefore, the efficiency calibration of the HPGe detector using $^{152}Eu$ calibration source is sufficient to obtain a low-uncertainty fitting of the analytical function. A source of interferences with these radionuclides gamma lines might arise from the presence of chemical impurities in target materials.

### Target materials

Thorough characterization of target materials should be performed in advance to determine the potential presence of impurities. The knowledge of these impurities is relevant in the context of gamma spectrometry, since activation products might emit gamma rays that could interfere with the characteristic gamma emissions of the radionuclides we are interested in.

### Data analysis

The project requires to experimentally determine the EoB activities $A_m$ of the radionuclides whose cross section is recommended by the International Atomic Energy Agency (IAEA). This value is, then, compared to EoB activies evaluated using the recommended monitor cross sections $A_e$, considering $\bar{E_p}$ and $\bar{I_p}$ as estimators of the proton energy and current, respectively. The analysis algorithm will focus on minimizing the following *Objective Function*:

```math
Q=\sum \frac{(A_{m,i} - A_{e,i})^2}{\sigma_i^2}
```

where $\sigma_i^2$ is the variance associated to $A_{m,i}$.\
The *Least Square Method Estimation* (LSME) cannot be employed as fitting method because it assumes that the function is __directly evaluable__, __differentiable__, or at least expressible in a form where an anlytical or numerical least-squares approach can be applied. There are two ways to evaluate the energy at each foil of the stack:

1. __Stopping power models__, e.g., Bethe-Bloch equation with finite differences.
1. __Monte Carlo simulations__, e.g. SRIM, Geant4, Fluka, etc.

Bethe-Bloch equation derived using the quantum mechanical treatment of the particle collision for charged hadrons:

```math
-\langle \frac{dE}{dx}\rangle=K\frac{Z}{A}\rho\frac{z^2}{\beta^2}\left[ \frac{1}{2}\ln\frac{2m_e c^2\beta^2\gamma^2T_{max}}{I^2} -\beta^2 - \frac{\delta(\beta\gamma)}{2} - \frac{C(\beta\gamma,I)}{Z} \right]
```

where $K=4\pi N_A r_e^2 m_e c^2= 0.307\, MeV cm^2/mol$ and $r_e=e^2/4\pi \epsilon_0 m_e c^2\approx 2.8 fm$, the classical radius of the electron, are __fundamental constants__; $z,\beta,\gamma$ are the charge (in units of $e$), the velocity, and the relativistic factor, respectively, which are properties of the __incident particle__; $Z$ atomic number, $A$ the atomic mass (g/mol), $I$ the average energy required to ionize the medium, $T_{max}$ the maximum energy in a head-on collision, $\delta$ the density correction (relevant for relativistic motion), $C/Z$ is an additional 'shell' correction for low values of $\beta$. However, the use of this equation still requires to use a discrete approach based on finite differences to determine the energy degradation of the projectiles in matter.

### General algorithm

The algorithm to fit the data in order to obtain the esitmators' values is the following:

1. Use *Beam Delivery SIMulation* (BDSIM), a Geant4-based program, to determine the initial guesses for the estimators.
1. Use *Stopping and Range of Ions in Matter* (SRIM) or *BDSIM* to simulate the proton energy degradation through the foil stack to obtain the proton energy at each foil.
1. Compute the EoB activity of the radionuclides in each foil.
1. Minimize the objective function $Q$ by means of an optmization algorithm to find the estimators' values.

### Optimization algorithms

The EoB acitivity is differentiable, however, the input proton energy must be computed for each foil iteratively. For this reason, optimization algorithm should satisfy the following requirements:

- __Handle black-box evaluations__. We have to compute MC simulations or finite differences.
- __Not require explicit derivatives__. We cannot derive an analytical expression for the total activity.
- __Efficiency__. The algorithm should be computationally efficient to find the optimum solution through the least number of evaluation steps.

The following methods could be potential candidates to solve our problem:

- __Adaptive Grid Search__.
- __Differential Evolution (DE) Optmization__.
- __Bayesian Optimization__.

| __Method__                | __Pros__                                                                 | __Cons__                                                        | __Best Use Case__                                             |
|---------------------------|------------------------------------------------------------------------|----------------------------------------------------------------|--------------------------------------------------------------|
| __Adaptive Grid Search__   | ✅ Simple to implement ✅ Guaranteed to find the global minimum <br> ✅ No need for derivatives | ❌ Computationally expensive if function evaluations are slow <br> ❌ Requires manually setting refinement levels | ✔️ When function evaluations (e.g., SRIM) are __fast__ <br> ✔️ When the objective function is __well-behaved__ |
| __Differential Evolution__ (DE) | ✅ Global optimization (avoids local minima) <br> ✅ No need for gradients <br> ✅ Works well with noisy data | ❌ Requires many function evaluations <br> ❌ Slower convergence than gradient-based methods | ✔️ When function evaluations are __moderate in cost__ <br> ✔️ When the function is __complex and non-differentiable__ |
| __Bayesian Optimization__  | ✅ Minimizes number of function evaluations <br> ✅ Provides uncertainty quantification <br> ✅ Works well for expensive black-box functions | ❌ Requires probabilistic modeling (Gaussian Process) <br> ❌ Computationally expensive in high dimensions | ✔️ When function evaluations (e.g., Monte Carlo simulations) are __very expensive__ <br> ✔️ When __function shape is unknown__ |
