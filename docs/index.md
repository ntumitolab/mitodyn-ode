# Glucose stimulation and Mitochondrial dynamics

## Difference from the original Fridlyand's model

- We  added an dilution term for pyruvate and a ADP Michaelis constant in Glyceraldehyde 3-phosphate dehydrogenase to increase the robustness in the ODE systems.
- We added adenylate kinase (Adk) equilibrium for the ATP-ADP-AMP pool to address the high energy phosphate distribution in the cytosol. Also, the adenylate pool size was bumped up to keep the steady state ATP/ADP ratio the same as before (~5) at the resting glucose level of 5mM.
- We altered the mathematical expressions of mitochondrial sodium calcium exchanger (NCLX). The  specific concentrations of sodium and calcium were on the exponent in the original model. We moved the specific concentrations down to multipliers was inconsistent with the original model of Nguyen et al.
- The mitochondrial fission rate is fixed to once per 10 minutes, and the fusion rate depends on the proton leak flux and the ATP synthase flux.

## Current Limitations

- Glycolysis, ETC, ATP synthesis, and cytosolic calcium dynamics were expressed in an empirical fashion rather than fully mechanistic. No signal transduction nor complex calcium dynamics.
- There is no role for reactive oxygen species (ROS),  potential signals for bioenergetics and mitochondrial dynamics.
- Mitochondrial content was assumed conserved and homogenous.

## General parameters

| Parameter     | Value         | Description                                                           |
| ------------- | ------------- | --------------------------------------------------------------------- |
| $V_{i}$       | 0.53          | Relative cytoplasmic volume                                           |
| $V_{m}$       | 0.06          | Relative mitochondrial volume                                         |
| $V_{mtx}$     | 0.0144        | Relative mitochondrial matrix volume (Adjustable)                     |
| $C_{mito}$    | 1.812 mM/V    | Mitochondrial membrane capacitance                                    |
| $F$           | 96484.6 C/mol | Faraday's constant                                                    |
| $\delta_{Ca}$ | 0.0003        | Fraction of free Ca in mitochondria                                   |
| $[Na^+]_c$    | 10 mM         | Cytoplasmic Na concentration                                          |
| $[Na^+]_m$    | 5 mM          | Mitochondrial Na concentration                                        |
| $T_v$         | 26.73 mV      | Thermal voltage (RT/F) (@37°C)                                        |
| $\Sigma A_c$  | 4.5 mM        | Cellular adenine nucleotides concentration  (Adjustable)              |
| $\Sigma N_m$  | 2.2 mM        | Free pyridine nucleotides concentration in<br/>mitochondrial matrix   |
| $\Sigma N_c$  | 2.0 mM        | Free pyridine nucleotides concentration in<br/>cytoplasm (Adjustable) |
| $k_{gpd}$     | 0.01/s        | Consumption rate of G3P                                               |
| $k_{NADHm}$   | 0.1/s         | Consumption rate of mito NADH                                         |
| $k_{NADHc}$   | 0.1/s         | Consumption rate of cyto NADH                                         |
| $k_{ATP}$     | 0.04/s        | Basal consumption rate of ATP                                         |
| $k_{ATPCa}$   | 90/mM/s       | Consumption rate of ATP activated by calcium                          |


## Conservation relationships

$$
\begin{aligned}
\Sigma N_{m} &= [NAD^+]_m + [NADH]_m         \\
\Sigma N_{c} &= [NAD^+]_c + [NADH]_c         \\
\Sigma A_{c} &= [ATP]_c + [ADP]_c + [AMP]_c  \\
1 &= X_1 + 2X_2 + 3X_3
\end{aligned}
$$

## Adenylate kinase

$$
\begin{aligned}
J_{ADK} &= k_f ([ADP]_c^2 - [ATP]_c [AMP]_c / K_{eq}^{AK})
\end{aligned}
$$

| Parameter     | Value     | Description                                                                                                                              |
| :------------ | --------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| $k_f$         | 1000 mM/s | Forward (AMP-forming) rate constant of adenylate kinase. <br />The parameters was set arbitrary large for equilibrium of adenylate pool. |
| $K_{eq}^{AK}$ | 0.931     | Equilibrium constant of adenylate kinase (AMP-forming).                                                                                  |

## Glucokinase (GK)

$$
J_{glu} = V_{m} \frac{[ATP]_c}{[ATP]_c + K_{ATP}} \frac{[Glc]^n}{[Glc]^n + K_{Glc}^{n}}
$$

| Parameter | Value      | Description                    |
| --------- | ---------- | ------------------------------ |
| $V_{m}$   | 0.011 mM/s | Max rate of glucokinase        |
| $K_{ATP}$ | 0.5 mM     | Michaelis constant for ATP     |
| $K_{Glc}$ | 7 mM       | Michaelis constant for glucose |
| n         | 1.7        | Cooperativity for glucose      |


## Glyceraldehyde 3-phosphate dehydrogenase (GPD)

$$
J_{gpd} = V_m \frac{[G3P]}{[G3P] + K_{G3P}} \frac{[NAD^+]_c}{[NAD^+]_c + K_{NAD}[NADH]_c}
$$

| Parameter | Value    | Description                                      |
| --------- | -------- | ------------------------------------------------ |
| $V_{m}$   | 0.5 mM/s | Max rate of GPD (Adjustable)                     |
| $K_{G3P}$ | 0.2 mM   | Michaelis constant for G3P                       |
| $K_{NAD}$ | 0.09     | Activation constant for cytosolic NAD/NADH ratio |

## Lactate production by lactate dehydrogenase (LDH)

The rate of lactate output is approximately 5% of the rate of glucose consumption when the glucose level is 8 mM.

$$
J_{LDH} = V_m \frac{[Pyr]}{[Pyr] + K_{Pyr}} \frac{[NADH]_c}{[NADH]_c + K_{NADH} [NAD^+]_c}
$$

| Parameter  | Value     | Description                                      |
| ---------- | --------- | ------------------------------------------------ |
| $V_{m}$    | 1.2 mM/s  | Max rate of LDH (Adjustable)                     |
| $K_{Pyr}$  | 0.0475 mM | Michaelis constant for pyruvate                  |
| $K_{NADH}$ | 1         | Activation constant for cytosolic NADH/NAD ratio |


## Steady-state cytosolic calcium levels

$$
[Ca^{2+}]_c = [Ca^{2+}]_R + k_{A}^{Ca} \frac{([ATP]_c)^n}{([ATP]_c)^n + (K_{ATP} [ADP]_c)^n}
$$

| Parameter     | Value  | Description                               |
| ------------- | ------ | ----------------------------------------- |
| $[Ca^{2+}]_R$ | 90 nM  | Resting cytoplasmic calcium concentration |
| $k_{A}^{Ca}$  | 250 nM | Maximal activated calcium concentration   |
| $K_{ATP}$     | 25     | Activation constant for ATP/ADP ratio     |
| $n$           | 4      | Cooperativity for ATP/ADP ratio           |

## Oscillating calcium levels

Use in simulations for Fig. 4 : oscillating calcium on mitochondrial bioenergetics and dynamics only.

$$
\begin{aligned}
[Ca^{2+}]_c &= [Ca^{2+}]_R + k_{A}^{Ca}(Axe^{1-Ax})^B \\
x           &= \frac{t}{T} - \lfloor \frac{t}{T} \rfloor
\end{aligned}
$$

| Parameter     | Value    | Description                               |
| ------------- | -------- | ----------------------------------------- |
| $[Ca^{2+}]_R$ | 90 nM    | Resting cytoplasmic calcium concentration |
| $k_{A}^{Ca}$  | 250 nM   | Maximal activated calcium concentration   |
| $A$           | 5        | Asymmetric factor                         |
| $B$           | 4        | Steepness factor                          |
| $T$           | 2 minute | Period of calcium oscillations.           |


## Pyruvate dehydrogenase (PDH)

We assume that pyruvate diffuses freely and fast across the inner mitochondrial membrane (IMM). Therefore, pyruvate is the same concentration in the cytosol and in the mitochondrial matrix.

$$
\begin{aligned}
J_{PDH} &= V_m \frac{ [Pyr] }{ [Pyr] + K_{Pyr}} \frac{[NAD^+]_m}{[NAD^+]_m + K_{NAD}} \frac{(1 + C)^2}{(1 + C)^2 (1 + u_2) + u_2 u_1} \\
C &= [Ca^{2+}]_m / K_{Ca}
\end{aligned}
$$

| Parameter | Value     | Description                                          |
| --------- | --------- | ---------------------------------------------------- |
| $V_{m}$   | 0.3 mM/s  | Max rate of PDH                                      |
| $K_{Pyr}$ | 0.0475 mM | Michaelis constant for pyruvate                      |
| $K_{NAD}$ | 81        | Activation constant for mitochondrial NAD/NADH ratio |
| $K_{Ca}$  | 50 nM     | Activation constant for mitochondrial Ca             |
| $u_1$     | 1.5       | Factor for calcium activation                        |
| $u_2$     | 1.1       | Factor for calcium activation                        |

## Electron transport chain (ETC)

$$
\begin{aligned}
J_{hr} &= V_m \frac{[NADH]_m}{[NADH]_m + K_{NADH}} \frac{1 + k_A \Delta \Psi_m}{1 + k_B \Delta \Psi_m} F_{O_2}
\end{aligned}
$$

| Parameter  | Value       | Description                    |
| ---------- | ----------- | ------------------------------ |
| $V_{m}$    | 22 mM/s     | Max rate of ETC                |
| $K_{NADH}$ | 3 mM        | Michaelis constant for NADH    |
| $k_A$      | -4.92 /Volt | thermodynamic potential factor |
| $k_B$      | -4.43 /Volt | thermodynamic potential factor |
| $F_{O_2}$  | 1           | Oxygen availability            |

## F1Fo ATPase (ATP synthase)

ATP synthase was lumped with ANT and depended on cytosolic ADP.

$$
\begin{aligned}
J_{hf} &= V_m f_{ADP} f_{\Psi} f_{Ca} \\
J_{ANT} &= J_{hf} / H_{ATP} \\
f_{ADP} &= \frac{[MgADP]_c^{n_A}}{[MgADP]_c^{n_A} + K_{ADP}^{n_A}}   \\
f_{\Psi} &= \frac{\Delta\Psi_m^{n_\Psi}}{\Delta\Psi_m^{n_\Psi} + K_{\Psi}^{n_\Psi}}   \\
f_{Ca} &=  1 - \exp(-[Ca^{2+}]_m / K_{Ca} )    \\
\text{[MgADP]}_c &= 0.055[ADP]_c
\end{aligned}
$$

| Parameter  | Value    | Description                                       |
| ---------- | -------- | ------------------------------------------------- |
| $V_{m}$    | 8 mM/s   | Max rate of ATP synthase (Adjustable)             |
| $K_{ADP}$  | 20 μM    | Apparent Michaelis constant for cytosolic MgADP   |
| $n_A$      | 2        | Cooperativity for MgADP                           |
| $n_\Psi$   | 8        | Cooperativity for mitochondrial potential         |
| $K_{\Psi}$ | 131.4 mV | Mid-activity constant for mitochondrial potential |
| $K_{Ca}$   | 0.165 μM | Activation constant for mitochondrial calcium     |
| $H_{ATP}$  | 3        | H:ATP ratio                                       |

## Proton leak

The basal leak approaches ~20% of the electron transport rate at $\Delta\Psi_m$ of 160 mV.

$$
J_{hl}  = P_{H}\exp(k_{lp} \Delta \Psi_m)
$$

| Parameter | Value       | Description                    |
| --------- | ----------- | ------------------------------ |
| $P_{H}$   | 0.0024 mM/s | leak coefficient               |
| $k_{lp}$  | 30.5/V      | membrane potential coefficient |

## NADH shuttles

$$
J_{TNADH} = T_{NADH} \frac{[NADH]_c}{[NADH]_c + [NAD^+]_c K_c} \frac{[NAD^+]_m}{[NAD^+]_m + [NADH]_m K_m}
$$

| Parameter  | Value     | Description                                      |
| ---------- | --------- | ------------------------------------------------ |
| $T_{NADH}$ | 0.05 mM/s | NADH transport rate                              |
| $K_c$      | 0.002     | Affinity coefficients for cytoplasmic NADH/NAD   |
| $K_m$      | 16.78     | Affinity coefficients for mitochondrial NAD/NADH |

## Mitochondrial calcium uniporter (MCU)

$$
\begin{aligned}
J_{uni} &= P_{Ca} \frac{\delta}{e^{\delta} - 1} (e^{\delta} \alpha_i [Ca^{2+}]_c - \alpha_m [Ca^{2+}]_m)  \\
\delta &= Z_{Ca}\Delta\Psi_m / V_T
\end{aligned}
$$

| Parameter  | Value | Description                       |
| ---------- | ----- | --------------------------------- |
| $P_{Ca}$   | 4 / s | Permeability of calcium           |
| $Z_{Ca}$   | 2     | Valence of calcium                |
| $\alpha_i$ | 0.341 | Activity of cytoplasmic calcium   |
| $\alpha_m$ | 0.2   | Activity of mitochondrial calcium |

## Mitochondrial Sodium-Calcium exchanger (NCLX)

We used the electron-neutral descriptor of NCLX since this model generated smooth and monotonous increment of mitochondrial calcium levels upon increasing glucose levels.

$$
\begin{aligned}
J_{NCLX} &= V_{m} (AB - PQ ) / D \\
D &= 1 + A + B + P + Q + AB + PQ \\
A &= ([Na^+]_c / K_{Na})^2  \\
B &= [Ca^{2+}]_m / K_{Ca}  \\
P &= ([Na^+]_m / K_{Na})^2  \\
Q &= [Ca^{2+}]_c / K_{Ca}  \\
\end{aligned}
$$

| Parameter | Value      | Description                 |
| --------- | ---------- | --------------------------- |
| $V_m$     | 0.075 mM/s | Max rate of NCLX            |
| $K_{Ca}$  | 8 μM       | Dissociation constant of Ca |
| $K_{Na}$  | 8.2 mM     | Dissociation constant of Na |

## Mitochondrial Dynamics

$$
\begin{align}
k_{fuse, 1} &= k_{0}^{fuse} \frac{ J_{ANT} }{ J_{HL} } \\
k_{fiss, 1} &= k_{0}^{fiss}                            \\
k_{fuse, 2} &= 0.1k_{fuse, 1}                          \\
k_{fiss, 2} &= 1.5k_{fiss, 1}                          \\
\end{align}
$$

| Parameter    | Value              | Description            |
| ------------ | ------------------ | ---------------------- |
| $k_0^{fuss}$ | $\frac{1}{600}$ Hz | The basal fusion rate  |
| $k_0^{fiss}$ | $\frac{1}{600}$ Hz | The basal fission rate |

## Ordinary differential equations

$$
\begin{aligned}
\frac{d}{dt}[G3P] &= \frac{1}{V_i}(2J_{glu} - J_{GPD}) - k_{g3p}[G3P] \\
\frac{d}{dt}[Pyr] &= \frac{1}{V_i + V_{mtx}}(J_{GPD} - J_{LDH} - J_{PDH}) - k_{pyr}[Pyr]   \\
\frac{d}{dt}[NADH]_c &= \frac{1}{V_i}(J_{GPD} - J_{LDH} - J_{NADHT}) - k_{nadhc}[NADH]_c   \\
\frac{d}{dt}[NADH]_m &= \frac{1}{V_{mtx}}(J_{NADHT} + 4.6J_{PDH} - 0.1J_{hr}) - k_{nadhm}[NADH]_m \\
\frac{d}{dt}\Delta\Psi_m &= \frac{1}{C_{mito}} (J_{hr} - J_{hf} - J_{ANT} - J_{hl} - 2J_{uni}) \\
\frac{d}{dt}[Ca^{2+}]_m &= \frac{f_m}{V_{mtx}} (J_{uni} - J_{NCLX}) \\
\frac{d}{dt}[ATP]_c &= \frac{1}{V_i} (J_{ANT} - 2J_{Glu} + 2J_{GPD} + J_{AdK}) - (k_{ATP} + k_{CaATP}[Ca^{2+}]_c)[ATP]_c \\
\frac{d}{dt}[ADP]_c &= -\frac{d}{dt}[ATP]_c - \frac{J_{AdK}}{V_i}  \\
\frac{d}{dt} X_2 &= k_{fuse, 1} X_1^2 - (k_{fiss, 1} + k_{fuse, 2}X_1) X_2 + k_{fiss, 2}X_3 \\
\frac{d}{dt} X_3 &= k_{fuse, 2}X_1X_2- k_{fiss, 2}X_3
\end{aligned}
$$


## Initial conditions

| State variable   | Value   | Description                                |
| ---------------- | ------- | ------------------------------------------ |
| $[G3P]$          | 2.8μM   | Glyceraldehyde-3-phosphate                 |
| $[Pyr]$          | 8.5μM   | Pyruvate                                   |
| $[NADH]_c$       | 1μM     | Cytosolic NADH                             |
| $[NADH]_{m}$     | 60μM    | Mitochondrial NADH                         |
| $[ATP]_c$        | 4mM     | Cytosolic ATP concentration                |
| $[ADP]_c$        | 0.5mM   | Cytosolic ADP concentration                |
| $[Ca^{2+}]_{m}$  | 0.250μM | Mitochondrial calcium concentration        |
| $\Delta\Psi_{m}$ | 100mV   | Mitochondrial membrane potential           |
| $X_2$            | 0.20    | Population of degree-2 mitochondrial nodes |
| $X_3$            | 0.05    | Population of degree-3 mitochondrial nodes |
