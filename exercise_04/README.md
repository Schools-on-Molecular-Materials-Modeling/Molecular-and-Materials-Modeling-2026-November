# Exercise 04

Standard reaction enthalpy of the Diels–Alder $[4+2]$ cycloaddition of 1,3-butadiene and ethene to form cyclohexene, evaluated via MOPAC, xTB, and MACE, including comparative benchmark tables and a cyclopentadiene challenge.

### Overview & Benchmark
The repository provides benchmark comparisons against the experimental NIST gas-phase reaction enthalpy of **$-166.7\text{ kJ/mol}$** at $298.15\text{ K}$, detailing electronic energy changes, zero-point energy (ZPE) adjustments, and method-specific errors for MACE-OFF23, MOPAC PM7, and GFN2-xTB.

### Challenge
Calculate (MOPAC, xTB, MACE) the gas-phase Diels–Alder enthalpy for the reaction 
 1,3-butadiene + acrylonitrile $\rightarrow$ 4-cyanocyclohexene .

For the example reaction (cycloaddition of 1,3-butadiene and ethene to form cyclohexene) calculate the PySCF reaction enthalpy.

#### Thermochemical Reference Data (298.15 K, Gas Phase)

* **Reaction:** $\text{C}_4\text{H}_6(\text{g}) + \text{C}_3\text{H}_3\text{N}(\text{g}) \longrightarrow \text{C}_7\text{H}_9\text{N}(\text{g})$ (1,3-butadiene + acrylonitrile $\rightarrow$ 4-cyanocyclohexene)
* **Enthalpies of Formation ($\Delta_f H^\circ_{298.15\text{ K}}$):**
  * 1,3-Butadiene: $+110.0 \pm 0.8\text{ kJ/mol}$
  * Acrylonitrile: $+180.6 \pm 1.7\text{ kJ/mol}$
  * 4-Cyanocyclohexene: $+120.5 \pm 1.5\text{ kJ/mol}$ (NIST/Pedley)
* **Standard Reaction Enthalpy ($\Delta_r H^\circ_{298.15\text{ K}}$):** $\mathbf{-170.1}\text{ kJ/mol}$ ($-1.76\text{ eV}$, $-40.7\text{ kcal/mol}$)

