# Exercise 03

Hydrogenation enthalpy of ethene, computed with MOPAC, xTB, and MLIPs (ALIGNN-FF, GAP, MACE).

### Ethene Hydrogenation Reference

$$\text{C}_2\text{H}_4(\text{g}) + \text{H}_2(\text{g}) \longrightarrow \text{C}_2\text{H}_6(\text{g})$$

The standard gas-phase hydrogenation enthalpy of ethene at $298.15\text{ K}$ is **$-136.3\text{ kJ/mol}$** (**$-1.41\text{ eV}$**, or **$-32.6\text{ kcal/mol}$**).

#### Thermochemical Breakdown

Calculated from standard enthalpies of formation ($\Delta_f H^\circ_{298.15\text{ K}}$) via NIST data:

- **Ethene ($\text{C}_2\text{H}_4$, g):** $+52.4 \pm 0.5\text{ kJ/mol}$
- **Hydrogen ($\text{H}_2$, g):** $0.00\text{ kJ/mol}$
- **Ethane ($\text{C}_2\text{H}_6$, g):** $-83.8 \pm 0.4\text{ kJ/mol}$

Subtracting the reactants from the product gives:

$$\Delta_r H^\circ(298.15\text{ K}) = -83.8 - (+52.4 + 0.00) = -136.2\text{ to }-136.3\text{ kJ/mol}$$



## Challenge

With the help of AI agents, adapt the scripts for **propene**:

$$\text{C}_3\text{H}_6(\text{g}) + \text{H}_2(\text{g}) \longrightarrow \text{C}_3\text{H}_8(\text{g})$$

The standard gas-phase hydrogenation enthalpy of propene at $298.15\text{ K}$ is **$-125.1\text{ kJ/mol}$** (**$-1.30\text{ eV}$**, or commonly rounded to **$-124$ to $-125\text{ kJ/mol}$**).

### Thermochemical Breakdown

Calculated from standard enthalpies of formation ($\Delta_f H^\circ_{298.15\text{ K}}$) via NIST data:

- **Propene ($\text{C}_3\text{H}_6$, g):** $+20.41\text{ kJ/mol}$
- **Hydrogen ($\text{H}_2$, g):** $0.00\text{ kJ/mol}$
- **Propane ($\text{C}_3\text{H}_8$, g):** $-104.7 \pm 0.5\text{ kJ/mol}$

This yields a reaction enthalpy of **$-125.11\text{ kJ/mol}$**, supported by early calorimetric data adjusted to $298.15\text{ K}$ ($\approx -124.2$ to $-125.1\text{ kJ/mol}$).

