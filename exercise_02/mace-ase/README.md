# Methane combustion


Because standard MACE-OFF23 is parameterized strictly for closed-shell singlet organic molecules and has no spin-state awareness for ground-state triplet \(\text{O}_{2}\), this script uses the universal foundation model mace_mp (which covers the periodic table, supports D3 dispersion, and is trained on spin-polarized Materials Project PBE/PBE+U calculations where ground-state \(\text{O}_{2}\) is naturally a triplet) while passing spin=1 for \(\text{O}_{2}\) into IdealGasThermo



Important Pedagogical Note for exercise_02

• Why mace_mp is used here instead of mace_off: MACE-OFF23 was trained exclusively on closed-shell singlet organic molecules (SPICE / QM7-X), so it has no representation of open-shell triplet \(\text{O}_{2}\) (\(^3\Sigma_g^-\)). By contrast, MACE-MP-0 is trained on spin-polarized VASP PBE/PBE+U calculations from the Materials Project, where gas-phase \(\text{O}_{2}\) was converged in its triplet ground state.
• Known PBE \(\text{O}_{2}\) overbinding: Because MACE-MP-0 inherits the underlying GGA PBE potential energy surface—which famously overbinds the \(\text{O}_{2}\) double bond by \(\sim 0.4\text{ eV}\) per \(\text{O}_{2}\) molecule (the same error that requires Wang–Maxisch–Ceder \(\text{O}_{2}\) corrections in Materials Project phase diagrams)—the predicted \(\Delta_c H^\circ(298\text{ K})\) will reflect that PBE \(\text{O}_{2}\) bias rather than matching experiment as closely as MACE-OFF23 did for closed-shell ethene hydrogenation.

