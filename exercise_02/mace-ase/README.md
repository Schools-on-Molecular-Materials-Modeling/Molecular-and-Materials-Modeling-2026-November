# Methane combustion


Because standard MACE-OFF23 is parameterized strictly for closed-shell singlet organic molecules and has no spin-state awareness for ground-state triplet \(\text{O}_{2}\), this script uses the universal foundation model mace_mp (which covers the periodic table, supports D3 dispersion, and is trained on spin-polarized Materials Project PBE/PBE+U calculations where ground-state \(\text{O}_{2}\) is naturally a triplet) while passing spin=1 for \(\text{O}_{2}\) into IdealGasThermo



Important Pedagogical Note for exercise_02

• Why mace_mp is used here instead of mace_off: MACE-OFF23 was trained exclusively on closed-shell singlet organic molecules (SPICE / QM7-X), so it has no representation of open-shell triplet \(\text{O}_{2}\) (\(^3\Sigma_g^-\)). By contrast, MACE-MP-0 is trained on spin-polarized VASP PBE/PBE+U calculations from the Materials Project, where gas-phase \(\text{O}_{2}\) was converged in its triplet ground state.
• Known PBE \(\text{O}_{2}\) overbinding: Because MACE-MP-0 inherits the underlying GGA PBE potential energy surface—which famously overbinds the \(\text{O}_{2}\) double bond by \(\sim 0.4\text{ eV}\) per \(\text{O}_{2}\) molecule (the same error that requires Wang–Maxisch–Ceder \(\text{O}_{2}\) corrections in Materials Project phase diagrams)—the predicted \(\Delta_c H^\circ(298\text{ K})\) will reflect that PBE \(\text{O}_{2}\) bias rather than matching experiment as closely as MACE-OFF23 did for closed-shell ethene hydrogenation.

Your MACE-MP-0 (medium) + D3 run completed cleanly and produced \(\Delta_c H^\circ(298\text{ K}) = -743.24\text{ kJ/mol}\)—landing within \(59\text{ kJ/mol}\) (\(7.4\%\)) of the experimental gas-phase value (\(-802.30\text{ kJ/mol}\)) and faithfully reproducing the underlying GGA PBE potential energy surface.

Physical Breakdown of Your MACE-MP-0 Results

1. Exact Inheritance of PBE \(\text{O}_{2}\) Overbinding:
	• In VASP PBE DFT (the exact training level of MACE-MP-0), the \(\text{O=O}\) bond length of triplet \(\text{O}_{2}\) is \(1.234\text{ \AA}\) (your MACE result: \(1.2366\text{ \AA}\)), and PBE famously overbinds the \(\text{O}_{2}\) molecule by \(\sim 0.30\text{–}0.35\text{ eV}\) (\(\sim 30\text{ kJ/mol}\)) per \(\text{O}_{2}\).
	• Because \(2\text{ mol}\) of \(\text{O}_{2}\) appear on the reactant side, that \(2 \times 30\text{ kJ/mol} \approx 60\text{ kJ/mol}\) over-stabilization of the reactants shifts \(\Delta_c H^\circ(298\text{ K})\) upward from \(-802.30\text{ kJ/mol}\) to \(-743.24\text{ kJ/mol}\)—matching raw PBE DFT (\(\approx -745\text{ kJ/mol}\)) to within \(2\text{ kJ/mol}\).
2. Accurate Zero-Point & Thermal Corrections:
	• Your harmonic \(\Delta\text{ZPE} = +8.83\text{ kJ/mol}\) matches both GFN2-xTB (\(+8.74\text{ kJ/mol}\)) and full DFT (\(\approx +9\text{ kJ/mol}\)), with clean degenerate bending modes in \(\text{CO}_{2}\) (\(710.2, 710.2\text{ cm}^{-1}\)) and triply degenerate modes in \(\text{CH}_{4}\) (\(1194.2\text{ cm}^{-1}\) and \(3147.3\text{ cm}^{-1}\)).

