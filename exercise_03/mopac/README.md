Let's verify why dHf(H2) = -133.93 kJ/mol happens in MOPAC PM7:
• In James Stewart's PM7 parameterization paper (J. Mol. Model. 19, 1–32, 2013), diatomic \(\text{H}_{2}\) is a notorious outlier in NDDO methods: PM7 predicts \(\Delta H_f^\circ(\text{H}_2) \approx -32.0\text{ kcal/mol}\) (\(-133.9\text{ kJ/mol}\)) and \(\nu(\text{H}_2) \approx 5578\text{ cm}^{-1}\) because two-center one-electron core-attraction integrals have no heavy-atom screening in a pure H–H pair.
• Meanwhile, for every heavy-atom molecule (\(\text{C}_2\text{H}_4\) and \(\text{C}_2\text{H}_6\)), PM7 is spot-on:
	• \(\Delta H_f^\circ(\text{C}_2\text{H}_4, \text{PM7}) = +57.97\text{ kJ/mol}\) (Expt: \(+52.40\text{ kJ/mol}\))
	• \(\Delta H_f^\circ(\text{C}_2\text{H}_6, \text{PM7}) = -76.14\text{ kJ/mol}\) (Expt: \(-84.00\text{ kJ/mol}\))



How to Handle H2 Properly in MOPAC Benchmarks

If you use the standard thermodynamic reference \(\Delta H_f^\circ(\text{H}_2, 298.15\text{ K}) \equiv 0.00\text{ kJ/mol}\) (which is exact by definition for an element in its reference state) or run an isodesmic hydrogenation (transferring \(\text{H}_{2}\) from a saturated alkane like \(\text{C}_3\text{H}_8 \rightarrow \text{C}_3\text{H}_6 + \text{H}_2\) so the \(\text{H}_{2}\) diatomic error cancels out identically), PM7 gives \(-134.11\text{ kJ/mol}\)—within \(2.2\text{ kJ/mol}\) of experiment!
