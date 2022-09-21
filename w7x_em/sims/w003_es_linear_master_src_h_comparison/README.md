### An initial set of simulations to benchmark ES stella (src h vs master)
- These are relatively low-resolution
- These were constructed using the src. h sims from the fprim tprim ky scans as the starting point
- The master sims used the same physical & non-physical parameters
- gradients were set to 3. ky=0.5
- Unfortunately, it looks like these non-physical parameters are insufficient to resolve the physics; all sims disagree
- For reference, the same case, but with ky=0.452, has (freq, gamma) = (-0.17365153E+00  0.51264636E-01)
## How the simulations perform
- master, explicit: (freq, gamma) = (0.31718741E+03  0.81150892E+01). Numerical instability?
- master, str + mirror implicit: (freq, gamma) = (-0.10463090E+03  0.88375408E+00). Numerical instability?
- src. h, explicit: (freq, gamma) = (0.41830835E+00  0.13323020E-01). Suspiciously stable
- src. h, str + mirror implicit: (freq, gamma) = (-0.10167023E+03  0.25218895E+01) (but gamma unconverged)
