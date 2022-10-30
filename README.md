Theory from "Highly accurate prediction of material optical properties based on density functional theory" by Mitsutoshi Nishiwaki and Hiroyuki Fujiwara in 'Computational Materials Science 172 (2020) 109315'

Provides function that, given equispaced arrays of (1) energy, (2) / (3) real / imaginary parts of the
dielectric function, and a band gap shift (delta EG in the paper), returns the updated 'PHS' dielectric
functions on the same input array

Use:
	eps1PHS, eps2PHS = PHS(arrLen, energy, eps1, eps2, delta, energyTol=1.e-3)
'eps*PHS' are the blue/red shifted (by 'delta') dielectric response functions of 'eps1', 'eps2' computed on the EQUISPACED grid 'energy').  All arrays are of length 'arrLen'.  

Implementation detail that may be important:
The principal value integral is implemented as a modified trapezoid method that works outward from the pole in the integrand, summing the two competing components before adding the result to the trapezoid sum (this prevents a total loss of accuracy associated with densly sampled dielectric functions and the P.V. integral).
NOTE: the PHS dielectric functions will be inaccurate near the end-points of the sample arrays.  If you want the function accurate between 1 and 3 eV, for example, just send in samples for 0.5 to 3.5 eV.

