#Theory from:
# "Highly accurate prediction of material optical properties based on density functional theory"
#by:
# Mitsutoshi Nishiwaki, Hiroyuki Fujiwara
#in:
# Computational Materials Science 172 (2020) 109315

#Provides functions that, given equispaced arrays of (1) energy, (2) / (3) real / imaginary parts of the
#dielectric function, and a band gap shift (delta EG in the paper), returns the updated 'PHS' dielectric
#functions on the same input array

#The principal value integral is implemented as a modified trapezoid method that works outward from the
#pole in the integrand, summing the two competing components before adding the result to the trapezoid sum
#(this prevents a total loss of accuracy associated with densly sampled dielectric functions and the P.V.
#integral).
#NOTE: this means that the PHS dielectric functions will be inaccurate near the end-points of the sample
#arrays.  If you want the function accurate between 1 and 3 eV, for example, just send in samples for 0.5
#to 3.5 eV.

#Integrand - helps below function look like less of a disaster
def _f_aux_phs(ep, e, eps2PHS_ep):
    return ep*eps2PHS_ep / (ep*ep - e*e)


#Returns eps1', eps2' in that order
def PHS(arrLen, energy, eps1, eps2, delta, energyTol=1.e-3):
    dE = energy[1] - energy[0]

    ##check to make sure that the energy is actually evenly spaced
    for i in range(2, arrLen):
        if (abs((energy[i] - energy[i - 1]) - dE) > energyTol):
            print("PHS(): energy sample spacings seem to be inconsistent!")
            return None, None


    ##Easy part: shift imaginary portion by delta
    ## ... well, it was simple until I took negative deltas into account ...
    di = int(round(delta/dE)) ##Set difference in index from E to deltaE
    eps2PHS = [0. for i in range(0, arrLen)]
    if(di >= 0):
        ###Main, exact portion of the algo
        for i in range(di, arrLen):
            eps2PHS[i] = (energy[i] - delta)/(energy[i] + energyTol) * eps2[i - di]
        ###Approximate portion that the algo can't deal with
        for i in range(0, di):
            eps2PHS[i] = eps2PHS[di]
    else:
        ###Main, exact portion of the algo
        for i in range(0, arrLen + di):
            eps2PHS[i] = (energy[i] - delta)/(energy[i] + energyTol) * eps2[i - di]
        ###Approximate portion that the algo can't deal with
        for i in range(arrLen + di, arrLen):
            eps2PHS[i] = eps2PHS[arrLen + di - 1]

    ##Hard part: PV integral for real portion
    ##Note that the multiplicitive factors for the function evaluations in the t-rule are always 2 here:
    ##the dielectric function is 0 at 0 eV and infinite eV, so there's no need to account for the 1's.  At
    ##the singularity, the function is +(infty) + (-infty) = 0, so we can ignore the 1's there too.

    ##TODO (if it ever actually matters ...):
    ##technically, itr should start between 1 and some index offset determined by an epsilon value related to
    ##NEDOS (higher NEDOS = closer evals to the pole needs a higher epsilon)
    tmp = 0.
    eps1PHS = [0. for i in range(0, arrLen)]
    for sampleInd in range(0, arrLen):
        for itr in range(1, max(sampleInd, arrLen - sampleInd)):
            if(sampleInd + itr >= arrLen): ##extend rhs to 0
                eps1PHS[sampleInd] += 2.*_f_aux_phs(ep=energy[sampleInd - itr],
                                                    e=energy[sampleInd],
                                                    eps2PHS_ep=eps2PHS[sampleInd - itr])
                continue

            if(sampleInd - itr < 0): ##extend lhs to 0
                eps1PHS[sampleInd] += 2.*_f_aux_phs(ep=energy[sampleInd + itr],
                                                    e=energy[sampleInd],
                                                    eps2PHS_ep=eps2PHS[sampleInd + itr])
                continue

            tmp = 2.*_f_aux_phs(ep=energy[sampleInd - itr], e=energy[sampleInd],
                                eps2PHS_ep=eps2PHS[sampleInd - itr]) + \
                  2.*_f_aux_phs(ep=energy[sampleInd + itr], e=energy[sampleInd],
                                eps2PHS_ep=eps2PHS[sampleInd + itr])
            eps1PHS[sampleInd] += tmp

        eps1PHS[sampleInd] *= 2./3.141592654 * dE/2. ##second part from the trap method
        eps1PHS[sampleInd] += 1.


    return eps1PHS, eps2PHS
