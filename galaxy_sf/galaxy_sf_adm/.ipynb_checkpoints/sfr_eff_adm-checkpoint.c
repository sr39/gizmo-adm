#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_eigen.h>
#include "../../allvars.h"
#include "../../proto.h"

/*!
 *  routines for star formation in cosmological/galaxy/single-star/black hole simulations
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel,
 *   but the physical modules for star formation and feedback have been
 *   replaced, and the algorithm is mostly new to GIZMO. Many additional modules
 *   added since, with significant contributions from Mike Grudic.
 */


#ifdef GALSF // master switch for compiling the routines below //
#ifdef ADM


/* Routine to actually determine the SFR assigned to an individual gas particle at each time */
double get_starformation_rate_adm(int i)
{
    double rateOfSF,tsfr,y; y=0; int flag=1, j, k; /* flag to proceed to SFR calc */
    if(P[i].Mass <= 0 || SphP[i].Density <= 0) {flag=0;} /* zero-mass elements [for deletion] not eligible for SF */
#ifdef GALSF_SUBGRID_WINDS
    if(SphP[i].DelayTime > 0) {flag=0;} /* 'decoupled' wind elements not eligible for SF */
#endif
    if(All.ComovingIntegrationOn && SphP[i].Density < All.OverDensThresh) {flag=0;} /* below overdensity threshold required for SF */
    if(SphP[i].Density*All.cf_a3inv < All.PhysDensThresh) {flag=0;} /* below physical density threshold */
#if defined(GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED)
    if(flag==0) {SphP[i].AlphaVirial_SF_TimeSmoothed=0;} /* for time-smoothed virial param, reset to nil if fall below threshold */
#endif
    tsfr = sqrt(All.PhysDensThresh / (SphP[i].Density * All.cf_a3inv)) * All.MaxSfrTimescale; /* set default SFR timescale to scale appropriately with the gas dynamical time */
    rateOfSF = P[i].Mass / tsfr; /* 'normal' sfr from density law above */
    if(tsfr<=0 || rateOfSF<=0 || flag==0) {return 0;} /* nonsense here, return 0 */

#ifdef GALSF_EFFECTIVE_EQS /* do the SFR calc for the Springel-Hernquist EOS, before any 'fancy' sf criteria, when above-threshold, or else risk incorrect entropies */
    double factorEVP = pow(SphP[i].Density * All.cf_a3inv / All.PhysDensThresh, -0.8) * All.FactorEVP; /* evaporation factor */
    double egyhot = All.EgySpecSN / (1 + factorEVP) + All.EgySpecCold; /* specific energy of hot [volume-filling] phase gas */
    double tcool = GetCoolingTime_adm(egyhot, SphP[i].Density * All.cf_a3inv, SphP[i].Ne, i); /* cooling time of two-phase mix */
    y = tsfr / tcool * egyhot / (All.FactorSN * All.EgySpecSN - (1 - All.FactorSN) * All.EgySpecCold); /* parameter */
    double cloudmass_fraction = (1 + 1 / (2 * y) - sqrt(1 / y + 1 / (4 * y * y))); /* quasi-equilibrium mass in cold phase */
    rateOfSF = (1 - All.FactorSN) * cloudmass_fraction * P[i].Mass / tsfr; /* SFR given by cold mass (less SNe-entrainment fraction) divided by tSFR */
    update_internalenergy_for_galsf_effective_eos(i,tcool,tsfr,cloudmass_fraction,rateOfSF); /* updates entropies for the effective equation-of-state */
#endif

    int exceeds_force_softening_threshold; exceeds_force_softening_threshold = 0; /* flag that notes if the density is so high such that gravity is non-Keplerian [inside of smallest force-softening limits] */
#if (SINGLE_STAR_SINK_FORMATION & 1024)
    if(DMIN(PPP[i].Hsml, 2.*Get_Particle_Size(i)) <= DMAX(All.MinHsml, 2.*All.ForceSoftening[0])) {exceeds_force_softening_threshold=1;}
    if(exceeds_force_softening_threshold) {return 1.e4 * rateOfSF;}
#endif

    /* compute various velocity-gradient terms which are potentially used in the various criteria below */
    double dv2abs=0, divv=0, gradv[9]={0}, cs_eff=0, vA=0, v_fast=0; /* calculate local velocity dispersion (including hubble-flow correction) in physical units */
    cs_eff=Get_Gas_thermal_soundspeed_i(i); vA=Get_Gas_Alfven_speed_i(i); /* specifically get the -thermal- soundspeed and Alfven speed [dont include terms like radiation pressure or cosmic ray pressure in the relevant speeds here] */
    v_fast=sqrt(cs_eff*cs_eff + vA*vA); /* calculate fast magnetosonic speed for use below */
    for(j=0;j<3;j++) {
        for(k=0;k<3;k++) {
            double vt = SphP[i].Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
            if(All.ComovingIntegrationOn) {if(j==k) {vt += All.cf_hubble_a;}} /* add hubble-flow correction */
            gradv[3*j + k]=vt; dv2abs+=vt*vt; if(j==k) {divv+=vt;} // save for possible use below
        }}
#if defined(SINGLE_STAR_SINK_DYNAMICS) && defined(SINGLE_STAR_SINK_FORMATION) && ((defined(COOLING) && !defined(COOL_LOWTEMP_THIN_ONLY)) || defined(EOS_GMC_BAROTROPIC)) // if we have to deal with optically-thick thermo
    double nHcgs = HYDROGEN_MASSFRAC * (SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS);
    if(nHcgs > 1e13) {v_fast = DMIN(v_fast, 0.2/UNIT_VEL_IN_KMS);} // limiter to permit sink formation in simulations that really resolve the opacity limit and bog down when an optically-thick core forms. Modify this if you want to follow first collapse more/less - scale as c_s ~ n^(1/5)
#endif

#if (SINGLE_STAR_SINK_FORMATION & 1) || defined(GALSF_SFR_VIRIAL_SF_CRITERION) /* apply standard virial-parameter criteria here. note that our definition of virial parameter here is ratio of [Kinetic+Internal Energy]/|Gravitational Energy| -- so <1=bound, >1=unbound, 1/2=bound-and-virialized, etc. */
    double k_cs = M_PI * v_fast / (Get_Particle_Size(i)*All.cf_atime), alpha_crit; alpha_crit = 1.0; /* effective wavenumber for thermal+B-field+CR+whatever internal energy support, and threshold virial parameter */
    double Mach_eff_2=0, cs2_contrib=2.*k_cs*k_cs; Mach_eff_2=dv2abs/cs2_contrib; dv2abs+=2.*k_cs*k_cs; // account for thermal+magnetic pressure with standard Jeans criterion (k^2*cs^2 vs 4pi*G*rho) //
    double alpha_vir = dv2abs / (8.*M_PI * All.G * SphP[i].Density * All.cf_a3inv); // coefficient comes from different density profiles, assuming a constant velocity gradient tensor: 22.6=constant-density cube, 8pi[approximate]=constant-density sphere, e.g. rho~exp(-r^n) n={4,8,16,32,64}->{17.1,22.1,24.1,24.9,25.1,25.15} [approaches uniform-density sphere as n->infinity]
#if defined(GALSF_SFR_VIRIAL_CRITERION_TIMEAVERAGED) /* compute and prepare to use our time-rolling average virial criterion */
    double dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); /* the physical time-step */
    double alpha_0=1./(1.+alpha_vir), dtau=DMIN(1.,DMAX(0.,exp(-(DMIN(DMAX(8.*dtime/tsfr,0.),20.))))); /* dimensionless units for below */
    SphP[i].AlphaVirial_SF_TimeSmoothed = DMIN(DMAX(SphP[i].AlphaVirial_SF_TimeSmoothed * dtau + alpha_0 * (1.-dtau) , 1.e-10), 1.); /* update rolling time-averaged virial parameter */
    alpha_vir = 1./SphP[i].AlphaVirial_SF_TimeSmoothed - 1.; /* use the rolling average below */
#endif
    if(exceeds_force_softening_threshold) {alpha_vir /= 10.;} /* account for gravitational softening effects here, making this threshold less steep */
#if (GALSF_SFR_VIRIAL_SF_CRITERION <= 0) && !defined(GALSF_SFR_VIRIAL_CONTINUOUS_THOLD) && !(SINGLE_STAR_SINK_FORMATION & 512) /* 'weakest' mode: reduce [do not zero] SFR if above alpha_crit, and not -too- dense */
    if((alpha_vir>alpha_crit) && (SphP[i].Density*All.cf_a3inv<100.*All.PhysDensThresh)) {rateOfSF *= 0.0015;} /* PFH: note the 100x threshold limit here is an arbitrary choice currently set -by hand- to prevent runaway densities from this prescription! */
#endif
#if (GALSF_SFR_VIRIAL_SF_CRITERION > 0) && !defined(GALSF_SFR_VIRIAL_CONTINUOUS_THOLD) && !(SINGLE_STAR_SINK_FORMATION & 512) /* 'normal' mode: zero SF if don't meet virial threshold */
    if(alpha_vir>alpha_crit) {rateOfSF=0;} /* simple 'hard' threshold here */
#endif
#if (SINGLE_STAR_SINK_FORMATION & 512) || defined(GALSF_SFR_VIRIAL_CONTINUOUS_THOLD) /* semi-continuous SF as a function of alpha_vir */
    //rateOfSF *= exp(-1.4 * DMIN(DMIN(DMAX(sqrt(DMAX(MIN_REAL_NUMBER,alpha_vir)), 1.e-4), 1.e10),22.)); /* continuous cutoff of rateOfSF with increasing virial parameter as ~exp[-1.4*sqrt(alpha_vir)], following fitting function from Padoan 2012 [limit the values of sqrt(alpha_vir) here since we'll take an exponential so don't want a nan] */
    Mach_eff_2 = DMIN(DMAX(1.e-5, Mach_eff_2/3.), 1.e4); if(!isfinite(Mach_eff_2)) {Mach_eff_2=1.e4;}
    double S_ln=log(1.+Mach_eff_2/4.), S_crit=log(alpha_vir*(1.+2.*Mach_eff_2*Mach_eff_2/(1.+Mach_eff_2*Mach_eff_2))); // Mach_eff_2 is determined by the ratio of the kinetic to the thermal terms in the virial parameter, corrected to the 1D dispersion here
    rateOfSF *= 0.5 * exp(3.*S_ln/8.) * (1. + erf((S_ln-S_crit)/sqrt(2.*S_ln))); // multi-free-fall model, as in e.g. Federrath+Klessen 2012/2013 ApJ 761,156; 763,51 (similar to that implemented in e.g. Kretschmer+Teyssier 2020), based on the analytic models in Hopkins MNRAS 2013, 430 1653, with correct virial parameter [K+T used a definition which gives the wrong value for thermally-supported clouds]
#endif
#endif

#if (SINGLE_STAR_SINK_FORMATION & 256) || defined(GALSF_SFR_MOLECULAR_CRITERION) /* scale SFR to fraction of 'molecular' gas in cell */
    double ne=1, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred; // pull various known thermal properties, prepare to extract others //
    temperature = ThermalProperties_adm(u0, rho, i, &mu_meanwt, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties, like neutral fraction, temperature, etc, that we will use below //
    rateOfSF *= Get_Gas_Molecular_Mass_Fraction(i, temperature, nh0, ne, 0.);
#endif

#if (SINGLE_STAR_SINK_FORMATION & 2) || (GALSF_SFR_VIRIAL_SF_CRITERION >= 4) /* restrict to convergent flows */
    if(divv >= 0) {rateOfSF=0;} /* diverging flow, no SF */
#endif

#if (SINGLE_STAR_SINK_FORMATION & 128) || (GALSF_SFR_VIRIAL_SF_CRITERION >= 4) /* check that the velocity gradient is negative-definite, ie. converging along all principal axes, which is much stricter than div v < 0 */
    for(j=0;j<3;j++){ // symmetrize the velocity gradient
      for(k=0;k<j;k++){double temp = gradv[3*j + k]; gradv[3*j + k] = 0.5*(gradv[3*j + k] + gradv[3*k + j]); gradv[3*k + j] = 0.5*(temp + gradv[3*k + j]);}}
    gsl_matrix_view M = gsl_matrix_view_array(gradv, 3, 3); gsl_vector *eval1 = gsl_vector_alloc(3);
    gsl_eigen_symm_workspace *v = gsl_eigen_symm_alloc(3); gsl_eigen_symm(&M.matrix, eval1,  v);
    if(exceeds_force_softening_threshold==0) {for(k=0;k<3;k++) if(gsl_vector_get(eval1,k) >= 0) {rateOfSF=0;}} /* cannot apply this criterion when we exceed the limits where gravity is treated as fully-Newtonian, it will severely suppress 'true' collapse */
    gsl_eigen_symm_free(v); gsl_vector_free(eval1);
#endif

#ifdef GALSF_SFR_TIDAL_HILL_CRITERION /* check that the tidal tensor is negative-definite, ie. converging along all principal axes, indicating that we're dominating our environment gravitationally and are living in our own Hill sphere */
    if(exceeds_force_softening_threshold==0) {for(k=0;k<3;k++) {if(P[i].tidal_tensorps[k][k] >= 0) {rateOfSF=0;}}} /* we've already diagonized this in gravtree.c, so this is a straightforward check. again should only be applied where force calculation is fully-reliable */
#endif

#if (SINGLE_STAR_SINK_FORMATION & 64) || (GALSF_SFR_VIRIAL_SF_CRITERION >= 3) /* check if Jeans mass is low enough for conceivable formation of 'stars' */
    double cs_touse=cs_eff, MJ_crit=DMAX(DMIN(1.e3, 1.*P[i].Mass*UNIT_MASS_IN_SOLAR), 100.); /* for galaxy-scale SF, default to large ~1000 Msun threshold */
    if(exceeds_force_softening_threshold) {MJ_crit = DMAX(1.e4 , 10.*P[i].Mass*UNIT_MASS_IN_SOLAR);}
#ifdef SINGLE_STAR_SINK_FORMATION
    cs_touse=v_fast; MJ_crit=DMIN(1.e4, DMAX(1.e-3 , 100.*P[i].Mass*UNIT_MASS_IN_SOLAR)); /* for single-star formation use un-resolved Jeans mass criterion, with B+thermal pressure */
#endif
    double MJ_solar = 2.*pow(cs_touse*UNIT_VEL_IN_KMS/0.2,3)/sqrt(SphP[i].Density*All.cf_a3inv*UNIT_DENSITY_IN_NHCGS / (HYDROGEN_MASSFRAC*1.0e3));
    if(MJ_solar > MJ_crit) {rateOfSF=0;} /* if too massive Jeans mass, go no further */
#endif

#if (SINGLE_STAR_SINK_FORMATION & 4) /* restrict to local density/potential maxima */
    if(SphP[i].Density_Relative_Maximum_in_Kernel > 0) {rateOfSF=0;}
#endif

#if (SINGLE_STAR_SINK_FORMATION & 8) /* restrict to cell which neither 'sees' or 'is seen by' a sink too close */
    if(P[i].BH_Ngb_Flag) {rateOfSF=0;} /* cell cannot be 'seen' by -any- sink as a potential interacting neighbor */
    if(P[i].min_dist_to_bh < 1.24*Get_Particle_Size(i)) {rateOfSF=0;} /* cell does not see a sink within a volume = 8x=2^3 times its cell volume [set coefficient =1.86 for 27x=3^3 its cell volume] */
#endif

#if (SINGLE_STAR_SINK_FORMATION & 16) /* restrict to cells which have a local SF time or free-fall time shorter than their free-fall time onto the nearest sink */
    if(DMIN(P[i].min_bh_approach_time, P[i].min_bh_freefall_time) < tsfr) {rateOfSF=0;}
#endif



#if defined(SINGLE_STAR_SINK_DYNAMICS) && defined(SINGLE_STAR_SINK_FORMATION)
    rateOfSF *= 1.0e20; /* make sink formation guaranteed to happen, where it can, by setting rate super-high if non-zero */
#endif
    return rateOfSF; /* finally, we have a SFR! */
}



#endif // ADM
#endif // GALSF
