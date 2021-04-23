#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
#ifdef PTHREADS_NUM_THREADS
#include <pthread.h>
#endif

/* Routines for models that require stellar evolution: luminosities, mass loss, SNe rates, etc.
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */
#ifdef GALSF

#ifdef ADM
/* this routine tells the feedback algorithms what to 'inject' when a stellar feedback event occurs.
    you must define the mass, velocity (which defines the momentum and energy), and metal content (yields)
    of the ejecta for the event[s] of interest. Mass [Msne] and velocity [SNe_v_ejecta] should
    be in code units. yields[k] should be defined for all metal species [k], and in dimensionless units
    (mass fraction of the ejecta in that species). */
void particle2in_addFB_fromstars_adm(struct addFB_evaluate_data_in_ *in, int i, int fb_loop_iteration)
{
#ifdef METALS
    //int k; for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=0.178*All.SolarAbundances[k]/All.SolarAbundances[0];} // assume a universal solar-type yield with ~2.63 Msun of metals
    //if(NUM_LIVE_SPECIES_FOR_COOLTABLES>=10) {in->yields[1] = 0.4;} // (catch for Helium, which the above scaling would give bad values for)
#endif
#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)
    in->adm = P[i].adm;  // Only ADM particles should passed in. So make the input structure adm.
    if(P[i].SNe_ThisTimeStep<=0) {in->Msne=0; return;} // no event
    // 'dummy' example model. Just assumes ADM experiences no FB.
    in->Msne = 0; // No SNe physics
    in->SNe_v_ejecta = 0; // No SNe physics
    return;
#endif
}


/* this routine calculates the event rates for different types of mechanical/thermal feedback
    algorithms. things like SNe rates, which determine when energy/momentum/mass are injected, should go here.
    you can easily modify this to accomodate any form of thermal or mechanical feedback/injection of various
    quantities from stars. */
double mechanical_fb_calculate_eventrates_adm(int i, double dt)
{

#ifdef GALSF_FB_THERMAL /* STELLAR-POPULATION version: pure thermal feedback: assumes AGORA model (Kim et al., 2016 ApJ, 833, 202) where everything occurs at 5Myr exactly */
    if(P[i].SNe_ThisTimeStep != 0) {P[i].SNe_ThisTimeStep=-1; return 0;} // already had an event, so this particle is "done"
    if(evaluate_stellar_age_Gyr(P[i].StellarAge) < 0.005) {return 0;} // enforce age limit of 5 Myr
    P[i].SNe_ThisTimeStep = P[i].Mass*UNIT_MASS_IN_SOLAR / 91.; // 1 event per 91 solar masses
    return 1;
#endif

#ifdef GALSF_FB_MECHANICAL /* dummmy mechanical feedback model. Just assumes that no ADM SNe occur, so return nothing and make sure particle is not SNe_ThisTimeStep */
    double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
    if(star_age < 0.03)
    {
        //double RSNe = 3.e-4; // assume a constant rate ~ 3e-4 SNe/Myr/solar mass for t = 0-30 Myr //
        //double p = RSNe * (P[i].Mass*UNIT_MASS_IN_SOLAR) * (dt*UNIT_TIME_IN_MYR); // unit conversion factor
        //double n_sn_0=(float)floor(p); p-=n_sn_0; if(get_random_number(P[i].ID+6) < p) {n_sn_0++;} // determine if SNe occurs
        //P[i].SNe_ThisTimeStep = n_sn_0; // assign to particle
        // return RSNe;
        P[i].SNe_ThisTimeStep = 0; // assign to particle
        return 0;
    } else {
	P[i].SNe_ThisTimeStep = 0; // assign to particle
    }
#endif

    return 0;
}

#endif // ADM
#endif /* GALSF */
