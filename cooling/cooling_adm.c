#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"
#include "./cooling.h"


/*
#######################
######## EDITS ########
#######################
This file will be separate from cooling.c It will contain all the necessary functions
for the cooling routine for a ADM particle. However, it will not go through
any particles. The actual function call will be done by cooling.c. This file only
contain the ADM cooling functions needed to cool the ADM particles. That's it.

Note: As of 09 April 2021, I have not included modifications for the metal line
cooling. I have not enabled those features for ADM
*/


/*
 * This file contains the routines for optically-thin cooling (generally aimed towards simulations of the ISM,
 *   galaxy formation, and cosmology). A wide range of heating/cooling processes are included, including
 *   free-free, metal-line, Compton, collisional, photo-ionization and recombination, and more. Some of these
 *   are controlled by individual modules that need to be enabled or disabled explicitly.
 *
 * This file was originally part of the GADGET3 code developed by Volker Springel. The code has been modified heavily by
 *   Phil Hopkins (phopkins@caltech.edu) for GIZMO; essentially everything has been re-written at this point */


#ifdef COOLING
#ifdef ADM
/* these are variables of the cooling tables. they are static but this shouldnt be a problem for shared-memory structure because
    they are only defined once in a global operation, then locked for particle-by-particle operations */
/* requires the cooling table TREECOOL, which is included in the GIZMO source in the cooling directory */
#define NCOOLTAB_ADM  2000 /* defines size of cooling table */

static double Tmin_adm = -1.0, Tmax_adm = 9.0, deltaT_adm; /* minimum/maximum temp, in log10(T/K) and temperature gridding: will be appropriately set in make_cooling_tables subroutine below */
static double *BetaH0_adm, *BetaHep_adm, *Betaff_adm, *AlphaHp_adm, *AlphaHep_adm, *Alphad_adm, *AlphaHepp_adm, *GammaeH0_adm, *GammaeHe0_adm, *GammaeHep_adm; // UV background parameters
#ifdef COOL_METAL_LINES_BY_SPECIES
/* if this is enabled, the cooling table files should be in a folder named 'spcool_tables' in the run directory.
 cooling tables can be downloaded at: http://www.tapir.caltech.edu/~phopkins/public/spcool_tables.tgz or on the Bitbucket site (downloads section) */
/* ############# THESE NEED TO BE ALTERED FOR METAL LINE COOLING!!! ###########
   NEED TO EDIT SpCoolTable */
static float *SpCoolTable0, *SpCoolTable1;
#endif
/* these are constants of the UV background at a given redshift: they are interpolated from TREECOOL but then not modified particle-by-particle */
static double J_UV_ADM = 0, gJH0_adm = 0, gJHep_adm = 0, gJHe0_adm = 0, epsH0_adm = 0, epsHep_adm = 0, epsHe0_adm = 0;


/* subroutine which actually sends the particle data to the cooling routine and updates the entropies */
void do_the_cooling_for_particle_adm(int i)
{
    double unew, dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);

    if((dtime>0)&&(P[i].Mass>0)&&(P[i].Type==0))  // upon start-up, need to protect against dt==0 //
    {
#ifdef COOL_MOLECFRAC_NONEQM
        update_explicit_molecular_fraction_adm(i, 0.5*dtime*UNIT_TIME_IN_CGS); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif
        double uold = DMAX(All.MinEgySpec, SphP[i].InternalEnergy);

#ifndef COOLING_OPERATOR_SPLIT
        /* do some prep operations on the hydro-step determined heating/cooling rates before passing to the cooling subroutine */
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
        /* calculate the contribution to the energy change from the mass fluxes in the gravitation field */
        double grav_acc; int k;
        for(k = 0; k < 3; k++)
        {
            grav_acc = All.cf_a2inv * P[i].GravAccel[k];
#ifdef PMGRID
            grav_acc += All.cf_a2inv * P[i].GravPM[k];
#endif
            SphP[i].DtInternalEnergy -= SphP[i].GravWorkTerm[k] * All.cf_atime * grav_acc;
        }
#endif
        /* limit the magnitude of the hydro dtinternalenergy */
        SphP[i].DtInternalEnergy = DMAX(SphP[i].DtInternalEnergy , -0.99*SphP[i].InternalEnergy/dtime ); // equivalent to saying this wouldn't lower internal energy to below 1% in one timestep
        SphP[i].DtInternalEnergy = DMIN(SphP[i].DtInternalEnergy ,  1.e8*SphP[i].InternalEnergy/dtime ); // equivalent to saying we cant massively enhance internal energy in a single timestep from the hydro work terms: should be big, since just numerical [shocks are real!]
        /* and convert to cgs before use in the cooling sub-routine */
        SphP[i].DtInternalEnergy *= (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (PROTONMASS/HYDROGEN_MASSFRAC);
#endif


#ifndef RT_COOLING_PHOTOHEATING_OLDFORMAT
        /* Call the actual COOLING subroutine! */
        unew = DoCooling_adm(uold, SphP[i].Density * All.cf_a3inv, dtime, SphP[i].Ne, i);
#else
        unew = uold + dtime * (rt_DoHeating(i, dtime) + rt_DoCooling(i, dtime));
#endif


#if defined(BH_THERMALFEEDBACK)
        if(SphP[i].Injected_BH_Energy) {unew += SphP[i].Injected_BH_Energy / P[i].Mass; SphP[i].Injected_BH_Energy = 0;}
#endif



#ifdef RT_INFRARED /* assume (for now) that all radiated/absorbed energy comes from the IR bin [not really correct, this should just be the dust term] */
        double nHcgs = HYDROGEN_MASSFRAC * UNIT_DENSITY_IN_CGS * SphP[i].Density / PROTONMASS;	/* hydrogen number dens in cgs units */
        double ratefact = nHcgs * nHcgs / (SphP[i].Density * UNIT_DENSITY_IN_CGS);
        double de_u = -SphP[i].LambdaDust * ratefact * dtime*UNIT_TIME_IN_CGS / (UNIT_SPECEGY_IN_CGS) * P[i].Mass; /* energy gained by gas needs to be subtracted from radiation */
        if(de_u<=-0.99*SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED]) {de_u=-0.99*SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED]; unew=DMAX(0.01*SphP[i].InternalEnergy , SphP[i].InternalEnergy-de_u/P[i].Mass);}
        SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED] += de_u; /* energy gained by gas is lost here */
        SphP[i].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED] = SphP[i].Rad_E_gamma[RT_FREQ_BIN_INFRARED]; /* updated drifted */
#if defined(RT_EVOLVE_INTENSITIES)
        int k_tmp; for(k_tmp=0;k_tmp<N_RT_INTENSITY_BINS;k_tmp++) {SphP[i].Rad_Intensity[RT_FREQ_BIN_INFRARED][k_tmp] += de_u/RT_INTENSITY_BINS_DOMEGA; SphP[i].Rad_Intensity_Pred[RT_FREQ_BIN_INFRARED][k_tmp] += de_u/RT_INTENSITY_BINS_DOMEGA;}
#endif
        double momfac = CRSOL_OVER_CTRUE_SQUARED_FOR_BEAMING * de_u / All.cf_atime; int kv; // add leading-order relativistic corrections here, accounting for gas motion in the addition/subtraction to the flux
#if defined(RT_EVOLVE_FLUX)
        for(kv=0;kv<3;kv++) {SphP[i].Rad_Flux[RT_FREQ_BIN_INFRARED][kv] += momfac*SphP[i].VelPred[kv]; SphP[i].Rad_Flux_Pred[RT_FREQ_BIN_INFRARED][kv] += momfac*SphP[i].VelPred[kv];}
#endif
        momfac = 1. - CRSOL_OVER_CTRUE_SQUARED_FOR_BEAMING * de_u / (P[i].Mass * C_LIGHT_CODE_REDUCED*C_LIGHT_CODE_REDUCED); // back-reaction on gas from emission
        for(kv=0;kv<3;kv++) {P[i].Vel[kv] *= momfac; SphP[i].VelPred[kv] *= momfac;}
#endif


        /* InternalEnergy, InternalEnergyPred, Pressure, ne are now immediately updated; however, if COOLING_OPERATOR_SPLIT
         is set, then DtInternalEnergy carries information from the hydro loop which is only half-stepped here, so is -not- updated.
         if the flag is not set (default), then the full hydro-heating is accounted for in the cooling loop, so it should be re-zeroed here */
        SphP[i].InternalEnergy = unew;
        SphP[i].InternalEnergyPred = SphP[i].InternalEnergy;
        SphP[i].Pressure = get_pressure(i);
#ifndef COOLING_OPERATOR_SPLIT
        SphP[i].DtInternalEnergy = 0;
#endif

#ifdef COOL_MOLECFRAC_NONEQM
        update_explicit_molecular_fraction_adm(i, 0.5*dtime*UNIT_TIME_IN_CGS); // if we're doing the H2 explicitly with this particular model, we update it in two half-steps before and after the main cooling step
#endif

    } // closes if((dt>0)&&(P[i].Mass>0)&&(P[i].Type==0)) check
}




/* returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
double DoCooling_adm(double u_old, double rho, double dt, double ne_guess, int target)
{
    double u, du; u=0; du=0;

#ifdef COOL_GRACKLE
#ifndef COOLING_OPERATOR_SPLIT
    /* because grackle uses a pre-defined set of libraries, we can't properly incorporate the hydro heating
     into the cooling subroutine. instead, we will use the approximate treatment below to split the step */
    du = dt * SphP[target].DtInternalEnergy / ( (UNIT_SPECEGY_IN_CGS/UNIT_TIME_IN_CGS) * (PROTONMASS/HYDROGEN_MASSFRAC));
    u_old += 0.5*du;
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
    /* now we attempt to correct for what the solution would have been if we had included the remaining half-step heating
     term in the full implicit solution. The term "r" below represents the exact solution if the cooling function has
     the form d(u-u0)/dt ~ -a*(u-u0)  around some u0 which is close to the "ufinal" returned by the cooling routine,
     to which we then add the heating term from hydro and compute the solution over a full timestep */
    double r=u/u_old; if(r>1) {r=1/r;} if(fabs(r-1)>1.e-4) {r=(r-1)/log(r);} r=DMAX(0,DMIN(r,1));
    du *= 0.5*r; if(du<-0.5*u) {du=-0.5*u;} u+=du;
#else
    /* with full operator splitting we just call grackle normally. note this is usually fine,
     but can lead to artificial noise at high densities and low temperatures, especially if something
     like artificial pressure (but not temperature) floors are used such that the temperature gets
     'contaminated' by the pressure terms */
    u = CallGrackle(u_old, rho, dt, ne_guess, target, 0);
#endif
    return DMAX(u,All.MinEgySpec);
#endif


    int iter=0, iter_upper=0, iter_lower=0, iter_condition = 0; double LambdaNet, ratefact, u_upper, u_lower;
#ifdef RT_INFRARED
    double LambdaDust;
#endif
    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u_old *= UNIT_SPECEGY_IN_CGS;
    dt *= UNIT_TIME_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho;

    u = u_old; u_lower = u; u_upper = u; /* initialize values */
    LambdaNet = CoolingRateFromU_adm(u, rho, ne_guess, target);

    /* bracketing */
    if(u - u_old - ratefact * LambdaNet * dt < 0)	/* heating */
    {
        u_upper *= sqrt(1.1); u_lower /= sqrt(1.1);
        while((iter_upper<MAXITER)&&(u_upper - u_old - ratefact * CoolingRateFromU_adm(u_upper, rho, ne_guess, target) * dt < 0))
        {
            u_upper *= 1.1; u_lower *= 1.1; iter_upper++;
        }

    }

    if(u - u_old - ratefact * LambdaNet * dt > 0) /* cooling */
    {
        u_lower /= sqrt(1.1); u_upper *= sqrt(1.1);
        while((iter_lower<MAXITER)&&(u_lower - u_old - ratefact * CoolingRateFromU_adm(u_lower, rho, ne_guess, target) * dt > 0))
        {
            u_upper /= 1.1; u_lower /= 1.1; iter_lower++;
        }
    }

    /* core iteration to convergence */
    do
    {
        u = 0.5 * (u_lower + u_upper);
#ifdef RT_INFRARED
        LambdaDust = SphP[target].LambdaDust;
#endif
        LambdaNet = CoolingRateFromU_adm(u, rho, ne_guess, target);
        if(u - u_old - ratefact * LambdaNet * dt > 0) {u_upper = u;} else {u_lower = u;}
        du = u_upper - u_lower;
        iter++;
        if(iter >= (MAXITER - 10)) {printf("u=%g u_old=%g u_upper=%g u_lower=%g ne_guess=%g dt=%g iter=%d \n", u,u_old,u_upper,u_lower,ne_guess,dt,iter);}

        iter_condition = ((fabs(du/u) > 3.0e-2)||((fabs(du/u) > 3.0e-4)&&(iter < 10)));
#ifdef RT_INFRARED
        iter_condition = iter_condition || (((fabs(LambdaDust - SphP[target].LambdaDust) > 1e-2*fabs(LambdaDust)) || (fabs(u - u_old - ratefact * LambdaNet * dt) > 0.01*fabs(u-u_old)))  && (iter < MAXITER-11));
#endif
        iter_condition = iter_condition &&  (iter < MAXITER); // make sure we don't iterate more than MAXITER times

    }
    while(iter_condition); /* iteration condition */
    /* crash condition */
    if(iter >= MAXITER) {printf("failed to converge in DoCooling(): u_in=%g rho_in=%g dt=%g ne_in=%g target=%d \n",u_old,rho,dt,ne_guess,target); endrun(10);}
    double specific_energy_codeunits_toreturn = u / UNIT_SPECEGY_IN_CGS;    /* in internal units */

#ifdef RT_CHEM_PHOTOION
    /* set variables used by RT routines; this must be set only -outside- of iteration, since this is the key chemistry update */
    double u_in=specific_energy_codeunits_toreturn, rho_in=SphP[target].Density*All.cf_a3inv, mu=1, temp, ne=1, nHI=SphP[target].HI, nHII=SphP[target].HII, nHeI=1, nHeII=0, nHeIII=0;
    temp = ThermalProperties_adm(u_in, rho_in, target, &mu, &ne, &nHI, &nHII, &nHeI, &nHeII, &nHeIII);
    SphP[target].HI = nHI; SphP[target].HII = nHII;
#ifdef RT_CHEM_PHOTOION_HE
    SphP[target].HeI = nHeI; SphP[target].HeII = nHeII; SphP[target].HeIII = nHeIII;
#endif
#endif

    /* safe return */
    return specific_energy_codeunits_toreturn;
}



/* returns cooling time.
 * NOTE: If we actually have heating, a cooling time of 0 is returned.
 */
double GetCoolingTime_adm(double u_old, double rho, double ne_guess, int target)
{
#if defined(COOL_GRACKLE) && !defined(GALSF_EFFECTIVE_EQS)
    double LambdaNet = CallGrackle(u_old, rho, 0.0, ne_guess, target, 1);
    if(LambdaNet >= 0) LambdaNet = 0.0;
    return LambdaNet / UNIT_TIME_IN_CGS;
#else
    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u_old *= UNIT_SPECEGY_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    double LambdaNet = CoolingRateFromU_adm(u_old, rho, ne_guess, target);
    if(LambdaNet >= 0) {return 0;} /* net heating due to UV background */
    return u_old / (-(nHcgs * nHcgs / rho) * LambdaNet) / UNIT_TIME_IN_CGS;
#endif
}


/* returns new internal energy per unit mass.
 * Arguments are passed in code units, density is proper density.
 */
double DoInstabilityCooling_adm(double m_old, double u, double rho, double dt, double fac, double ne_guess, int target)
{
    if(fac <= 0) {return 0.01*m_old;} /* the hot phase is actually colder than the cold reservoir! */
    double m, dm, m_lower, m_upper, ratefact, LambdaNet;
    int iter = 0;

    rho *= UNIT_DENSITY_IN_CGS;	/* convert to physical cgs units */
    u *= UNIT_SPECEGY_IN_CGS;
    dt *= UNIT_TIME_IN_CGS;
    fac /= UNIT_SPECEGY_IN_CGS;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    ratefact = nHcgs * nHcgs / rho * fac;
    m = m_old; m_lower = m; m_upper = m;
    LambdaNet = CoolingRateFromU_adm(u, rho, ne_guess, target);

    /* bracketing */
    if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt < 0)	/* heating */
    {
        m_upper *= sqrt(1.1); m_lower /= sqrt(1.1);
        while(m_upper - m_old - m_upper * m_upper / m_old * ratefact * CoolingRateFromU_adm(u, rho * m_upper / m_old, ne_guess, target) * dt < 0)
        {
            m_upper *= 1.1; m_lower *= 1.1;
        }
    }
    if(m - m_old - m_old * ratefact * LambdaNet * dt > 0)
    {
        m_lower /= sqrt(1.1); m_upper *= sqrt(1.1);
        while(m_lower - m_old - m_lower * m_lower / m_old * ratefact * CoolingRateFromU_adm(u, rho * m_lower / m_old, ne_guess, target) * dt > 0)
        {
            m_upper /= 1.1; m_lower /= 1.1;
        }
    }

    do
    {
        m = 0.5 * (m_lower + m_upper);
        LambdaNet = CoolingRateFromU_adm(u, rho * m / m_old, ne_guess, target);
        if(m - m_old - m * m / m_old * ratefact * LambdaNet * dt > 0) {m_upper = m;} else {m_lower = m;}
        dm = m_upper - m_lower;
        iter++;
        if(iter >= (MAXITER - 10)) {printf("->m= %g\n", m);}
    }
    while(fabs(dm / m) > 1.0e-6 && iter < MAXITER);
    if(iter >= MAXITER) {printf("failed to converge in DoInstabilityCooling(): m_in=%g u_in=%g rho=%g dt=%g fac=%g ne_in=%g target=%d \n",m_old,u,rho,dt,fac,ne_guess,target); endrun(11);}
    return m;
}








/* this function determines the electron fraction, and hence the mean molecular weight. With it arrives at a self-consistent temperature.
 * Ionization abundances and the rates for the emission are also computed */
double convert_u_to_temp_adm(double u, double rho, int target, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess, double *mu_guess)
{
    int iter = 0;
    double temp, temp_old, temp_old_old = 0, temp_new, prefac_fun_old, prefac_fun, fac, err_old, err_new, T_bracket_errneg = 0, T_bracket_errpos = 0, T_bracket_min = 0, T_bracket_max = 1.e20, bracket_sign = 0; // double max = 0;
    double u_input = u, rho_input = rho, temp_guess;
    double T_0 = u * PROTONMASS / BOLTZMANN; // this is the dimensional temperature, which since u is fixed is -frozen- in this calculation: we can work dimensionlessly below
    temp_guess = (GAMMA(target)-1) * T_0; // begin assuming mu ~ 1
    *mu_guess = Get_Gas_Mean_Molecular_Weight_mu(temp_guess, rho, nH0_guess, ne_guess, 0., target); // get mu with that temp
    prefac_fun = (GAMMA(target)-1) * (*mu_guess); // dimensionless pre-factor determining the temperature
    err_new = prefac_fun - temp_guess / T_0; // define initial error from this iteration
    if(err_new < 0) {T_bracket_errneg = temp_guess;} else {T_bracket_errpos = temp_guess;}
    temp = prefac_fun * T_0; // re-calculate temo with the new mu

    do
    {
        //qfun_old = *ne_guess; // guess for ne
        //qfun_old = *mu_guess; // guess for mu
        prefac_fun_old = prefac_fun;
        err_old = err_new; // error from previous timestep
        find_abundances_and_rates_adm(log10(temp), rho, target, -1, 0, ne_guess, nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu_guess); // all the thermo variables for this T
        prefac_fun = (GAMMA(target)-1) * (*mu_guess); // new value of the dimensionless pre-factor we need to solve
        temp_old = temp; // guess for T we just used
        temp_new = prefac_fun * T_0; // updated temp using the new values from the iteration of find_abundances_and_rates above
        err_new = (temp_new - temp_old) / T_0; // new error
        if(T_bracket_errpos == 0) {if(err_new > 0) {T_bracket_errpos = temp_old;} else {T_bracket_errneg = temp_old;}} // update the bracket values to the new T while its error still reflects here
        if(T_bracket_errneg == 0) {if(err_new < 0) {T_bracket_errneg = temp_old;} else {T_bracket_errpos = temp_old;}} // update the bracket values to the new T while its error still reflects here
        if(T_bracket_errneg > 0 && T_bracket_errpos > 0)
        {
            if(bracket_sign == 0) {if(T_bracket_errpos > T_bracket_errneg) {bracket_sign=1;} else {bracket_sign=-1;}}
            if(err_new > 0) {
                if(bracket_sign > 0) {T_bracket_errpos = DMIN(T_bracket_errpos, temp_old); /* Tpos>Tneg */} else {T_bracket_errpos = DMAX(T_bracket_errpos, temp_old); /* Tpos<Tneg */}
            } else {
                if(bracket_sign > 0) {T_bracket_errneg = DMAX(T_bracket_errneg, temp_old); /* Tpos>Tneg */} else {T_bracket_errneg = DMIN(T_bracket_errneg, temp_old); /* Tpos<Tneg */}
            } /* update bracket values if we can */
            if(bracket_sign > 0) {T_bracket_max=T_bracket_errpos; T_bracket_min=T_bracket_errneg;} else {T_bracket_max=T_bracket_errneg; T_bracket_min=T_bracket_errpos;}
        }

        //max = DMAX(max, temp_new * (*mu_guess) * HYDROGEN_MASSFRAC * fabs((*ne_guess - qfun_old) / (temp_new - temp_old + 1.0))); // old iteration: hardwired assumption that ne is only varying quanity in mu, and that Tmin ~ 1e4 or so
        //max = DMAX(max , temp_new / (*mu_guess) * fabs(*mu_guess - qfun_old) / (fabs(temp_new - temp_old) + 1.e-4*(All.MinGasTemp+0.1))); // newer - more flexible mu, and dimensionless T dependence
        //temp = temp_old + (temp_new - temp_old) / (1 + max);

        if(fabs(prefac_fun-prefac_fun_old) < 1.e-4 || fabs(temp_new-temp_old)/(temp_new+temp_old) < 1.e-4) {break;} // break pre-emptively if we'll trigger a nan below
        fac = (prefac_fun-prefac_fun_old)*T_0 / (temp_old-temp_old_old); // numerical derivative factor: want to use this to limit for convergence

        if(fac > 1) {fac = 1;} // don't allow us to move in the opposite direction from the new evaluation (should 'guess' in the direction of T_new-T_old) -- this tells us Newton-Raphson/Secant-type method fails here, so we simply follow the iteration to t_new
        if(fac > 0.9) {fac=0.9;} // don't allow us to 'jump' by a factor >10 times the temperature difference (arbitrary choice, slows convergence a bit but helps limit bad overshoot)
        if(fac < -9999.) {fac=-9999.;} // don't allow smaller step than 1e-4 times the temperature difference (since that's below our error tolerance anyways)

        temp = temp_old + (temp_new - temp_old) * 1./(1. - fac); // standard Newton-Raphson-type (technically Secant method) iteration
        if(temp < 0.5*temp_old) {temp = 0.5*temp_old;} // limiter to prevent un-physical overshoot before we have bracketing established
        if(temp > 3.0*temp_old) {temp = 3.0*temp_old;} // limiter to prevent un-physical overshoot before we have bracketing established

        temp = temp_old + (temp_new - temp_old) * 1./(1. - fac); // standard Newton-Raphson iteration

        if(T_bracket_errneg > 0 && T_bracket_errpos > 0) // if have bracketing and this wants to go outside brackets, revert to bisection
        {
            if(temp >= T_bracket_max || temp <= T_bracket_min) {temp = sqrt(T_bracket_min*T_bracket_max);} // bisect (in log-space)
        }
#ifndef RT_INFRARED
        if(fabs(temp-temp_old_old)/(temp+temp_old_old) < 1.e-3) {double wt=get_random_number(12*iter+340*ThisTask+5435*target); temp=(wt*temp_old + (1.-wt)*temp_new);}
#endif
        temp_old_old = temp_old;
        iter++;
        if(iter > (MAXITER - 10)) {printf("-> temp_next/new/old/oldold=%g/%g/%g/%g ne=%g mu=%g rho=%g iter=%d target=%d err_new/prev=%g/%g gamma_minus_1_mu_new/prev=%g/%g Brackets: Error_bracket_positive=%g Error_bracket_negative=%g T_bracket_Min/Max=%g/%g fac_for_SecantDT=%g \n", temp,temp_new,temp_old,temp_old_old,*ne_guess, (*mu_guess) ,rho,iter,target,err_new,err_old,prefac_fun,prefac_fun_old,T_bracket_errpos,T_bracket_errneg,T_bracket_min,T_bracket_max,fac); fflush(stdout);}
    }
    while(
#ifdef RT_INFRARED
        (fabs(temp - temp_old) > 1e-3 * temp) && iter < MAXITER);
#else
          ((fabs(temp - temp_old) > 0.25 * temp) ||
           ((fabs(temp - temp_old) > 0.1 * temp) && (temp > 20.)) ||
           ((fabs(temp - temp_old) > 0.05 * temp) && (temp > 200.)) ||
           ((fabs(temp - temp_old) > 0.01 * temp) && (temp > 200.) && (iter<100)) ||
           ((fabs(temp - temp_old) > 1.0e-3 * temp) && (temp > 200.) && (iter<10))) && iter < MAXITER);
#endif
    if(iter >= MAXITER) {printf("failed to converge in convert_u_to_temp(): u_input= %g rho_input=%g n_elec_input=%g target=%d\n", u_input, rho_input, *ne_guess, target); endrun(12);}

    if(temp<=0) temp=pow(10.0,Tmin_adm);
    if(log10(temp)<Tmin_adm) temp=pow(10.0,Tmin_adm);
    return temp;
}




/* this function computes the actual ionization states, relative abundances, and returns the ionization/recombination rates if needed */
double find_abundances_and_rates_adm(double logT, double rho, int target, double shieldfac, int return_cooling_mode,
                                 double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess,
                                 double *mu_guess)
{
    int j, niter;
    double Tlow, Thi, flow, fhi, t, gJH0ne, gJHe0ne, gJHepne, logT_input, rho_input, ne_input, neold, nenew;
    double bH0, bHep, bff, aHp, aHep, aHepp, ad, geH0, geHe0, geHep, EPSILON_SMALL=1.e-40;
    double n_elec, nH0, nHe0, nHp, nHep, nHepp; /* ionization states */
    logT_input = logT; rho_input = rho; ne_input = *ne_guess; /* save inputs (in case of failed convergence below) */
    if(!isfinite(logT)) {logT=Tmin_adm;}    /* nan trap (just in case) */
    if(!isfinite(rho)) {logT=Tmin_adm;}

    if(logT <= Tmin_adm)		/* everything neutral */
    {
        nH0 = 1.0; nHe0 = yhelium(target); nHp = 0; nHep = 0; nHepp = 0; n_elec = 0;
        *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec;
        *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, 0, target);
        return 0;
    }
    if(logT >= Tmax_adm)		/* everything is ionized */
    {
        nH0 = 0; nHe0 = 0; nHp = 1.0; nHep = 0; nHepp = yhelium(target); n_elec = nHp + 2.0 * nHepp;
        *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec;
        *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, 1.e3, target);
        return 0;
    }

    /* initialize quantities needed for iteration below */
    t = (logT - Tmin_adm) / deltaT_adm;
    j = (int) t;
    if(j<0){j=0;}
    if(j>NCOOLTAB_ADM){
        PRINT_WARNING("j>NCOOLTAB_ADM : j=%d t %g Tlow %g Thi %g logT %g Tmin_adm %g deltaT_adm %g \n",j,t,Tmin_adm+deltaT_adm*j,Tmin_adm+deltaT_adm*(j+1),logT,Tmin_adm,deltaT_adm);fflush(stdout);
        j=NCOOLTAB_ADM;
    }
    Tlow = Tmin_adm + deltaT_adm * j;
    Thi = Tlow + deltaT_adm;
    fhi = t - j;
    flow = 1 - fhi;
    if(*ne_guess == 0) /* no guess provided, try to start from something sensible */
    {
        *ne_guess = 1.0;
        if(logT < 3.8) {*ne_guess = 0.1;}
        if(logT < 2) {*ne_guess = 1.e-10;}
    }
    /* CAFG: this is the density that we should use for UV background threshold */
    double local_gammamultiplier = return_local_gammamultiplier_adm(target); // account for local UVB terms in some expressions below
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    if(shieldfac < 0) {shieldfac = return_uvb_shieldfac_adm(target, local_gammamultiplier*gJH0_adm/1.0e-12, nHcgs, logT);} // if < 0, that's a key to tell us this needs to be recalculated
    n_elec = *ne_guess; if(!isfinite(n_elec)) {n_elec=1;}
    neold = n_elec; niter = 0;
    double dt = 0, fac_noneq_cgs = 0, necgs = n_elec * nHcgs; /* more initialized quantities */
    if(target >= 0) {dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target);} // dtime [code units]
    fac_noneq_cgs = (dt * UNIT_TIME_IN_CGS) * necgs; // factor needed below to asses whether timestep is larger/smaller than recombination time

#if defined(RT_CHEM_PHOTOION)
    double c_light_ne=0, Sigma_particle=0, abs_per_kappa_dt=0;
    if(target >= 0)
    {
        double L_particle = Get_Particle_Size(target)*All.cf_atime; // particle effective size/slab thickness
        double cx_to_kappa = HYDROGEN_MASSFRAC / PROTONMASS * UNIT_MASS_IN_CGS; // pre-factor for converting cross sections into opacities
        Sigma_particle = cx_to_kappa * P[target].Mass / (M_PI*L_particle*L_particle); // effective surface density through particle
        abs_per_kappa_dt = cx_to_kappa * C_LIGHT_CODE_REDUCED * (SphP[target].Density*All.cf_a3inv) * dt; // fractional absorption over timestep
        nH0 = SphP[target].HI; // need to initialize a value for the iteration below
#ifdef RT_CHEM_PHOTOION_HE
        nHe0 = SphP[target].HeI; nHep = SphP[target].HeII; // need to intialize a value for the iteration below
#endif
    }
#endif

    /* evaluate number densities iteratively (cf KWH eqns 33-38) in units of nH */
    do
    {
        niter++;

        aHp = flow * AlphaHp_adm[j] + fhi * AlphaHp_adm[j + 1];
        aHep = flow * AlphaHep_adm[j] + fhi * AlphaHep_adm[j + 1];
        aHepp = flow * AlphaHepp_adm[j] + fhi * AlphaHepp_adm[j + 1];
        ad = flow * Alphad_adm[j] + fhi * Alphad_adm[j + 1];
        geH0 = flow * GammaeH0_adm[j] + fhi * GammaeH0_adm[j + 1];
        geH0 = DMAX(geH0, EPSILON_SMALL);
        geHe0 = flow * GammaeHe0_adm[j] + fhi * GammaeHe0_adm[j + 1];
        geHe0 = DMAX(geHe0, EPSILON_SMALL);
        geHep = flow * GammaeHep_adm[j] + fhi * GammaeHep_adm[j + 1];
        geHep = DMAX(geHep, EPSILON_SMALL);
        fac_noneq_cgs = (dt * UNIT_TIME_IN_CGS) * necgs; // factor needed below to asses whether timestep is larger/smaller than recombination time
        if(necgs <= 1.e-25 || J_UV_ADM == 0)
        {
            gJH0ne = gJHe0ne = gJHepne = 0;
        }
        else
        {
            /* account for self-shielding in calculating UV background effects */
            gJH0ne = gJH0_adm * local_gammamultiplier / necgs * shieldfac; // check units, should be = c_light * n_photons_vol * rt_ion_sigma_HI[0] / necgs;
            gJH0ne = DMAX(gJH0ne, EPSILON_SMALL); if(!isfinite(gJH0ne)) {gJH0ne=0;} // need traps here b/c very small numbers assigned in some newer TREECOOL versions cause a nan underflow
            gJHe0ne = gJHe0_adm * local_gammamultiplier / necgs * shieldfac;
            gJHe0ne = DMAX(gJHe0ne, EPSILON_SMALL); if(!isfinite(gJHe0ne)) {gJHe0ne=0;}
            gJHepne = gJHep_adm * local_gammamultiplier / necgs * shieldfac;
            gJHepne = DMAX(gJHepne, EPSILON_SMALL); if(!isfinite(gJHepne)) {gJHepne=0;}
        }
#if defined(RT_DISABLE_UV_BACKGROUND)
        gJH0ne = gJHe0ne = gJHepne = 0;
#endif
#if defined(RT_CHEM_PHOTOION)
        /* add in photons from explicit radiative transfer (on top of assumed background) */
        if(target >= 0)
        {
            int k;
            c_light_ne = C_LIGHT / ((MIN_REAL_NUMBER + necgs) * UNIT_LENGTH_IN_CGS); // want physical cgs units for quantities below
            double gJH0ne_0=gJH0_adm * local_gammamultiplier / (MIN_REAL_NUMBER + necgs), gJHe0ne_0=gJHe0_adm * local_gammamultiplier / (MIN_REAL_NUMBER + necgs), gJHepne_0=gJHep_adm * local_gammamultiplier / (MIN_REAL_NUMBER + necgs); // need a baseline, so we don't over-shoot below
            gJH0ne = DMAX(gJH0ne, EPSILON_SMALL); if(!isfinite(gJH0ne)) {gJH0ne=0;} // need traps here b/c very small numbers assigned in some newer TREECOOL versions cause a nan underflow
            gJHe0ne = DMAX(gJHe0ne, EPSILON_SMALL); if(!isfinite(gJHe0ne)) {gJHe0ne=0;}
            gJHepne = DMAX(gJHepne, EPSILON_SMALL); if(!isfinite(gJHepne)) {gJHepne=0;}
#if defined(RT_DISABLE_UV_BACKGROUND)
            gJH0ne_0=gJHe0ne_0=gJHepne_0=MAX_REAL_NUMBER;
#endif
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if((k==RT_FREQ_BIN_H0)||(k==RT_FREQ_BIN_He0)||(k==RT_FREQ_BIN_He1)||(k==RT_FREQ_BIN_He2))
                {
                    double c_ne_time_n_photons_vol = c_light_ne * rt_return_photon_number_density(target,k); // gives photon flux
                    double cross_section_ion, dummy, thold=1.0e20;
#ifdef GALSF
                    if(All.ComovingIntegrationOn) {thold=1.0e10;}
#endif
                    if(rt_ion_G_HI[k] > 0)
                    {
                        cross_section_ion = nH0 * rt_ion_sigma_HI[k];
                        dummy = rt_ion_sigma_HI[k] * c_ne_time_n_photons_vol;// egy per photon x cross section x photon flux (w attenuation factors already included in flux/energy update:) * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt);
                        if(dummy > thold*gJH0ne_0) {dummy = thold*gJH0ne_0;}
                        gJH0ne += dummy;
                    }
#ifdef RT_CHEM_PHOTOION_HE
                    if(rt_ion_G_HeI[k] > 0)
                    {
                        cross_section_ion = nHe0 * rt_ion_sigma_HeI[k];
                        dummy = rt_ion_sigma_HeI[k] * c_ne_time_n_photons_vol;// * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt);
                        if(dummy > thold*gJHe0ne_0) {dummy = thold*gJHe0ne_0;}
                        gJHe0ne += dummy;
                    }
                    if(rt_ion_G_HeII[k] > 0)
                    {
                        cross_section_ion = nHep * rt_ion_sigma_HeII[k];
                        dummy = rt_ion_sigma_HeII[k] * c_ne_time_n_photons_vol;// * slab_averaging_function(cross_section_ion * Sigma_particle); // * slab_averaging_function(cross_section_ion * abs_per_kappa_dt);
                        if(dummy > thold*gJHepne_0) {dummy = thold*gJHepne_0;}
                        gJHepne += dummy;
                    }
#endif
                }
            }
        }
#endif


        nH0 = aHp / (MIN_REAL_NUMBER + aHp + geH0 + gJH0ne);	/* eqn (33) */
#ifdef RT_CHEM_PHOTOION
        if(target >= 0) {nH0 = (SphP[target].HI + fac_noneq_cgs * aHp) / (1 + fac_noneq_cgs * (aHp + geH0 + gJH0ne));} // slightly more general formulation that gives linear update but interpolates to equilibrium solution when dt >> dt_recombination
#endif
        nHp = 1.0 - nH0;		/* eqn (34) */

        if( ((gJHe0ne + geHe0) <= MIN_REAL_NUMBER) || (aHepp <= MIN_REAL_NUMBER) ) 	/* no ionization at all */
        {
            nHep = 0.0;
            nHepp = 0.0;
            nHe0 = yhelium(target);
        }
        else
        {
            nHep = yhelium(target) / (1.0 + (aHep + ad) / (geHe0 + gJHe0ne) + (geHep + gJHepne) / aHepp);	/* eqn (35) */
            nHe0 = nHep * (aHep + ad) / (geHe0 + gJHe0ne);	/* eqn (36) */
            nHepp = nHep * (geHep + gJHepne) / aHepp;	/* eqn (37) */
        }
#if defined(RT_CHEM_PHOTOION) && defined(RT_CHEM_PHOTOION_HE)
        if(target >= 0)
        {
            double yHe = yhelium(target); // will use helium fraction below
            nHep = SphP[target].HeII + yHe * fac_noneq_cgs * (geHe0 + gJHe0ne) - SphP[target].HeIII * (fac_noneq_cgs*(geHe0 + gJHe0ne - aHepp) / (1.0 + fac_noneq_cgs*aHepp));
            nHep /= 1.0 + fac_noneq_cgs*(geHe0 + gJHe0ne + aHep + ad + geHep + gJHepne) + (fac_noneq_cgs*(geHe0 + gJHe0ne - aHepp) / (1.0 + fac_noneq_cgs*aHepp)) * fac_noneq_cgs*(geHep + gJHepne);
            if(nHep < 0) {nHep=0;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            if(nHep > yHe) {nHep=yHe;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            nHepp = (SphP[target].HeIII + SphP[target].HeII * fac_noneq_cgs*(geHep + gJHepne)) / (1. + fac_noneq_cgs*aHepp);
            if(nHepp < 0) {nHepp=0;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            if(nHepp > yHe-nHep) {nHepp=yHe-nHep;} // check if this exceeded valid limits (can happen in 'overshoot' during iteration)
            nHe0 = yHe - (nHep + nHepp); // remainder is neutral
        }
#endif
        if(!isfinite(n_elec)) {printf("target=%d niter=%d logT=%g n_elec/old=%g/%g nHp/nHep/nHepp=%g/%g/%g nHcgs=%g yHe=%g dt=%g shieldfac/local_gammamult=%g/%g aHp/aHep/aHepp=%g/%g/%g geH0/geHe0/geHep=%g/%g/%g gJH0ne/gJHe0ne/gJHepne=%g/%g/%g \n",target,niter,logT,n_elec,neold,nHp,nHep,nHepp,nHcgs,yhelium(target),dt,shieldfac,local_gammamultiplier,aHp,aHep,aHepp,geH0,geHe0,geHep,gJH0ne,gJHe0ne,gJHepne);}

        neold = n_elec;
        n_elec = nHp + nHep + 2 * nHepp;	/* eqn (38) */
        necgs = n_elec * nHcgs;

        if(J_UV_ADM == 0) break;

        nenew = 0.5 * (n_elec + neold);
        n_elec = nenew;
        if(!isfinite(n_elec)) {n_elec=1;}
        necgs = n_elec * nHcgs;

        double dneTHhold = DMAX(n_elec*0.01 , 1.0e-4);
        if(fabs(n_elec - neold) < dneTHhold) break;

        if(niter > (MAXITER - 10)) {printf("n_elec= %g/%g/%g yh=%g nHcgs=%g niter=%d\n", n_elec,neold,nenew, yhelium(target), nHcgs, niter);}
    }
    while(niter < MAXITER);

    if(niter >= MAXITER) {printf("failed to converge in find_abundances_and_rates(): logT_input=%g  rho_input=%g  ne_input=%g target=%d shieldfac=%g cooling_return=%d", logT_input, rho_input, ne_input, target, shieldfac, return_cooling_mode); endrun(13);}

    bH0 = flow * BetaH0_adm[j] + fhi * BetaH0_adm[j + 1];
    bHep = flow * BetaHep_adm[j] + fhi * BetaHep_adm[j + 1];
    bff = flow * Betaff_adm[j] + fhi * Betaff_adm[j + 1];
    *nH0_guess=nH0; *nHe0_guess=nHe0; *nHp_guess=nHp; *nHep_guess=nHep; *nHepp_guess=nHepp; *ne_guess=n_elec; /* write to send back */
    *mu_guess=Get_Gas_Mean_Molecular_Weight_mu(pow(10.,logT), rho, nH0_guess, ne_guess, sqrt(shieldfac)*(gJH0_adm/2.29e-10), target);
    if(target >= 0) /* if this is a cell, update some of its thermodynamic stored quantities */
    {
        SphP[target].Ne = n_elec;
#if defined(OUTPUT_MOLECULAR_FRACTION)
        SphP[target].MolecularMassFraction = Get_Gas_Molecular_Mass_Fraction(target, pow(10.,logT), nH0, n_elec, sqrt(shieldfac)*(gJH0_adm/2.29e-10));
#endif
    }

    /* now check if we want to return the ionization/recombination heating/cooling rates calculated with all the above quantities */
    if(return_cooling_mode==1)
    {
        /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */
        double LambdaExcH0 = bH0 * n_elec * nH0;
        double LambdaExcHep = bHep * n_elec * nHep;
        double LambdaExc = LambdaExcH0 + LambdaExcHep;	/* collisional excitation */

        double LambdaIonH0 = 2.18e-11 * geH0 * n_elec * nH0;
        double LambdaIonHe0 = 3.94e-11 * geHe0 * n_elec * nHe0;
        double LambdaIonHep = 8.72e-11 * geHep * n_elec * nHep;
        double LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* collisional ionization */

        double T_lin = pow(10.0, logT);
        double LambdaRecHp = 1.036e-16 * T_lin * n_elec * (aHp * nHp);
        double LambdaRecHep = 1.036e-16 * T_lin * n_elec * (aHep * nHep);
        double LambdaRecHepp = 1.036e-16 * T_lin * n_elec * (aHepp * nHepp);
        double LambdaRecHepd = 6.526e-11 * ad * n_elec * nHep;
        double LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd; /* recombination */

        double LambdaFF = bff * (nHp + nHep + 4 * nHepp) * n_elec; /* free-free (Bremsstrahlung) */

        double Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF; /* sum all of the above */
        return Lambda; /* send it back */
    }
    return 0;
} // end of find_abundances_and_rates() //



/*  this function first computes the self-consistent temperature and abundance ratios, and then it calculates (heating rate-cooling rate)/n_h^2 in cgs units */
double CoolingRateFromU_adm(double u, double rho, double ne_guess, int target)
{
    double nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu; nH0_guess = DMAX(0,DMIN(1,1.-ne_guess/1.2));
    double temp = convert_u_to_temp_adm(u, rho, target, &ne_guess, &nH0_guess, &nHp_guess, &nHe0_guess, &nHep_guess, &nHepp_guess, &mu);
    return CoolingRate_adm(log10(temp), rho, ne_guess, target);
}





extern FILE *fd;



/*  Calculates (heating rate-cooling rate)/n_h^2 in cgs units
 */
double CoolingRate_adm(double logT, double rho, double n_elec_guess, int target)
{
    double n_elec=n_elec_guess, nH0, nHe0, nHp, nHep, nHepp, mu; /* ionization states [computed below] */
    double Lambda, Heat, LambdaFF, LambdaCompton, LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
    double LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd, redshift, T, shieldfac, LambdaMol, LambdaMetal;
    double nHcgs = HYDROGEN_MASSFRAC * rho / PROTONMASS;	/* hydrogen number dens in cgs units */
    LambdaMol=0; LambdaMetal=0; LambdaCompton=0;
    if(logT <= Tmin_adm) {logT = Tmin_adm + 0.5 * deltaT_adm;}	/* floor at Tmin_adm */
    if(!isfinite(rho)) {return 0;}
    T = pow(10.0, logT);

    /* some blocks below to define useful variables before calculation of cooling rates: */

#ifdef COOL_METAL_LINES_BY_SPECIES
    double *Z;
    if(target>=0)
    {
        Z = P[target].Metallicity;
    } else { /* initialize dummy values here so the function doesn't crash, if called when there isn't a target particle */
        int k; double Zsol[NUM_METAL_SPECIES]; for(k=0;k<NUM_METAL_SPECIES;k++) {Zsol[k]=All.SolarAbundances[k];}
        Z = Zsol;
    }
#endif
    double local_gammamultiplier = return_local_gammamultiplier_adm(target);
    shieldfac = return_uvb_shieldfac_adm(target, local_gammamultiplier*gJH0_adm/1.0e-12, nHcgs, logT);

#if defined(COOL_LOW_TEMPERATURES)
    double Tdust = 30., LambdaDust = 0.; /* set variables needed for dust heating/cooling. if dust cooling not calculated, default to 0 */
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) && (FLAG_NOT_IN_PUBLIC_CODE > 2)) || defined(SINGLE_STAR_SINK_DYNAMICS)
    Tdust = get_equilibrium_dust_temperature_estimate_adm(target, shieldfac);
#endif
#endif


#if defined(RT_CHEM_PHOTOION) || defined(RT_PHOTOELECTRIC)
    double Sigma_particle = 0, abs_per_kappa_dt = 0, cx_to_kappa = 0;
    if(target >= 0)
    {
        double L_particle = Get_Particle_Size(target)*All.cf_atime; // particle effective size/slab thickness
        double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(target); // dtime [code units]
        Sigma_particle = P[target].Mass / (M_PI*L_particle*L_particle); // effective surface density through particle
        abs_per_kappa_dt = C_LIGHT_CODE_REDUCED * (SphP[target].Density*All.cf_a3inv) * dt; // fractional absorption over timestep
        cx_to_kappa = HYDROGEN_MASSFRAC / PROTONMASS * UNIT_MASS_IN_CGS; // pre-factor for converting cross sections into opacities
    }
#endif
    if(logT < Tmax_adm)
    {
        /* get ionization states for H and He with associated ionization, collision, recombination, and free-free heating/cooling */
        Lambda = find_abundances_and_rates_adm(logT, rho, target, shieldfac, 1, &n_elec, &nH0, &nHp, &nHe0, &nHep, &nHepp, &mu);

        LambdaCompton = evaluate_Compton_heating_cooling_rate_adm(target,T,nHcgs,n_elec,shieldfac); /* note this can have either sign: heating or cooling */
        if(LambdaCompton > 0) {Lambda += LambdaCompton;}

#ifdef COOL_METAL_LINES_BY_SPECIES
        /* can restrict to low-densities where not self-shielded, but let shieldfac (in ne) take care of this self-consistently */
        if((J_UV_ADM != 0)&&(logT > 4.00))
        {
            /* cooling rates tabulated for each species from Wiersma, Schaye, & Smith tables (2008) */
            LambdaMetal = GetCoolingRateWSpecies_adm(nHcgs, logT, Z); //* nHcgs*nHcgs;
            /* tables normalized so ne*ni/(nH*nH) included already, so just multiply by nH^2 */
            /* (sorry, -- dont -- multiply by nH^2 here b/c that's how everything is normalized in this function) */
            LambdaMetal *= n_elec;
            /* (modified now to correct out tabulated ne so that calculated ne can be inserted; ni not used b/c it should vary species-to-species */
            Lambda += LambdaMetal;
#if defined(OUTPUT_COOLRATE_DETAIL)
            if(target >= 0) {SphP[target].MetalCoolingRate = LambdaMetal;}
#endif
        }
#endif

#ifdef COOL_LOW_TEMPERATURES
        if(logT <= 5.3)
        {
            /* approx to cooling function for solar metallicity and nH=1 cm^(-3) -- want to do something
             much better, definitely, but for now use this just to get some idea of system with cooling to very low-temp */
            LambdaMol = 2.8958629e-26/(pow(T/125.21547,-4.9201887)+pow(T/1349.8649,-1.7287826)+pow(T/6450.0636,-0.30749082));
            LambdaMol *= (1-shieldfac) / (1. + nHcgs/700.); // above the critical density, cooling rate suppressed by ~1/n; use critical density of CO[J(1-0)] as a proxy for this
            double Z_sol=1, truncation_factor=1; /* if don't have actual metallicities, we'll assume solar */
            if(logT>4.5) {double dx=(logT-4.5)/0.20; truncation_factor *= exp(-DMIN(dx*dx,40.));} /* continuous cutoff here just to avoid introducing artificial features in temperature-density */
#ifdef COOL_METAL_LINES_BY_SPECIES
            Z_sol = Z[0] / All.SolarAbundances[0]; /* use actual metallicity for this */
#endif
            LambdaMol *= (1+Z_sol)*(0.001 + 0.1*nHcgs/(1.+nHcgs) + 0.09*nHcgs/(1.+0.1*nHcgs) + Z_sol*Z_sol/(1.0+nHcgs)); // gives very crude estimate of metal-dependent terms //
            LambdaMol *= truncation_factor; // cutoff factor from above for where the tabulated rates take over at high temperatures
            if(!isfinite(LambdaMol)) {LambdaMol=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            Lambda += LambdaMol;

            /* now add the dust cooling/heating terms */
            LambdaDust = 1.116e-32 * (Tdust-T) * sqrt(T)*(1.-0.8*exp(-75./T)) * Z_sol;  // Meijerink & Spaans 2005; Hollenbach & McKee 1979,1989 //
#ifdef RT_INFRARED
            if(target >= 0) {LambdaDust = get_rt_ir_lambdadust_effective(T, rho, &nH0, &n_elec, target);} // call our specialized subroutine, because radiation and gas energy fields are co-evolving and tightly-coupled here //
#endif
            if(T>3.e5) {double dx=(T-3.e5)/2.e5; LambdaDust *= exp(-DMIN(dx*dx,40.));} /* needs to truncate at high temperatures b/c of dust destruction */
            LambdaDust *= truncation_factor; // cutoff factor from above for where the tabulated rates take over at high temperatures
#ifdef RT_INFRARED
            SphP[target].LambdaDust = LambdaDust;
#endif
            if(!isfinite(LambdaDust)) {LambdaDust=0;} // here to check vs underflow errors since dividing by some very small numbers, but in that limit Lambda should be negligible
            if(LambdaDust<0) {Lambda -= LambdaDust;} /* add the -positive- Lambda-dust associated with cooling */
        }
#endif



        Heat = 0;  /* Now, collect heating terms */

        if(J_UV_ADM != 0) {Heat += local_gammamultiplier * (nH0 * epsH0_adm + nHe0 * epsHe0_adm + nHep * epsHep_adm) / nHcgs * shieldfac;} // shieldfac allows for self-shielding from background
#if defined(RT_DISABLE_UV_BACKGROUND)
        Heat = 0;
#endif
#if defined(RT_CHEM_PHOTOION)
        /* add in photons from explicit radiative transfer (on top of assumed background) */
        if((target >= 0) && (nHcgs > MIN_REAL_NUMBER))
        {
            int k; double c_light_nH = C_LIGHT / (nHcgs * UNIT_LENGTH_IN_CGS) * UNIT_ENERGY_IN_CGS; // want physical cgs units for quantities below
            for(k = 0; k < N_RT_FREQ_BINS; k++)
            {
                if((k==RT_FREQ_BIN_H0)||(k==RT_FREQ_BIN_He0)||(k==RT_FREQ_BIN_He1)||(k==RT_FREQ_BIN_He2))
                {
                    double c_nH_time_n_photons_vol = c_light_nH * rt_return_photon_number_density(target,k); // gives photon flux
                    double cross_section_ion, kappa_ion, dummy;
                    if(rt_ion_G_HI[k] > 0)
                    {
                        cross_section_ion = nH0 * rt_ion_sigma_HI[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HI[k] * cross_section_ion * c_nH_time_n_photons_vol;// (egy per photon x cross section x photon flux) :: attenuation factors [already in flux/energy update]: * slab_averaging_function(kappa_ion * Sigma_particle); // egy per photon x cross section x photon flux (w attenuation factors) // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);
                        Heat += dummy;
                    }
                    if(rt_ion_G_HeI[k] > 0)
                    {
                        cross_section_ion = nHe0 * rt_ion_sigma_HeI[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HeI[k] * cross_section_ion * c_nH_time_n_photons_vol;// * slab_averaging_function(kappa_ion * Sigma_particle); // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);
                        Heat += dummy;
                    }
                    if(rt_ion_G_HeII[k] > 0)
                    {
                        cross_section_ion = nHep * rt_ion_sigma_HeII[k];
                        kappa_ion = cx_to_kappa * cross_section_ion;
                        dummy = rt_ion_G_HeII[k] * cross_section_ion * c_nH_time_n_photons_vol;// * slab_averaging_function(kappa_ion*Sigma_particle); // * slab_averaging_function(kappa_ion * abs_per_kappa_dt);
                        Heat += dummy;
                    }
                }
            }
        }
#endif


#ifdef COOL_LOW_TEMPERATURES
        /* if COSMIC_RAYS is not enabled, but low-temperature cooling is on, we account for the CRs as a heating source using
         a more approximate expression (assuming the mean background of the Milky Way clouds) */
        if(logT <= 5.2)
        {
            double prefac_CR=1.; if(All.ComovingIntegrationOn) {
                double rhofac = rho / (1000.*COSMIC_BARYON_DENSITY_CGS);
                if(rhofac < 0.2) {prefac_CR=0;} else {if(rhofac > 200.) {prefac_CR=1;} else {prefac_CR=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off CR heating for any gas with density unless it's >1000 times the cosmic mean density
            double cr_zeta=1.e-16, e_per_cr_ioniz=8.8e-12; // either high background (zeta=1e-16), with softer spectrum (5.5eV per ionization), following e.g. van Dishoeck & Black (1986); or equivalently lower rate with higher ~20eV per ionization per Goldsmith & Langer (1978); this is formally degenerate here. however both produce ~3-10x higher rates than more modern estimates (basically in both cases, assuming a CR energy density of ~2-5 eV/cm^3, instead of more modern ~0.5-2 eV/cm^3
            Heat += prefac_CR * cr_zeta * (1. + 1.68*n_elec*HYDROGEN_MASSFRAC) / (1.e-2 + nHcgs) * e_per_cr_ioniz; // final result
        }
#endif

#if defined(COOL_LOW_TEMPERATURES)
        if(LambdaDust>0) {Heat += LambdaDust;} /* Dust collisional heating (Tdust > Tgas) */
#endif

        if(LambdaCompton<0) {Heat -= LambdaCompton;} /* Compton heating rather than cooling */

#if defined(FLAG_NOT_IN_PUBLIC_CODE) || defined(RT_PHOTOELECTRIC)
        /* Photoelectric heating following Bakes & Thielens 1994 (also Wolfire 1995); now with 'update' from Wolfire 2005 for PAH [fudge factor 0.5 below] */
        if((target >= 0) && (T < 1.0e6))
        {
#ifdef RT_PHOTOELECTRIC
            double photoelec = SphP[target].Rad_E_gamma[RT_FREQ_BIN_PHOTOELECTRIC] * (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; // convert to Habing field //
            if(photoelec > 0) {if(photoelec > 1.e4) {photoelec = 1.e4;}}
#endif
            if(photoelec > 0)
            {
                double LambdaPElec = 1.3e-24 * photoelec / nHcgs * P[target].Metallicity[0]/All.SolarAbundances[0];
                double x_photoelec = photoelec * sqrt(T) / (0.5 * (1.0e-12+n_elec) * nHcgs);
                LambdaPElec *= 0.049/(1+pow(x_photoelec/1925.,0.73)) + 0.037*pow(T/1.0e4,0.7)/(1+x_photoelec/5000.);
                Heat += LambdaPElec;
            }
        }
#endif
    }
  else				/* here we're outside of tabulated rates, T>Tmax_adm K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present.  Assumes no heating. */
      Heat = LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep = LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;
      nHp = 1.0; nHep = 0; nHepp = yhelium(target); n_elec = nHp + 2.0 * nHepp; /* very hot: H and He both fully ionized */

      LambdaFF = 1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * n_elec; // free-free
      LambdaCompton = evaluate_Compton_heating_cooling_rate_adm(target,T,nHcgs,n_elec,shieldfac); // Compton

      Lambda = LambdaFF + DMAX(LambdaCompton,0);
    }

    double Q = Heat - Lambda;
#if defined(OUTPUT_COOLRATE_DETAIL)
    if (target>=0){SphP[target].CoolingRate = Lambda; SphP[target].HeatingRate = Heat;}
#endif

#if defined(COOL_LOW_TEMPERATURES) && !defined(COOL_LOWTEMP_THIN_ONLY)
    /* if we are in the optically thick limit, we need to modify the cooling/heating rates according to the appropriate limits;
        this flag does so by using a simple approximation. we consider the element as if it were a slab, with a column density
        calculated from the simulation properties and the Sobolev approximation. we then assume it develops an equilibrium internal
        temperature structure on a radiative diffusion timescale much faster than the dynamical time, and so the surface radiation
        from a photosphere can be simply related to the local density by the optical depth to infinity. the equations here follow
        Rafikov, 2007 (ApJ, 662, 642):
            denergy/dt/dArea = sigma*T^4 / fc(tau)
            fc(tau) = tau^eta + 1/tau (taking chi, phi~1; the second term describes the optically thin limit, which is calculated above
                more accurately anyways - that was just Kirchoff's Law; so we only need to worry about the first term)
            eta = 4*(gamma-1) / [gamma*(1+alpha+beta*(gamma-1)/gamma)], where gamma=real polytropic index, and alpha/beta follow
                an opacity law kappa=kappa_0 * P^alpha * T^beta. for almost all the regimes of interest, however, eta~1, which is also
                what is obtained for a convectively stable slab. so we will use this.
            now, this gives sigma*T^4/tau * Area_eff / nHcgs as the 'effective' cooling rate in our units of Heat or Lambda above.
                the nHcgs just puts it in the same volumetric terms. The Area_eff must be defined as ~m_particle/surface_density
                to have the same meaning for a slab as assumed in Rafikov (and to integrate correctly over all particles in the slab,
                if/when the slab is resolved). We estimate this in our usual fashion with the Sobolev-type column density
            tau = kappa * surface_density; we estimate kappa ~ 5 cm^2/g * (0.001+Z/Z_solar), as the frequency-integrated kappa for warm
                dust radiation (~150K), weighted by the dust-to-gas ratio (with a floor for molecular absorption). we could make this
                temperature-dependent, though, fairly easily - for this particular problem it won't make much difference
        This rate then acts as an upper limit to the net heating/cooling calculated above (restricts absolute value)
     */
    if( (nHcgs > 0.1) && (target >= 0) )  /* don't bother at very low densities, since youre not optically thick, and protect from target=-1 with GALSF_EFFECTIVE_EQS */
    {
        double surface_density = evaluate_NH_from_GradRho(SphP[target].Gradients.Density,PPP[target].Hsml,SphP[target].Density,PPP[target].NumNgb,1,target);
        surface_density *= 0.2 * UNIT_SURFDEN_IN_CGS; // converts to cgs; 0.2 is a tuning factor so that the Masunaga & Inutsuka 2000 solution is reproduced
        double effective_area = 2.3 * PROTONMASS / surface_density; // since cooling rate is ultimately per-particle, need a particle-weight here
        double kappa_eff; // effective kappa, accounting for metal abundance, temperature, and density //
        if(T < 1500.)
        {
            if(T < 150.) {kappa_eff=0.0027*T*sqrt(T);} else {kappa_eff=5.;}
            kappa_eff *= P[target].Metallicity[0]/All.SolarAbundances[0];
            if(kappa_eff < 0.1) {kappa_eff=0.1;}
        } else {
            /* this is an approximate result for high-temperature opacities, but provides a pretty good fit from 1.5e3 - 1.0e9 K */
            double k_electron = 0.2 * (1. + HYDROGEN_MASSFRAC); //0.167 * n_elec; /* Thompson scattering (non-relativistic) */
            double k_molecular = 0.1 * P[target].Metallicity[0]; /* molecular line opacities */
            double k_Hminus = 1.1e-25 * sqrt(P[target].Metallicity[0] * rho) * pow(T,7.7); /* negative H- ion opacity */
            double k_Kramers = 4.0e25 * (1.+HYDROGEN_MASSFRAC) * (P[target].Metallicity[0]+0.001) * rho / (T*T*T*sqrt(T)); /* free-free, bound-free, bound-bound transitions */
            double k_radiative = k_molecular + 1./(1./k_Hminus + 1./(k_electron+k_Kramers)); /* approximate interpolation between the above opacities */
            double k_conductive = 2.6e-7 * n_elec * T*T/(rho*rho); //*(1+pow(rho/1.e6,0.67) /* e- thermal conductivity can dominate at low-T, high-rho, here it as expressed as opacity */
            kappa_eff = 1./(1./k_radiative + 1./k_conductive); /* effective opacity including both heat carriers (this is exact) */
        }
        double tau_eff = kappa_eff * surface_density;
        double Lambda_Thick_BlackBody = 5.67e-5 * (T*T*T*T) * effective_area / ((1.+tau_eff) * nHcgs);
        if(Q > 0) {if(Q > Lambda_Thick_BlackBody) {Q=Lambda_Thick_BlackBody;}} else {if(Q < -Lambda_Thick_BlackBody) {Q=-Lambda_Thick_BlackBody;}}
    }
#endif

#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target>=0){SphP[target].NetHeatingRateQ = Q;}
#endif
#ifdef OUTPUT_MOLECULAR_FRACTION
    if(target>0) {SphP[target].MolecularMassFraction = Get_Gas_Molecular_Mass_Fraction(target, T, nH0, n_elec, sqrt(shieldfac)*(gJH0_adm/2.29e-10));}
#endif

#ifndef COOLING_OPERATOR_SPLIT
    /* add the hydro energy change directly: this represents an additional heating/cooling term, to be accounted for in the semi-implicit solution determined here. this is more accurate when tcool << tdynamical */
    if(target >= 0) {Q += SphP[target].DtInternalEnergy / nHcgs;}
#if defined(OUTPUT_COOLRATE_DETAIL)
    if(target >= 0) {SphP[target].HydroHeatingRate = SphP[target].DtInternalEnergy / nHcgs;}
#endif
#endif

  return Q;
} // ends CoolingRate





void InitCoolMemory_adm(void)
{
    BetaH0_adm = (double *) mymalloc("BetaH0_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    BetaHep_adm = (double *) mymalloc("BetaHep_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    AlphaHp_adm = (double *) mymalloc("AlphaHp_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    AlphaHep_adm = (double *) mymalloc("AlphaHep_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    Alphad_adm = (double *) mymalloc("Alphad_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    AlphaHepp_adm = (double *) mymalloc("AlphaHepp_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    GammaeH0_adm = (double *) mymalloc("GammaeH0_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    GammaeHe0_adm = (double *) mymalloc("GammaeHe0_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    GammaeHep_adm = (double *) mymalloc("GammaeHep_adm", (NCOOLTAB_ADM + 1) * sizeof(double));
    Betaff_adm = (double *) mymalloc("Betaff_adm", (NCOOLTAB_ADM + 1) * sizeof(double));

#ifdef COOL_METAL_LINES_BY_SPECIES
    long i_nH=41; long i_T=176; long kspecies=(long)NUM_LIVE_SPECIES_FOR_COOLTABLES;
    SpCoolTable0 = (float *) mymalloc("SpCoolTable0",(kspecies*i_nH*i_T)*sizeof(float));
    if(All.ComovingIntegrationOn) {SpCoolTable1 = (float *) mymalloc("SpCoolTable1",(kspecies*i_nH*i_T)*sizeof(float));}
#endif
}



void MakeCoolingTable_adm(void)
     /* Set up interpolation tables in T for cooling rates given in KWH, ApJS, 105, 19
        Hydrogen, Helium III recombination rates and collisional ionization cross-sections are updated */
{
    int i; double T,Tfact;
    if(All.MinGasTemp > 0.0) {Tmin_adm = log10(All.MinGasTemp);} else {Tmin_adm=-1.0;} // set minimum temperature in this table to some very low value if zero, where none of the cooling approximations above make sense
    deltaT_adm = (Tmax_adm - Tmin_adm) / NCOOLTAB_ADM;
    /* minimum internal energy for neutral gas */
    for(i = 0; i <= NCOOLTAB_ADM; i++)
    {
        BetaH0_adm[i] = BetaHep_adm[i] = Betaff_adm[i] = AlphaHp_adm[i] = AlphaHep_adm[i] = AlphaHepp_adm[i] = Alphad_adm[i] = GammaeH0_adm[i] = GammaeHe0_adm[i] = GammaeHep_adm[i] = 0;
        T = pow(10.0, Tmin_adm + deltaT_adm * i);
        Tfact = 1.0 / (1 + sqrt(T / 1.0e5));
        if(118348. / T < 70.) {BetaH0_adm[i] = 7.5e-19 * exp(-118348 / T) * Tfact;}
        if(473638. / T < 70.) {BetaHep_adm[i] = 5.54e-17 * pow(T, -0.397) * exp(-473638 / T) * Tfact;}

        Betaff_adm[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));
        //AlphaHp_adm[i] = 8.4e-11 * pow(T / 1000, -0.2) / (1. + pow(T / 1.0e6, 0.7)) / sqrt(T);	/* old Cen92 fit */
        //AlphaHep_adm[i] = 1.5e-10 * pow(T, -0.6353); /* old Cen92 fit */
        //AlphaHepp_adm[i] = 4. * AlphaHp_adm[i];	/* old Cen92 fit */
        AlphaHp_adm[i] = 7.982e-11 / ( sqrt(T/3.148) * pow((1.0+sqrt(T/3.148)), 0.252) * pow((1.0+sqrt(T/7.036e5)), 1.748) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        AlphaHep_adm[i]= 9.356e-10 / ( sqrt(T/4.266e-2) * pow((1.0+sqrt(T/4.266e-2)), 0.2108) * pow((1.0+sqrt(T/3.676e7)), 1.7892) ); /* Verner & Ferland (1996) [more accurate than Cen92] */
        AlphaHepp_adm[i] = 2. * 7.982e-11 / ( sqrt(T/(4.*3.148)) * pow((1.0+sqrt(T/(4.*3.148))), 0.252) * pow((1.0+sqrt(T/(4.*7.036e5))), 1.748) ); /* Verner & Ferland (1996) : ~ Z*AlphaHp_adm[1,T/Z^2] */

        if(470000.0 / T < 70) {Alphad_adm[i] = 1.9e-3 * pow(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));}
        if(157809.1 / T < 70) {GammaeH0_adm[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;}
        if(285335.4 / T < 70) {GammaeHe0_adm[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;}
        if(631515.0 / T < 70) {GammaeHep_adm[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;}
    }
}


#ifdef COOL_METAL_LINES_BY_SPECIES

void LoadMultiSpeciesTables_adm(void)
{
    if(All.ComovingIntegrationOn) {
        int i;
        double z;
        if(All.Time==All.TimeBegin) {
            All.SpeciesTableInUse=48;
            ReadMultiSpeciesTables_adm(All.SpeciesTableInUse);
        }
        z=log10(1/All.Time)*48;
        i=(int)z;
        if(i<48) {
            if(i<All.SpeciesTableInUse) {
                All.SpeciesTableInUse=i;
                ReadMultiSpeciesTables_adm(All.SpeciesTableInUse);
            }}
    } else {
        if(All.Time==All.TimeBegin) ReadMultiSpeciesTables_adm(0);
    }
}

void ReadMultiSpeciesTables_adm(int iT)
{
    /* read table w n,T for each species */
    long i_nH=41; long i_Temp=176; long kspecies=(long)NUM_LIVE_SPECIES_FOR_COOLTABLES; long i,j,k,r;
    /* int i_He=7;  int l; */
    FILE *fdcool; char *fname;

    fname=GetMultiSpeciesFilename_adm(iT,0);
    if(ThisTask == 0) printf(" ..opening Cooling Table %s \n",fname);
    if(!(fdcool = fopen(fname, "r"))) {
        printf(" Cannot read species cooling table in file `%s'\n", fname); endrun(456);}
    for(i=0;i<kspecies;i++) {
        for(j=0;j<i_nH;j++) {
            for(k=0;k<i_Temp;k++) {
                r=fread(&SpCoolTable0[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                if(r!=1) {printf(" Reached Cooling EOF! \n");
                }
            }}}
    fclose(fdcool);
    /*
     GetMultiSpeciesFilename(iT,&fname,1);
     if(!(fdcool = fopen(fname, "r"))) {
     printf(" Cannot read species (He) cooling table in file `%s'\n", fname); endrun(456);}
     for(i=0;i<2;i++)
     for(j=0;j<i_nH;j++)
     for(k=0;k<i_Temp;k++)
     for(l=0;l<i_He;l++)
     fread(&SpCoolTable0_He[i][j][k][l],sizeof(float),1,fdcool);
     fclose(fdcool);
     */
    if (All.ComovingIntegrationOn && i<48) {
        fname=GetMultiSpeciesFilename_adm(iT+1,0);
        if(ThisTask == 0) printf(" ..opening (z+) Cooling Table %s \n",fname);
        if(!(fdcool = fopen(fname, "r"))) {
            printf(" Cannot read species 1 cooling table in file `%s'\n", fname); endrun(456);}
        for(i=0;i<kspecies;i++) {
            for(j=0;j<i_nH;j++) {
                for(k=0;k<i_Temp;k++) {
                    r=fread(&SpCoolTable1[i*i_nH*i_Temp + j*i_Temp + k],sizeof(float),1,fdcool);
                    if(r!=1) {printf(" Reached Cooling EOF! \n");
                    }
                }}}
        fclose(fdcool);
        /*
         GetMultiSpeciesFilename(iT+1,&fname,1);
         if(!(fdcool = fopen(fname, "r"))) {
         printf(" Cannot read species 1 (He) cooling table in file `%s'\n", fname); endrun(456);}
         for(i=0;i<2;i++)
         for(j=0;j<i_nH;j++)
         for(k=0;k<i_Temp;k++)
         for(l=0;l<i_He;l++)
         fread(&SpCoolTable1_He[i][j][k][l],sizeof(float),1,fdcool);
         fclose(fdcool);
         */
    }
}

char *GetMultiSpeciesFilename_adm(int i, int hk)
{
    static char fname[100];
    if(i<0) i=0; if(i>48) i=48;
    if(hk==0) {
        sprintf(fname,"./spcool_tables/spcool_%d",i);
    } else {
        sprintf(fname,"./spcool_tables/spcool_He_%d",i);
    }
    return fname;
}

#endif



/* table input (from file TREECOOL) for ionizing parameters */
#define JAMPL_ADM	1.0		/* amplitude factor relative to input table */
#define TABLESIZE_ADM 250		/* Max # of lines in TREECOOL */
static float inlogz_adm[TABLESIZE_ADM];
static double gH0_adm[TABLESIZE_ADM], gHe_adm[TABLESIZE_ADM], gHep_adm[TABLESIZE_ADM]; // upgrade from float to double, should read fine
static double eH0_adm[TABLESIZE_ADM], eHe_adm[TABLESIZE_ADM], eHep_adm[TABLESIZE_ADM]; // upgrade from float to double, should read fine
static int nheattab_adm;		/* length of table */


void ReadIonizeParams_adm(char *fname)
{
    int i; FILE *fdcool;
    if(!(fdcool = fopen(fname, "r"))) {printf(" Cannot read ionization table in file `%s'. Make sure the correct TREECOOL file is placed in the code run-time directory, and that any leading comments (e.g. lines preceded by ##) are deleted from the file.\n", fname); endrun(456);}
    for(i=0; i<TABLESIZE_ADM; i++) {gH0_adm[i]=0;}
    for(i=0; i<TABLESIZE_ADM; i++) {if(fscanf(fdcool, "%g %lg %lg %lg %lg %lg %lg", &inlogz_adm[i], &gH0_adm[i], &gHe_adm[i], &gHep_adm[i], &eH0_adm[i], &eHe_adm[i], &eHep_adm[i]) == EOF) {break;}}
    fclose(fdcool);
    for(i=0, nheattab_adm=0; i<TABLESIZE_ADM; i++) {if(gH0_adm[i] != 0.0) {nheattab_adm++;} else {break;}} /*  nheattab_adm is the number of entries in the table */
    if(ThisTask == 0) printf(" ..read ionization table [TREECOOL] with %d non-zero UVB entries in file `%s'. Make sure to cite the authors from which the UV background was compiled! (See user guide for the correct references).\n", nheattab_adm, fname);
}


void IonizeParams_adm(void)
{
    IonizeParamsTable_adm();
}



void IonizeParamsTable_adm(void)
{
    int i, ilow;
    double logz, dzlow, dzhi;
    double redshift;

    if(All.ComovingIntegrationOn)
        {redshift = 1 / All.Time - 1;}
    else
    {
        /* in non-cosmological mode, still use, but adopt z=0 background */
        redshift = 0;
        /*
         gJHe0_adm = gJHep = gJH0_adm = epsHe0_adm = epsHep_adm = epsH0_adm = J_UV_ADM = 0;
         return;
         */
    }

    logz = log10(redshift + 1.0);
    ilow = 0;
    for(i=0; i<nheattab_adm; i++) {if(inlogz_adm[i] < logz) {ilow = i;} else {break;}}
    dzlow = logz - inlogz_adm[ilow];
    dzhi = inlogz_adm[ilow + 1] - logz;

    if(logz > inlogz_adm[nheattab_adm - 1] || gH0_adm[ilow] == 0 || gH0_adm[ilow + 1] == 0 || nheattab_adm == 0)
    {
        gJHe0_adm = gJHep_adm = gJH0_adm = 0; epsHe0_adm = epsHep_adm = epsH0_adm = 0; J_UV_ADM = 0;
        return;
    }
    else {J_UV_ADM = 1.e-21;}		/* irrelevant as long as it's not 0 */

    gJH0_adm = JAMPL_ADM * pow(10., (dzhi * log10(gH0_adm[ilow]) + dzlow * log10(gH0_adm[ilow + 1])) / (dzlow + dzhi));
    gJHe0_adm = JAMPL_ADM * pow(10., (dzhi * log10(gHe_adm[ilow]) + dzlow * log10(gHe_adm[ilow + 1])) / (dzlow + dzhi));
    gJHep_adm = JAMPL_ADM * pow(10., (dzhi * log10(gHep_adm[ilow]) + dzlow * log10(gHep_adm[ilow + 1])) / (dzlow + dzhi));
    epsH0_adm = JAMPL_ADM * pow(10., (dzhi * log10(eH0_adm[ilow]) + dzlow * log10(eH0_adm[ilow + 1])) / (dzlow + dzhi));
    epsHe0_adm = JAMPL_ADM * pow(10., (dzhi * log10(eHe_adm[ilow]) + dzlow * log10(eHe_adm[ilow + 1])) / (dzlow + dzhi));
    epsHep_adm = JAMPL_ADM * pow(10., (dzhi * log10(eHep_adm[ilow]) + dzlow * log10(eHep_adm[ilow + 1])) / (dzlow + dzhi));

    return;
}


void SetZeroIonization_adm(void)
{
    gJHe0_adm = gJHep_adm = gJH0_adm = 0; epsHe0_adm = epsHep_adm = epsH0_adm = 0; J_UV_ADM = 0;
}


void IonizeParamsFunction_adm(void)
{
    int i, nint;
    double a0, planck, ev, e0_H, e0_He, e0_Hep;
    double gint, eint, t, tinv, fac, eps;
    double at, beta, s;
    double pi;

#define UVALPHA         1.0
    double Jold = -1.0;
    double redshift;

    J_UV_ADM = 0.;
    gJHe0_adm = gJHep_adm = gJH0_adm = 0.;
    epsHe0_adm = epsHep_adm = epsH0_adm = 0.;


    if(All.ComovingIntegrationOn)	/* analytically compute params from power law J_nu */
    {
        redshift = 1 / All.Time - 1;

        if(redshift >= 6) {J_UV_ADM = 0.;}
        else
        {
            if(redshift >= 3) {J_UV_ADM = 4e-22 / (1 + redshift);}
            else
            {
                if(redshift >= 2) {J_UV_ADM = 1e-22;}
                else {J_UV_ADM = 1.e-22 * pow(3.0 / (1 + redshift), -3.0);}
            }
        }
        if(J_UV_ADM == Jold) {return;}
        Jold = J_UV_ADM;
        if(J_UV_ADM == 0) {return;}


        a0 = 6.30e-18;
        planck = 6.6262e-27;
        ev = 1.6022e-12;
        e0_H = 13.6058 * ev;
        e0_He = 24.59 * ev;
        e0_Hep = 54.4232 * ev;

        gint = 0.0;
        eint = 0.0;
        nint = 5000;
        at = 1. / ((double) nint);

        for(i = 1; i <= nint; i++)
        {
            t = (double) i;
            t = (t - 0.5) * at;
            tinv = 1. / t;
            eps = sqrt(tinv - 1.);
            fac = exp(4. - 4. * atan(eps) / eps) / (1. - exp(-2. * M_PI / eps)) * pow(t, UVALPHA + 3.);
            gint += fac * at;
            eint += fac * (tinv - 1.) * at;
        }

        gJH0_adm = a0 * gint / planck;
        epsH0_adm = a0 * eint * (e0_H / planck);
        gJHep_adm = gJH0_adm * pow(e0_H / e0_Hep, UVALPHA) / 4.0;
        epsHep_adm = epsH0_adm * pow((e0_H / e0_Hep), UVALPHA - 1.) / 4.0;

        at = 7.83e-18;
        beta = 1.66;
        s = 2.05;

        gJHe0_adm = (at / planck) * pow((e0_H / e0_He), UVALPHA) *
        (beta / (UVALPHA + s) + (1. - beta) / (UVALPHA + s + 1));
        epsHe0_adm = (e0_He / planck) * at * pow(e0_H / e0_He, UVALPHA) *
        (beta / (UVALPHA + s - 1) + (1 - 2 * beta) / (UVALPHA + s) - (1 - beta) / (UVALPHA + s + 1));

        pi = M_PI;
        gJH0_adm *= 4. * pi * J_UV_ADM;
        gJHep_adm *= 4. * pi * J_UV_ADM;
        gJHe0_adm *= 4. * pi * J_UV_ADM;
        epsH0_adm *= 4. * pi * J_UV_ADM;
        epsHep_adm *= 4. * pi * J_UV_ADM;
        epsHe0_adm *= 4. * pi * J_UV_ADM;
    }
}




void InitCool_adm(void)
{
    if(ThisTask == 0) printf("Initializing adm cooling ...\n");

    All.Time = All.TimeBegin;
    set_cosmo_factors_for_current_time();

#ifdef COOL_GRACKLE
    InitGrackle();
#endif

    InitCoolMemory_adm();
    MakeCoolingTable_adm();
    ReadIonizeParams_adm("TREECOOL");
    IonizeParams_adm();
#ifdef COOL_METAL_LINES_BY_SPECIES
    LoadMultiSpeciesTables_adm();
#endif
}



#ifdef COOL_METAL_LINES_BY_SPECIES
double GetCoolingRateWSpecies_adm(double nHcgs, double logT, double *Z)
{
    double ne_over_nh_tbl=1, Lambda=0;
    int k, N_species_active = (int)NUM_LIVE_SPECIES_FOR_COOLTABLES;

    /* pre-calculate the indices for density and temperature, then we just need to call the tables by species */
    int ixmax=40, iymax=175;
    int ix0, iy0, ix1, iy1;
    double dx, dy, dz, mdz;
    long i_T=iymax+1, inHT=i_T*(ixmax+1);
    if(All.ComovingIntegrationOn && All.SpeciesTableInUse<48) {dz=log10(1/All.Time)*48; dz=dz-(int)dz; mdz=1-dz;} else {dz=0; mdz=1;}

    dx = (log10(nHcgs)-(-8.0))/(0.0-(-8.0))*ixmax;
    dy = (logT-2.0)/(9.0-2.0)*iymax;
    if(dx<0) {dx=0;} else {if(dx>ixmax) {dx=ixmax;}}
    ix0=(int)dx; ix1=ix0+1; if(ix1>ixmax) {ix1=ixmax;}
    dx=dx-ix0;
    if(dy<0) {dy=0;} else {if(dy>iymax) {dy=iymax;}}
    iy0=(int)dy; iy1=iy0+1; if(iy1>iymax) {iy1=iymax;}
    dy=dy-iy0;
    long index_x0y0=iy0+ix0*i_T, index_x0y1=iy1+ix0*i_T, index_x1y0=iy0+ix1*i_T, index_x1y1=iy1+ix1*i_T;

    ne_over_nh_tbl = GetLambdaSpecies_adm(0,index_x0y0,index_x0y1,index_x1y0,index_x1y1,dx,dy,dz,mdz);
    if(ne_over_nh_tbl > 0)
    {
        double zfac = 0.0127 / All.SolarAbundances[0];
        for (k=1; k<N_species_active; k++)
        {
            long k_index = k * inHT;
            Lambda += GetLambdaSpecies_adm(k_index,index_x0y0,index_x0y1,index_x1y0,index_x1y1,dx,dy,dz,mdz) * Z[k+1]/(All.SolarAbundances[k+1]*zfac);
        }
        Lambda /= ne_over_nh_tbl;
    }
    return Lambda;
}


double GetLambdaSpecies_adm(long k_index, long index_x0y0, long index_x0y1, long index_x1y0, long index_x1y1, double dx, double dy, double dz, double mdz)
{
    long x0y0 = index_x0y0 + k_index;
    long x0y1 = index_x0y1 + k_index;
    long x1y0 = index_x1y0 + k_index;
    long x1y1 = index_x1y1 + k_index;
    double i1, i2, j1, j2, w1, w2, u1;
    i1 = SpCoolTable0[x0y0];
    i2 = SpCoolTable0[x0y1];
    j1 = SpCoolTable0[x1y0];
    j2 = SpCoolTable0[x1y1];
    if(dz > 0)
    {
        i1 = mdz * i1 + dz * SpCoolTable1[x0y0];
        i2 = mdz * i2 + dz * SpCoolTable1[x0y1];
        j1 = mdz * j1 + dz * SpCoolTable1[x1y0];
        j2 = mdz * j2 + dz * SpCoolTable1[x1y1];
    }
    w1 = i1*(1-dy) + i2*dy;
    w2 = j1*(1-dy) + j2*dy;
    u1 = w1*(1-dx) + w2*dx;
    return u1;
}

#endif // COOL_METAL_LINES_BY_SPECIES






/* subroutine to update the molecular fraction using our implicit solver for a simple --single-species-- network (just H2) */
void update_explicit_molecular_fraction_adm(int i, double dtime_cgs)
{
    if(dtime_cgs <= 0) {return;}
#ifdef COOL_MOLECFRAC_NONEQM
    // first define a number of environmental variables that are fixed over this update step
    double fH2_initial = SphP[i].MolecularMassFraction_perNeutralH; // initial molecular fraction per H atom, entering this subroutine, needed for update below
    double xn_e=1, nh0=0, nHe0, nHepp, nhp, nHeII, temperature, mu_meanwt=1, rho=SphP[i].Density*All.cf_a3inv, u0=SphP[i].InternalEnergyPred;
    temperature = ThermalProperties_adm(u0, rho, i, &mu_meanwt, &xn_e, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties [will assume fixed value of fH2 at previous update value]
    double T=1, Z_Zsol=1, urad_G0=1, xH0=1, x_e=0, nH_cgs=rho*UNIT_DENSITY_IN_NHCGS; // initialize definitions of some variables used below to prevent compiler warnings
    Z_Zsol=1; urad_G0=1; // initialize metal and radiation fields. will assume solar-Z and spatially-uniform Habing field for incident FUV radiation unless reset below.
    if(temperature > 3.e5) {SphP[i].MolecularMassFraction_perNeutralH=SphP[i].MolecularMassFraction=0; return;} else {T=temperature;} // approximations below not designed for high temperatures, should simply give null
    xH0 = DMIN(DMAX(nh0, 0.),1.); // get neutral fraction [given by call to this program]
    if(xH0 <= MIN_REAL_NUMBER) {SphP[i].MolecularMassFraction_perNeutralH=SphP[i].MolecularMassFraction=0; return;} // no neutral gas, no molecules!
    x_e = DMIN(DMAX(xn_e, 0.),2.); // get free electron ratio [number per H nucleon]
    double gamma_12=return_local_gammamultiplier_adm(i)*gJH0_adm/1.0e-12, shieldfac=return_uvb_shieldfac_adm(i,gamma_12,nH_cgs,log10(T)), urad_from_uvb_in_G0=sqrt(shieldfac)*(gJH0_adm/2.29e-10); // estimate UVB contribution if we have partial shielding, to full photo-dissociation rates //
#ifdef METALS
    Z_Zsol = P[i].Metallicity[0]/All.SolarAbundances[0]; // metallicity in solar units [scale to total Z, since this mixes dust and C opacity], and enforce a low-Z floor to prevent totally unphysical behaviors at super-low Z [where there is still finite opacity in reality; e.g. Kramer's type and other opacities enforce floor around ~1e-3]
#endif
    /* get incident radiation field from whatever module we are using to track it */
#if defined(RT_PHOTOELECTRIC) || defined(RT_LYMAN_WERNER)
    int whichbin = RT_FREQ_BIN_LYMAN_WERNER;
#if !defined(RT_LYMAN_WERNER)
    whichbin = RT_FREQ_BIN_PHOTOELECTRIC; // use photo-electric bin as proxy (very close) if don't evolve LW explicitly
#endif
    urad_G0 = SphP[i].Rad_E_gamma[whichbin] * (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_CGS / 3.9e-14; // convert to Habing field //
#endif
    urad_G0 += urad_from_uvb_in_G0; // include whatever is contributed from the meta-galactic background, fed into this routine
    urad_G0 = DMIN(DMAX( urad_G0 , 1.e-10 ) , 1.e10 ); // limit values, because otherwise exponential self-shielding approximation easily artificially gives 0 incident field
    // define a number of variables needed in the shielding module
    double dx_cell = Get_Particle_Size(i) * All.cf_atime; // cell size
    double surface_density_H2_0 = 5.e14 * PROTONMASS, x_exp_fac=0.00085, w0=0.2; // characteristic cgs column for -molecular line- self-shielding
    w0 = 0.035; // actual calibration from Drain, Gnedin, Richings, others: 0.2 is more appropriate as a re-calibration for sims doing local eqm without ability to resolve shielding at higher columns
    //double surface_density_local = xH0 * SphP[i].Density * All.cf_a3inv * dx_cell * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth through the local cell/slab. note G0 is -already- attenuated in the pre-processing by dust.
    double surface_density_local = xH0 * evaluate_NH_from_GradRho(P[i].GradRho,PPP[i].Hsml,SphP[i].Density,PPP[i].NumNgb,1,i) * UNIT_SURFDEN_IN_CGS; // this is -just- the [neutral] depth to infinity with our Sobolev-type approximation. Note G0 is already attenuated by dust, but we need to include H2 self-shielding, for which it is appropriate to integrate to infinity.
    double v_thermal_rms = 0.111*sqrt(T); // sqrt(3*kB*T/2*mp), since want rms thermal speed of -molecular H2- in kms
    double dv2=0; int j,k; for(j=0;j<3;j++) {for(k=0;k<3;k++) {double vt = SphP[i].Gradients.Velocity[j][k]*All.cf_a2inv; /* physical velocity gradient */
        if(All.ComovingIntegrationOn) {if(j==k) {vt += All.cf_hubble_a;}} /* add hubble-flow correction */
        dv2 += vt*vt;}} // calculate magnitude of the velocity shear across cell from || grad -otimes- v ||^(1/2)
    double dv_turb=sqrt(dv2)*dx_cell*UNIT_VEL_IN_KMS; // delta-velocity across cell
    double x00 = surface_density_local / surface_density_H2_0, x01 = x00 / (sqrt(1. + 3.*dv_turb*dv_turb/(v_thermal_rms*v_thermal_rms)) * sqrt(2.)*v_thermal_rms), y_ss, x_ss_1, x_ss_sqrt, fH2_tmp, fH2_max, fH2_min, Q_max, Q_min, Q_initial, Q_0, Q_1, fH2_0, fH2_1, fH2_new; // variable needed below. note the x01 term corrects following Gnedin+Draine 2014 for the velocity gradient at the sonic scale, assuming a Burgers-type spectrum [their Eq. 3]
    double b_time_Mach = 0.5 * dv_turb / (v_thermal_rms/sqrt(3.)); // cs_thermal for molecular [=rms v_thermal / sqrt(3)], dv_turb to full inside dx, assume "b" prefactor for compressive-to-solenoidal ratio corresponding to the 'natural mix' = 0.5. could further multiply by 1.58 if really needed to by extended dvturb to 2h = H, and vthermal from molecular to atomic for the generating field, but not as well-justified
    double clumping_factor = 1. + b_time_Mach*b_time_Mach; // this is the exact clumping factor for a standard lognormal PDF with S=ln[1+b^2 Mach^2] //
    double clumping_factor_3 = clumping_factor*clumping_factor*clumping_factor; // clumping factor N for <rho^n>/<rho>^n = clumping factor^(N*(N-1)/2) //

    /* evolve dot[nH2]/nH0 = d_dt[fH2[neutral]] = (1/nH0) * (a_H2*rho_dust*nHI [dust formation] + a_GP*nHI*ne [gas-phase formation] + b_3B*nHI*nHI*(nHI+nH2/8) [3-body collisional form] - b_H2HI*nHI*nH2 [collisional dissociation]
        - b_H2H2*nH2*nH2 [collisional mol-mol dissociation] - Gamma_H2^LW * nH2 [photodissociation] - Gamma_H2^+ [photoionization] - xi_H2*nH2 [CR ionization/dissociation] ) */
    double fH2=0, sqrt_T=sqrt(T), nH0=xH0*nH_cgs, n_e=x_e*nH_cgs, EXPmax=40.; int iter=0; // define some variables for below, including neutral H number density, free electron number, etc.
    double a_Z  = (9.e-19 * T / (1. + 0.04*sqrt_T + 0.002*T + 8.e-6*T*T)) * (0.5*Z_Zsol) * nH_cgs * clumping_factor; // dust formation
    double a_GP = (1.833e-21 * pow(T,0.88)) * n_e * clumping_factor; // gas-phase formation
    double b_3B = (6.0e-32/sqrt(sqrt_T) + 2.0e-31/sqrt_T) * nH0 * nH0 * clumping_factor_3; // 3-body collisional formation
    double b_H2HI = (7.073e-19 * pow(T,2.012) * exp(-DMIN(5.179e4/T,EXPmax)) / pow(1. + 2.130e-5*T , 3.512)) * (nH0/2.) * clumping_factor; // collisional dissociation
    double b_H2H2 = (5.996e-30 * pow(T,4.1881) * exp(-DMIN(5.466e4/T,EXPmax)) / pow(1. + 6.761e-6*T , 5.6881)) * (nH0/2.) * (1./2.) * clumping_factor; // collisional mol-mol dissociation
    double G_LW = 3.3e-11 * urad_G0 * (1./2.); // photo-dissociation (+ionization); note we're assuming a spectral shape identical to the MW background mean, scaling by G0
    double xi_cr_H2 = (7.525e-16) * (1./2.), prefac_CR=1.;; // CR dissociation (+ionization)
    if(All.ComovingIntegrationOn) {double rhofac = (rho*UNIT_DENSITY_IN_CGS)/(1000.*COSMIC_BARYON_DENSITY_CGS);
        if(rhofac < 0.2) {prefac_CR=0;} else {if(rhofac > 200.) {prefac_CR=1;} else {prefac_CR=exp(-1./(rhofac*rhofac));}}} // in cosmological runs, turn off CR heating for any gas with density unless it's >1000 times the cosmic mean density
    xi_cr_H2 *= prefac_CR;

    // want to solve the implicit equation: f_f = f_0 + g[f_f]*dt, where g[f_f] = df_dt evaluated at f=f_f, so root-find: dt*g[f_f] + f_0-f_f = 0
    // can write this as a quadtratic: x_a*f^2 - x_b_0*f - xb_LW*f + x_c = 0, where xb_LW is a non-linear function of f accounting for the H2 self-shielding terms
    double G_LW_dt_unshielded = G_LW * dtime_cgs; // LW term without shielding, multiplied by timestep for dimensions needed below
    double x_a = (b_3B + b_H2HI - b_H2H2) * dtime_cgs; // terms quadratic in f -- this term can in principle be positive or negative, usually positive
    double x_b_0 = (a_GP + a_Z + 2.*b_3B + b_H2HI + xi_cr_H2) * dtime_cgs + 1.; // terms linear in f [note sign, pulling the -sign out here] -- positive-definite
    double x_c = (a_GP + a_Z + b_3B) * dtime_cgs + fH2_initial; // terms independent of f -- positive-definite
    double y_a = x_a / (x_c + MIN_REAL_NUMBER), x_b, y_b, z_a; // convenient to convert to dimensionless variable needed for checking definite-ness
    // use the previous-timestep value of fH2 to guess the shielding term and then compute the resulting fH2
    fH2_tmp=fH2_initial; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
    Q_initial = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the previous value of fH2
    z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant

    /* now comes the tricky bit -- need to account for non-linear part of the solution for the molecular line self-shielding [depends on fH2, not just the dust external shielding already accounted for */
    if((fH2 > 1.e-10) && (fH2 < 1) && (G_LW_dt_unshielded > 0.1*x_b_0)) // fH2 is non-trivial, and the radiation term is significant, so to get an accurate update we need to invoke a non-linear solver here
    {
        // set updated guess values
        fH2_tmp=fH2; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
        fH2_1 = fH2; Q_1 = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the updated value of fH2

        // set lower values for bracketing
        x_b=x_b_0+G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // recalculate all terms that depend on the shielding
        fH2_min = DMAX(0,DMIN(1,fH2)); // this serves as a lower-limit for fH2
        fH2_tmp=fH2_min; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
        Q_min = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the updated value of fH2

        // set upper values for bracketing
        fH2_tmp=1.; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // recalculate all terms that depend on the shielding
        z_a=4.*y_a/(y_b*y_b + MIN_REAL_NUMBER); if(z_a>1.) {fH2=1.;} else {if(fabs(z_a)<0.1) {fH2=(1.+0.25*z_a*(1.+0.5*z_a))/(y_b + MIN_REAL_NUMBER);} else {fH2=(2./(y_b + MIN_REAL_NUMBER))*(1.-sqrt(1.-z_a))/z_a;}} // calculate f assuming the shielding term is constant
        fH2_max = DMAX(0,DMIN(1,fH2)); // this serves as an upper-limit for fH2
        fH2_tmp=fH2_max; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
        Q_max = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the updated value of fH2

        if(fH2_1 < fH2_min) {fH2 = fH2_min;} // hitting lower bound already, set to that value and exit
        else if(fH2_1 > fH2_max) {fH2 = fH2_max;} // hitting upper bound already, set to that value and exit
        else if(Q_min*Q_max >= 0) // bracketing indicates that in this timestep, we will move fully to one or the other limit -- so do that, and don't need to iterate!
            {if(fabs(Q_min) < fabs(Q_max)) {if(fabs(Q_min) < fabs(Q_1)) {fH2 = fH2_min;}} else {if(fabs(Q_max) < fabs(Q_1)) {fH2 = fH2_max;}}} // decide if Qmin/max corresponds more closely to desired zero, so move to fH2 matching that value
        else if(fH2_max > 1.01*fH2_min) // worth attempting to iterate
        {
            Q_0 = Q_initial; fH2_0 = fH2_initial; // first iteration step is already done in all the above
            fH2 = exp( (log(fH2_0)*Q_1 - log(fH2_1)*Q_0) / (Q_1-Q_0) ); // do a Newton-Raphson step in log[f_H2] space now that we have good initial brackets
            if(fH2 > fH2_max) // ok we overshot the upper limit, test if bracketing works between fH2 and fH2_max, otherwise just use fH2_max
                {if(Q_1*Q_max < 0) {fH2 = exp( (log(fH2_max)*Q_1 - log(fH2_1)*Q_max) / (Q_1-Q_max) );} else {fH2=fH2_max;}}
            else if(fH2 < fH2_min)
                {if(Q_1*Q_min < 0) {fH2 = exp( (log(fH2_min)*Q_1 - log(fH2_1)*Q_min) / (Q_1-Q_min) );} else {fH2=fH2_min;}}
            else while(1)
            {
                    fH2_tmp=fH2_1; x_ss_1=1.+fH2_tmp*x01; x_ss_sqrt=sqrt(1.+fH2_tmp*x00); y_ss=(1.-w0)/(x_ss_1*x_ss_1) + w0/x_ss_sqrt*exp(-DMIN(EXPmax,x_exp_fac*x_ss_sqrt)); x_b=x_b_0+y_ss*G_LW_dt_unshielded; y_b=x_b/(x_c + MIN_REAL_NUMBER); // calculate all the terms we need to solve for the zeros of this function
                    Q_1 = 1 + y_a*fH2_tmp*fH2_tmp - y_b*fH2_tmp; // value of the function we are trying to zero, with the updated value of fH2
                    //if(Q_1*Q_0 >= 0) {break;} // no longer bracketing, end while loop
                    fH2_new = exp( (log(fH2_0)*Q_1 - log(fH2_1)*Q_0) / (Q_1-Q_0) ); fH2_0=fH2_1; Q_0=Q_1; fH2_1=fH2_new; // update guess and previous values //
                    iter++; // count iterations
                    if(fabs(fH2_1-fH2_0) < 0.01*(0.5*(fH2_1+fH2_0))) {break;} // converged well enough for our purposes!
                    if((y_ss > 0.95) || (y_ss*G_LW_dt_unshielded < 0.1*x_b)) {break;} // negligible shielding, or converged to point where external LW is not dominant dissociator so no further iteration needed
                    if((fH2 > 0.95*fH2_max) || (fH2 > 0.99) || (fH2 < 1.e-10) || (fH2 < 1.05*fH2_min) || (iter > 10)) {break;} // approached physical limits or bounds of validity, or end of iteration cycle
            } // end of convergence iteration to find solution for fmol
        } // opened plausible iteration clause
    } // opened self-shielding clause [attempting to bracket]
    if(!isfinite(fH2)) {fH2=0;} else {if(fH2>1) {fH2=1;} else if(fH2<0) {fH2=0;}} // check vs nans, valid values
    SphP[i].MolecularMassFraction_perNeutralH = fH2; // record variable -- this will be used for the next update, meanwhile the total fraction will be used in various routines through the code
    SphP[i].MolecularMassFraction = xH0 * SphP[i].MolecularMassFraction_perNeutralH; // record variable -- this is largely what is needed below
#endif
}





/* simple subroutine to estimate the dust temperatures in our runs without detailed tracking of these individually [more detailed chemistry models do this] */
double get_equilibrium_dust_temperature_estimate_adm(int i, double shielding_factor_for_exgalbg)
{   /* simple three-component model [can do fancier] with cmb, dust, high-energy photons */
#if defined(RT_INFRARED)
    if(i >= 0) {return SphP[i].Dust_Temperature;} // this is pre-computed -- simply return it
#endif
    double e_CMB=0.262*All.cf_a3inv/All.cf_atime, T_cmb=2.73/All.cf_atime; // CMB [energy in eV/cm^3, T in K]
    double e_IR=0.31, Tdust_ext=DMAX(30.,T_cmb); // Milky way ISRF from Draine (2011), assume peak of dust emission at ~100 microns
    double e_HiEgy=0.66, T_hiegy=5800.; // Milky way ISRF from Draine (2011), assume peak of stellar emission at ~0.6 microns [can still have hot dust, this effect is pretty weak]
    if(i >= 0)
    {
#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
        e_HiEgy=0; e_IR = 0; int k; double E_tot_to_evol_eVcgs = (SphP[i].Density*All.cf_a3inv/P[i].Mass) * UNIT_PRESSURE_IN_EV;
        for(k=0;k<N_RT_FREQ_BINS;k++) {e_HiEgy+=SphP[i].Rad_E_gamma_Pred[k];}
#if defined(RT_INFRARED)
        e_IR += SphP[i].Rad_E_gamma_Pred[RT_FREQ_BIN_INFRARED]; Tdust_ext = SphP[i].Radiation_Temperature; // note IR [irrelevant b/c of call above, but we'll keep this as a demo]
#endif
        e_HiEgy -= e_IR; // don't double-count the IR component flagged above //
        e_IR *= E_tot_to_evol_eVcgs; e_HiEgy *= E_tot_to_evol_eVcgs;
#endif
    }
    e_HiEgy += shielding_factor_for_exgalbg * 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // this comes from the cosmic optical+UV backgrounds. small correction, so treat simply, and ignore when self-shielded.
    double Tdust_eqm = 10.; // arbitrary initial value //
    if(Tdust_ext*e_IR < 1.e-10 * (T_cmb*e_CMB + T_hiegy*e_HiEgy)) { // IR term is totally negligible [or zero exactly], use simpler expression assuming constant temperature for it to avoid sensitivity to floating-pt errors //
        Tdust_eqm = 2.92 * pow(Tdust_ext*e_IR + T_cmb*e_CMB + T_hiegy*e_HiEgy, 1./5.); // approximate equilibrium temp assuming Q~1/lambda [beta=1 opacity law], assuming background IR temp is a fixed constant [relevant in IR-thin limit, but we don't know T_rad, so this is a guess anyways]
    } else { // IR term is not vanishingly small. we will assume the IR radiation temperature is equal to the local Tdust. lacking any direct evolution of that field, this is a good proxy, and exact in the locally-IR-optically-thick limit. in the locally-IR-thin limit it slightly under-estimates Tdust, but usually in that limit the other terms dominate anyways, so this is pretty safe //
        double T0=2.92, q=pow(T0*e_IR,0.25), y=(T_cmb*e_CMB + T_hiegy*e_HiEgy)/(T0*e_IR*q); if(y<=1) {Tdust_eqm=T0*q*(0.8+sqrt(0.04+0.1*y));} else {double y5=pow(y,0.2), y5_3=y5*y5*y5, y5_4=y5_3*y5; Tdust_eqm=T0*q*(1.+15.*y5_4+sqrt(1.+30.*y5_4+25.*y5_4*y5_4))/(20.*y5_3);} // this gives an extremely accurate and exactly-joined solution to the full quintic equation assuming T_rad_IR=T_dust
    }
    return DMAX(DMIN(Tdust_eqm , 2000.) , 1.); // limit at sublimation temperature or some very low temp //
}




/* this function evaluates Compton heating+cooling rates and synchrotron cooling for thermal gas populations, accounting for the
    explicitly-evolved radiation field if it is evolved (otherwise assuming a standard background), and B-fields if they
    are evolved, as well as the proper relativistic or non-relativistic effects and two-temperature plasma effects. */
double evaluate_Compton_heating_cooling_rate_adm(int target, double T, double nHcgs, double n_elec, double shielding_factor_for_exgalbg)
{
    double Lambda = 0;
    double compton_prefac_eV = 2.16e-35 / nHcgs; // multiply field in eV/cm^3 by this and temperature difference to obtain rate

    double e_CMB_eV=0.262*All.cf_a3inv/All.cf_atime, T_cmb = 2.73/All.cf_atime; // CMB [energy in eV/cm^3, T in K]
    Lambda += compton_prefac_eV * n_elec * e_CMB_eV * (T-T_cmb);

    double e_UVB_eV = shielding_factor_for_exgalbg * 7.8e-3 * pow(All.cf_atime,3.9)/(1.+pow(DMAX(-1.+1./All.cf_atime,0.001)/1.7,4.4)); // this comes from the cosmic optical+UV backgrounds. small correction, so treat simply, and ignore when self-shielded.
    Lambda += compton_prefac_eV * n_elec * e_UVB_eV * (T-2.e4); // assume very crude approx Compton temp ~2e4 for UVB

#if defined(RADTRANSFER) || defined(RT_USE_GRAVTREE_SAVE_RAD_ENERGY) // use actual explicitly-evolved radiation field, if possible
    if(target >= 0)
    {
        int k; double E_tot_to_evol_eVcgs = (SphP[target].Density*All.cf_a3inv/P[target].Mass) * UNIT_PRESSURE_IN_EV;
        for(k=0;k<N_RT_FREQ_BINS;k++)
        {
            double e_tmp = SphP[target].Rad_E_gamma_Pred[k] * E_tot_to_evol_eVcgs, Teff = 0;

#if defined(RT_INFRARED) /* special mid-through-far infrared band, which includes IR radiation temperature evolution */
            if(k==RT_FREQ_BIN_INFRARED) {Teff=SphP[target].Dust_Temperature;}
#endif
#if defined(RT_OPTICAL_NIR) /* Optical-NIR approximate spectra for stars as used in the FIRE (Hopkins et al.) models; from 0.41-3.4 eV */
            if(k==RT_FREQ_BIN_OPTICAL_NIR) {Teff=2800.;}
#endif
#if defined(RT_NUV) /* Near-UV approximate spectra (UV/optical spectra, sub-photo-electric, but high-opacity) for stars as used in the FIRE (Hopkins et al.) models; from 3.4-8 eV */
            if(k==RT_FREQ_BIN_NUV) {Teff=12000.;}
#endif
#if defined(RT_PHOTOELECTRIC) /* photo-electric bands (8-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
            if(k==RT_FREQ_BIN_PHOTOELECTRIC) {Teff=24400.;}
#endif
#if defined(RT_LYMAN_WERNER) /* lyman-werner bands (11.2-13.6 eV, specifically): below is from integrating the spectra from STARBURST99 with the Geneva40 solar-metallicity + lower tracks */
            if(k==RT_FREQ_BIN_LYMAN_WERNER) {Teff=28800.;}
#endif
#if defined(RT_CHEM_PHOTOION) /* Hydrogen and Helium ionizing bands: H0 here */
            if(k==RT_FREQ_BIN_H0) {Teff=2340.*rt_nu_eff_eV[k];}
#endif
#if defined(RT_PHOTOION_MULTIFREQUENCY) /* Hydrogen and Helium ionizing bands: He bands */
            if(k==RT_FREQ_BIN_He0 || k==RT_FREQ_BIN_He1 || k==RT_FREQ_BIN_He2) {Teff=2340.*rt_nu_eff_eV[k];}
#endif
#if defined(RT_SOFT_XRAY) /* soft and hard X-rays for e.g. compton heating by X-ray binaries */
            if(k==RT_FREQ_BIN_SOFT_XRAY) {Teff=3.6e6;}
#endif
#if defined(RT_HARD_XRAY) /* soft and hard X-rays for e.g. compton heating by X-ray binaries */
            if(k==RT_FREQ_BIN_HARD_XRAY) {Teff=1.7e7;}
#endif
            if(Teff < 3.e4) {e_tmp *= n_elec;} // low-energy radiation acts inefficiently on neutrals here
            Lambda += compton_prefac_eV * e_tmp * (T - Teff); // add to compton heating/cooling terms
        }
    }
#else // no explicit RHD terms evolved, so assume a MW-like ISRF instead
    double e_IR_eV=0.31, T_IR=DMAX(30.,T_cmb); // Milky way ISRF from Draine (2011), assume peak of dust emission at ~100 microns
    double e_OUV_eV=0.66, T_OUV=5800.; // Milky way ISRF from Draine (2011), assume peak of stellar emission at ~0.6 microns [can still have hot dust, this effect is pretty weak]
    Lambda += compton_prefac_eV * n_elec * (e_IR_eV*(T-T_IR) + e_OUV_eV*(T-T_OUV));
#endif


#ifdef MAGNETIC /* include sychrotron losses as well as long as we're here, since these scale more or less identically just using the magnetic instead of radiation energy */
    if(target >= 0)
    {
        double b_muG = get_cell_Bfield_in_microGauss(target), U_mag_ev=0.0248342*b_muG*b_muG;
        Lambda += compton_prefac_eV * U_mag_ev * T; // synchrotron losses proportional to temperature (non-relativistic here), as inverse compton, just here without needing to worry about "T-T_eff", as if T_eff->0
    }
#endif

    double T_eff_for_relativistic_corr = T; /* used below, but can be corrected */
    if(Lambda > 0) /* per CAFG's calculations, we should note that at very high temperatures, the rate-limiting step may be the Coulomb collisions moving energy from protons to e-; which if slow will prevent efficient e- cooling */
    {
        double Lambda_limiter_var = 1.483e34 * Lambda*Lambda*T; /* = (Lambda/2.6e-22)^2 * (T/1e9): if this >> 1, follow CAFG and cap at cooling rate assuming equilibrium e- temp from Coulomb exchange balancing compton */
        if(Lambda_limiter_var > 1) {Lambda_limiter_var = 1./pow(Lambda_limiter_var,0.2); Lambda*=Lambda_limiter_var; T_eff_for_relativistic_corr*=Lambda_limiter_var;}
    }
    if(T_eff_for_relativistic_corr > 3.e7) {Lambda *= (T_eff_for_relativistic_corr/1.5e9) / (1-exp(-T_eff_for_relativistic_corr/1.5e9));} /* relativistic correction term, becomes important at >1e9 K, enhancing rate */

    return Lambda;
}







/*  this function computes the self-consistent temperature and electron fraction */
double ThermalProperties_adm(double u, double rho, int target, double *mu_guess, double *ne_guess, double *nH0_guess, double *nHp_guess, double *nHe0_guess, double *nHep_guess, double *nHepp_guess)
{
    if(target >= 0) {*ne_guess=SphP[target].Ne; *nH0_guess = DMAX(0,DMIN(1,1.-( *ne_guess / 1.2 )));} else {*ne_guess=1.; *nH0_guess=0.;}
    rho *= UNIT_DENSITY_IN_CGS; u *= UNIT_SPECEGY_IN_CGS;   /* convert to physical cgs units */
    double temp = convert_u_to_temp_adm(u, rho, target, ne_guess, nH0_guess, nHp_guess, nHe0_guess, nHep_guess, nHepp_guess, mu_guess);
    return temp;
}



/* function to return the local multiplier relative to the UVB model to account in some local RHD models for local ionizing sources */
double return_local_gammamultiplier_adm(int target)
{
    return 1;
}


/* function to attenuate the UVB to model self-shielding in optically-thin simulations */
double return_uvb_shieldfac_adm(int target, double gamma_12, double nHcgs, double logT)
{
#ifdef GALSF_EFFECTIVE_EQS
    return 1; // self-shielding is implicit in the sub-grid model already //
#endif
    double NH_SS_z, NH_SS = 0.0123; /* CAFG: H number density above which we assume no ionizing bkg (proper cm^-3) */
    if(gamma_12>0) {NH_SS_z = NH_SS*pow(gamma_12,0.66)*pow(10.,0.173*(logT-4.));} else {NH_SS_z = NH_SS*pow(10.,0.173*(logT-4.));}
    double q_SS = nHcgs/NH_SS_z;
#ifdef COOLING_SELFSHIELD_TESTUPDATE_RAHMATI
    return 0.98 / pow(1 + pow(q_SS,1.64), 2.28) + 0.02 / pow(1 + q_SS*(1.+1.e-4*nHcgs*nHcgs*nHcgs*nHcgs), 0.84); // from Rahmati et al. 2012: gives gentler cutoff at high densities. but we need to modify it with the extra 1+(nHcgs/10)^4 denominator term since at very high nH, this cuts off much too-slowly (as nH^-0.84), which means UVB heating can be stronger than molecular cooling even at densities >> 1e4
#else
    return 1./(1.+q_SS*(1.+q_SS/2.*(1.+q_SS/3.*(1.+q_SS/4.*(1.+q_SS/5.*(1.+q_SS/6.*q_SS)))))); // this is exp(-q) down to ~1e-5, then a bit shallower, but more than sufficient approximation here //
#endif
}


#endif //ADM
#endif // COOLING
