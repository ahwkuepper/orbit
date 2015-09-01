/***************************************************************************
 *   Copyright (C) 2015 by Andreas H.W. Kuepper                            *
 *   ahwkuepper@gmail.com                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/***************************************************************************
 *   Compile using the command: cc -o orbit orbit.c -lm                    *
 ***************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>


//constants
#define Pi 3.14159265
#define PI 3.14159265
#define G  0.0043009211           //gravitational constant in [km^2/s^2/Msun*pc]
#define SMALL 1.0E-3
#define SUPERSMALL -1.E50


//functions
int rk4_drv(double *t, double tmax, double dtout, double mdiff, double *x, double *v, double vorz, double *parameter, double *fit, int newspaper);
int rk4_drvtail(double *t, double tmax, double mdiff, double *x, double *v, double vorz, double *xc, double *vc, double *parameter);
void getforce(double *x, double *v, double *a, double *parameter);
void getforcetail(double *x, double *v, double *a, double *xc, double *parameter);
void do_step(double dt, double *x, double *v, double *parameter);
void do_steptail(double dt, double *x, double *v, double *xc, double *vc, double *parameter);
double get_gauss(void);
void convert(double *x, double *v, double *dsun, double *vrsun, double *vr, double *l, double *b, double *lcosb, double *RA, double *DEC, double *mu_alpha, double *mu_alphacosdelta, double *mu_delta, double *mutemp, double *PAtemp, int coordtype_coco, int vcoordtype, int radiococo, double vLSRtemp, double rgalsun);
double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);



//IMPORTANT PARAMETERS
double const dtout = 10.0;              //time step for output [Myr]
double tpast = -6000.0;                 //time at beginning of integration [Myr]
double const tstart = 0.0;              //time at input of cluster coordinates [Myr], usually today, i.e. 0.0
double const tfuture = 0.0;             //time at end of integration [Myr], change if you want to stop the integration at a differnt time than today
double const mass_cluster = 20000.0;    //mass of cluster at the present day [Msun]
double const mass_loss_rate = 0.0;      //mass loss rate [Msun/Myr]
double const R4 = 20.0;                 //cluster plummer radius [pc]
double const vrexcess = 2.0;            //mean escape velocity ~ 0.5-2 times velocity dispersion [km/s]
double const rexcess = 0.25;            //displacement from Lagrange point [rtide]
int const nperstep = 1;                 //tail particles generated per timestep

double xinitial[] = {8330.0, 0.0, 0.0}; //Initial position (at present day) of satellite [pc]
double vinitial[] = {0.0, 239.5, 0.0};  //Initial velocity (at present day) of satellite [km/s]



//integration parameters
int const radio = 0;                //say what you're doing
double const mdiff = 1.E-4;         //precission
double const dt0 = 1.E-4;			//initial time-step [Myr]
double const dtmax = 50.0;          //maximum time-step [Myr]
double const Rgalmin = 10.0;        //minimum galactocentric radius
double const Rgalmax = 1.0e10;      //maximum galactocentric radius
double sigma[] = {0.0, 0.0, 0.0};   //velocity scatter in x, y, z [km/s]

//potential parameters
int const gpot = 3;             //type of Galactic potential (1= Allen & Santillan (1991), 2= flattened log-halo (Koposov et al. 2010), 3= flattened NFW with Johnston et al. (2005) disk & bulge, 4= triaxial log-halo from Law & Majewski (2010) with Johnston et al. (2005) disk & bulge

//Allen & Santillan potential constants
double const b1 = 230.0;        //I12 //[pc]
double const M1 = 9.4888e09;    //I12 //[solar masses]
double const a2 = 4220.0;       //I12 //[pc]
double const b2 = 292.0;        //I12 //[pc]
double const M2 = 6.62592e10;   //I12 //[solar masses]
double const a3 = 6834.07;      //I12 //[pc]
double const M3 = 2.36176e10;   //I12 //[solar masses]
double const qz = 1.0;          //halo flattening along z axis

//Law, Majewski & Johnston (2009) potential constants
double const b1_LMJ = 700.0;        //[pc]
double const M1_LMJ = 3.4e10;       //[solar masses]
double const a2_LMJ = 6500.0;       //[pc]
double const b2_LMJ = 260.0;        //[pc]
double const M2_LMJ = 1.0e11;       //[solar masses]
double const q_halo = 0.946292;     //flattening of NFW halo
double const r_halo = 37857.300000; //scale radius of NFW halo
double const mass_halo = 1.584460e+12; //scale mass of NFW halo

//Law & Majweski (2010) potential parameters
double c_lm = 12000.0;                 //log-halo scale radius [pc]
double vh2_lm = 172.333236*172.333236; //log-halo scale velocity (squared) [(km/s)^2]
double q1_lm = 1.38;                   //flattenings
double q2_lm = 1.0;
double qz_lm = 1.36;
double phi_lm = 1.692969;              //position angle [rad]

//Galactic North Pole parameters
double const alphaGNP = 192.859508; //Galactic north pole in J2000 coordinates [deg]
double const deltaGNP = 27.128336;
double const PAGNP = 122.932;       //Position angle with respect to equatorial pole [deg]

//solar parameters
double const rgalsun = 8330.0;      //solar Galactocentric radius [pc] (standard = 8330.0; Gillessen et al. 2009)
double const vLSR = 239.5;          //rotational velocity of local standard of rest (Reid & Brunthaler 2004 - vysun)
double const vxsun = 11.1;          //+0.69/−0.75 - solar motion with respect to the LSR from Schoenrich, Binney & Dehnen (2010) [km/s]
double const vysun = 12.24;         //+0.47−0.47 //(24.0 = Bovy et al. (2012))
double const vzsun = 7.25;          //+0.37−0.36

//cluster parameters
int const tails = 1;                //integrate tidal tail test particles (0= no, 1= yes);
double const Rstop = 0.0;           //increase r_edge if test particles gets inside r < Rstop, set 0 for no redge parameterscan, else e.g. 20 pc
double const redge = 20.0;          //cluster edge radius = minimum cluster radius a particle is allowed to escape from [pc]
double const rtidemax = 350.0;      //maximum value for rtide





int main() {

    double *fit;
    int newspaper = 1;
    int err = 0;

    //round simulation time to multiple of output time
    int ratio;
    ratio = (int) 1.0*tpast/dtout;
    tpast = 1.0*ratio*dtout;
    
    //put model parameters in an array for passing it to the integrator
    double parameter[13];
    parameter[0] = mass_cluster;            //mass of cluster
    parameter[1] = 0.0;                     //dummy
    parameter[2] = 0.0;                     //dummy
    parameter[3] = mass_halo;               //mass of DM halo
    parameter[4] = 0.0;                     //dummy
    parameter[6] = rgalsun;                 //distance sun-galactic center
    parameter[7] = vLSR;                    //y-velocity of LSR
    parameter[8] = mass_loss_rate;          //cluster mass loss rate
    parameter[10] = q_halo;                 //flattening of dark halo
    parameter[11] = r_halo;                 //concentration of dark halo
    parameter[12] = tpast;                  //integration time
    
    //present-day coordinates of cluster in Cartesian frame
	double x[3], v[3];
    x[0] = xinitial[0];
    x[1] = xinitial[1];
    x[2] = xinitial[2];
    v[0] = vinitial[0];
    v[1] = vinitial[1];
    v[2] = vinitial[2];
    
    //get position of cluster at t = -tpast
    double sign = -1.0;
    double tmax = tpast;
    double dtoutt = tpast;
    double t = tstart;
    double tspan = sqrt(pow(tpast-tstart,2));
    rk4_drv(&t,tmax,dtoutt,mdiff,&x[0],&v[0],sign,parameter,fit,newspaper);
    
    //save whatever initial conditions for possible restarts
    double xrestart[3], vrestart[3], trestart;
    xrestart[0] = x[0];
    xrestart[1] = x[1];
    xrestart[2] = x[2];
    vrestart[0] = v[0];
    vrestart[1] = v[1];
    vrestart[2] = v[2];
    trestart = t;
    
    //insert restart loop here if necessary
    x[0] = xrestart[0];
    x[1] = xrestart[1];
    x[2] = xrestart[2];
    v[0] = vrestart[0];
    v[1] = vrestart[1];
    v[2] = vrestart[2];
    t = trestart;
    
    //integrate cluster orbit forwards from t = -tint till t = tstart+tfuture
    sign = 1.0;
    dtoutt = dtout;
    tmax = tfuture;
    if (radio) printf("\nt = %f\tdtout = %f\n",t,dtout);
    err = rk4_drv(&t,tmax,dtoutt,mdiff,&x[0],&v[0],sign,parameter,fit,newspaper);
    
	return 0;
}


/* --------------- extrapolation method --------------- */
int rk4_drv(double *t, double tmax, double dtout, double mdiff, double *x, double *v, double sign, double *parameter, double *fit, int newspaper){
	double tout,diff,dt, mdifft,dttemp = 0.0;
	double atemp[3], xe1[3], xe2[3], ve1[3], ve2[3], xt[3], vt[3], vmaxwell,vtemp, omega[3], omegat, xc[3], vc[3], ac[3];
	double actualclustermass;
    double rgalsun = parameter[6];
    double dsuntemp, vrsuntemp, vrtemp, ltemp, btemp, lcosbtemp, RAtemp, DECtemp, mu_alphatemp, mu_alphacosdeltatemp, mu_deltatemp, mutemp, PAtemp, vLSRtemp;
    vLSRtemp = parameter[7];
	int k,i;
	int err;
	double tt, r;
	mdifft = mdiff;
    char name[50], name_orbit[50];
    FILE *fz, *fz2;
    sprintf(name, "streaklinemodel.txt");
    sprintf(name_orbit, "clusterorbit.txt");
    double tpast;
    tpast = parameter[12];
    if (newspaper) fz = fopen(name,"w");
    if (newspaper) fz2 = fopen(name_orbit,"w");

	//initialize output/insertion of tail particles
	tout = *t;							/* time of next output/insertion */
	dt = sign*dt0;                /* initial time step */
    int NMAX = nperstep*2.0*(tstart-tpast)/(1.0*dtout);
    
	int columns = 5;
	double **star;
	star = (double **)calloc(NMAX,sizeof(double *));
	for (i=0;i<NMAX;i++){
		star[i] = (double *)calloc(columns,sizeof(double));
		if (star[i] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
    double rapo = 0.0;
    double rperi = 1.e10;

	int starcount = 0;
	

	//integrate cluster
	do {
		
		/***********
		 * CLUSTER *
		 ***********/
        
		if (sign**t>=sign*(tout)) {              /* output */
			
            //prepare coordinate transformation with convert, if output in galactic or equatorial coordinates
			//coordtype: 1 = equatorial to galactic and cartesian, 2 = galactic to equatorial and cartesian, 3 = cartesian to equatorial and galactic
			//vcoordtype: 1 = mu & position angle; 2 = mu_alpha & mu_delta or mu_alphacos(delta) & mu_delta; 3 = cartesian velocities

            //store coordinates and radial velocity of cluster particle
            convert(x, v, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSRtemp, rgalsun);
                
            if (ltemp>180.0) {
                ltemp -= 360.0;
                lcosbtemp = ltemp*cos(btemp/180.0*PI);
            }
            if (newspaper) fprintf(fz2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", x[0], x[1], x[2], v[0], v[1], v[2], ltemp, btemp, lcosbtemp, vrsuntemp, RAtemp, DECtemp);

            
			//integrate test tail particle for leading and trailing tail from current cluster position to time of observation, i.e. today
			if ((tails) && (sign>0.0) && (sign**t<=sign*tmax) && (sign**t<tstart))  {
				
				actualclustermass = parameter[0]+sqrt(*t**t)*parameter[8];
				parameter[9] = sqrt(*t**t);
				
				if (radio) printf("M(t) = %f\n", actualclustermass);
				

				/****************
				 * LEADING TAIL *
				 ****************/
				
				r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
				vtemp = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
				omega[0] = (x[1]*v[2]-x[2]*v[1])/(r*r);
				omega[1] = (x[2]*v[0]-x[0]*v[2])/(r*r);
				omega[2] = (x[0]*v[1]-x[1]*v[0])/(r*r);
				omegat = sqrt(omega[0]*omega[0]+omega[1]*omega[1]+omega[2]*omega[2]);
				
				double rtide, rtidef, dphi;
				
				double at[3], at2[3]; //force evaluation for dphi
				for (i=0;i<3;i++) {
					xt[i] = x[i]/(1.0*r)*(r-20.0);
					vt[i] = 0.0;
				}
				getforce(xt,vt,at,parameter);
				for (i=0;i<3;i++) xt[i] = x[i]/(1.0*r)*(r+20.0);
				getforce(xt,vt,at2,parameter);
                
				dphi = (sqrt(at[0]*at[0]+at[1]*at[1]+at[2]*at[2])-sqrt(at2[0]*at2[0]+at2[1]*at2[1]+at2[2]*at2[2]))/40.0;
				rtide = pow(G*parameter[0]/sqrt(pow((dphi+omegat*omegat),2)),1.0/3.0);
								
				if (rtide>parameter[5]) rtidef = rtide; //in case edge radius is wanted
				else rtidef = parameter[5];
                
                if (!radio) printf("T = %.0f\n", *t);
                
                
                //generate tail particles
                for (k=0;k<nperstep;k++) {
                    do {
                        tt = *t;
                        for (i=0;i<3;i++) {
                            xc[i] = x[i];
                            vc[i] = v[i];
                        }
        
                        if (radio) printf("tt = %.0f\trtide = %.1f\tomega = %.5f\tr = %.1f\t\t", tt, rtidef,omegat,r);
                        
                        vmaxwell = 1.0/3.0*sqrt(pow(get_gauss()*vrexcess,2)+pow(get_gauss()*vrexcess,2)+pow(get_gauss()*vrexcess,2));
                        
                        for (i=0;i<3;i++) {
                            xt[i] = x[i]/(1.0*r)*(r-rtidef)              +rexcess*rtidef*get_gauss();
                            vt[i] = v[i]/vtemp*(vtemp-omegat*rtidef)     -(vmaxwell*x[i]/(1.0*r))  + get_gauss()*sigma[i] ;
                        }
                        
                        err = rk4_drvtail(&tt,tstart,mdifft,&xt[0],&vt[0],sign,xc,vc,parameter);
                        
                        if (err==1) {
                            for (i=0;i<NMAX;i++) free (star[i]);
                            free(star);
                            if (newspaper) fclose(fz);
                            return 1;
                        } else if (err ==2) {
                            rtidef += 10.0;
                            if (rtidef>rtidemax) {
                                for (i=0;i<NMAX;i++) free (star[i]);
                                free(star);
                                return 1;
                            }
                        }
                    } while (err);
                    //store coordinates and radial velocity of tail particle
                    //coordtype: 1 = equatorial to galactic and cartesian, 2 = galactic to equatorial and cartesian, 3 = cartesian to equatorial and galactic
                    //vcoordtype: 1 = mu & position angle; 2 = mu_alpha & mu_delta or mu_alphacos(delta) & mu_delta; 3 = cartesian velocities
                    convert(xt, vt, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSRtemp,rgalsun);
                    
                    if (ltemp>180.0) {
                        ltemp -= 360.0;
                        lcosbtemp = ltemp*cos(btemp/180.0*PI);
                    }
                    if (newspaper) fprintf(fz,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", xt[0], xt[1], xt[2], vt[0], vt[1], vt[2], ltemp, btemp, lcosbtemp, vrsuntemp, RAtemp, DECtemp);
                    star[starcount][0] = ltemp;
                    star[starcount][1] = btemp;
                    star[starcount][3] = vrsuntemp;
                    starcount++;
//                    if (newspaper) fprintf(fz,"%.0lf %lf %lf %lf %lf %lf %lf\n", *t, xt[0], xt[1], xt[2], vt[0], vt[1], vt[2]);
                    
                }
                
				
                
                
				/*****************
				 * TRAILING TAIL *
				 *****************/
                
                for (k=0;k<nperstep;k++) {
                    do {
                        tt = *t;
                        for (i=0;i<3;i++) {
                            xc[i] = x[i];
                            vc[i] = v[i];
                        }
                        
                        vmaxwell = 1.0/3.0*sqrt(pow(get_gauss()*vrexcess,2)+pow(get_gauss()*vrexcess,2)+pow(get_gauss()*vrexcess,2));
                        
                        for (i=0;i<3;i++) {
                            xt[i] = x[i]/r*(r+rtidef)+rexcess*rtidef*get_gauss();
                            vt[i] = v[i]/vtemp*(vtemp+omegat*rtidef)+(vmaxwell*x[i]/(1.0*r))+ get_gauss()*sigma[i];
                        }
                        
                        err = rk4_drvtail(&tt,tstart,mdifft,&xt[0],&vt[0],sign,xc,vc,parameter);
                        
                        if (err==1) {
                            for (i=0;i<NMAX;i++) free (star[i]);
                            free(star);
                            if (newspaper) fclose(fz);
                            return 1;
                        } else if (err ==2) {
                            rtidef += 10.0;
                            if (rtidef>rtidemax) {
                                for (i=0;i<NMAX;i++) free (star[i]);
                                free(star);
                                return 1;
                            }
                        }
                    } while (err);
                    
                    //store coordinates and radial velocity of tail particle
                    //coordtype: 1 = equatorial to galactic and cartesian, 2 = galactic to equatorial and cartesian, 3 = cartesian to equatorial and galactic
                    //vcoordtype: 1 = mu & position angle; 2 = mu_alpha & mu_delta or mu_alphacos(delta) & mu_delta; 3 = cartesian velocities
                     convert(xt, vt, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSRtemp, rgalsun);
                    
                    if (ltemp>180.0) {
                        ltemp -= 360.0;
                        lcosbtemp = ltemp*cos(btemp/180.0*PI);
                    }
                    if (newspaper) fprintf(fz,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", xt[0], xt[1], xt[2], vt[0], vt[1], vt[2], ltemp, btemp, lcosbtemp, vrsuntemp, RAtemp, DECtemp);
                    star[starcount][0] = ltemp;
                    star[starcount][1] = btemp;
                    star[starcount][3] = vrsuntemp;
                    starcount++;
//                    if (newspaper) fprintf(fz,"%.0lf %lf %lf %lf %lf %lf %lf\n", *t, xt[0], xt[1], xt[2], vt[0], vt[1], vt[2]);
                }
            }
			
			tout+=dtout;              /* increase time of output/next insertion */
		}
		
		
        
        //advance cluster particle
		int count = 0;
        int laststep = 0;
		do {
			if (sign*(*t+dt) > sign*tout) {
				dt = tout-*t;
                laststep = 1;
			}
			for (k=0;k<3;k++) {
				xe1[k]=x[k];
				xe2[k]=x[k];
				ve1[k]=v[k];
				ve2[k]=v[k];
			}
			do_step(dt,xe1,ve1,parameter);      /* One full step */
			
			do_step(0.5*dt,xe2,ve2,parameter);  /* Two half steps */
			do_step(0.5*dt,xe2,ve2,parameter);
			
			diff = sqrt(pow(*xe1 - *xe2,2) + pow(*(xe1+1) - *(xe2+1),2) + pow(*(xe1+2) - *(xe2+2),2));
			
			if (diff<mdiff) {         /* Is difference below accuracy threshold? */
				*t+=dt;
				dttemp = dt;
				
				for (k=0;k<3;k++) {
					x[k]=xe2[k];          /* If yes -> continue and double step size */
					v[k]=ve2[k];
				}
				dt = dt*2.0;
			} else {
				dt = dt/2.0;
			}
			
			if (sign*dt < 0.01*dt0 && !laststep) {
				printf("Aborted... dt = %lf (>%lf, %lf)\n", dt, dt0, sign);
				return 1;
			}
			count++;
			
			if (sign*dt > dtmax) {
				dt = sign*dtmax;
			}
		} while (diff>mdiff);       /* Go through loop once and only repeat if difference is too large */
		if ((sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) < Rgalmin) || (sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) > Rgalmax)  ) return 1;  /* Abort if r is too extreme */
        if (sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) > rapo) rapo = sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));
        if (sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) < rperi) rperi = sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));
        
    } while (sign**t<sign*(tmax));

    
    //final output of cluster orbit
    if (sign>0.0) {
        //store coordinates and radial velocity of cluster particle
        //coordtype: 1 = equatorial to galactic and cartesian, 2 = galactic to equatorial and cartesian, 3 = cartesian to equatorial and galactic
        //vcoordtype: 1 = mu & position angle; 2 = mu_alpha & mu_delta or mu_alphacos(delta) & mu_delta; 3 = cartesian velocities
        convert(x, v, &dsuntemp, &vrsuntemp, &vrtemp, &ltemp, &btemp, &lcosbtemp, &RAtemp, &DECtemp, &mu_alphatemp, &mu_alphacosdeltatemp, &mu_deltatemp, &mutemp, &PAtemp, 3, 3, 0, vLSRtemp, rgalsun);
        
        if (ltemp>180.0) {
            ltemp -= 360.0;
            lcosbtemp = ltemp*cos(btemp/180.0*PI);
        }
        if (newspaper) fprintf(fz2,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", x[0], x[1], x[2], v[0], v[1], v[2], ltemp, btemp, lcosbtemp, vrsuntemp, RAtemp, DECtemp);
    }

    
    if (newspaper) fclose(fz);
    if (newspaper) fclose(fz2);
    
	for (i=0;i<NMAX;i++) free (star[i]);
	free(star);
	
	return 0;
	
}

/* --------------- extrapolation method tail --------------- */
int rk4_drvtail(double *t, double tmax, double mdiff, double *x, double *v, double sign, double *xc, double *vc, double *parameter){
	double diff,dt;
	double xe1[3], xe2[3], ve1[3], ve2[3];
	double xce1[3], xce2[3], vce1[3],vce2[3];
	int k;
    int laststep = 0;
    
	dt = sign*dt0;		/* initial time step */
    
	while (sign**t<sign*tmax) {
		if (sign*(*t+dt) > sign*tmax) {
			dt = tmax-*t;
            laststep = 1;
		}
		
		do {
			for (k=0;k<3;k++) {
				xe1[k]=x[k];
				xe2[k]=x[k];
				ve1[k]=v[k];
				ve2[k]=v[k];
				xce1[k]=xc[k];
				xce2[k]=xc[k];
				vce1[k]=vc[k];
				vce2[k]=vc[k];
			}
			
			parameter[9] = sqrt(*t**t);
			
			do_steptail(dt,xe1,ve1,xce1,vce1, parameter);
			
			do_steptail(0.5*dt,xe2,ve2,xce2,vce2, parameter); 
			do_steptail(0.5*dt,xe2,ve2,xce2,vce2, parameter);
			
			diff = sqrt(pow(*xe1 - *xe2,2) + pow(*(xe1+1) - *(xe2+1),2) + pow(*(xe1+2) - *(xe2+2),2));
			
			if (diff<mdiff) {       
				*t+=dt;
				
				for (k=0;k<3;k++) {
					x[k]=xe2[k];        
					v[k]=ve2[k];
					xc[k]=xce2[k];
					vc[k]=vce2[k];
				}
				dt = dt*2.0;
                
			} else {
				dt = dt/2.0;
			}
			
			if ((sign*dt < 0.01*dt0) && !laststep) {
				printf("Aborted... dt = %lf\n", dt);
				*t = tmax*2.0;
				diff = mdiff/2;
				return 1;
			}
			
		} while (diff>mdiff);     
		if (sqrt(pow(x[0]-xc[0],2)+pow(x[1]-xc[1],2)+pow(x[2]-xc[2],2)) < Rstop) return 2; //increase Rtide
		else if (sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2)) < Rgalmin) return 1;
	}
	
	return 0;
	
}

/* ----------- force ----------- */
void getforce(double *x, double *v, double *a, double *parameter){
	double r1, r2, r3;
	double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z;
	
	if (gpot == 1) {
		//Allen & Santillan (1991) potential w updated values from Irrgang et al. (2013)
		
		//Point mass
		r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2) + b1*b1);
		
		a1x = -G*M1/(r1*r1*r1)**x;
		a1y = -G*M1/(r1*r1*r1)**(x+1);
		a1z = -G*M1/(r1*r1*r1)**(x+2);
		
		//Miyamato disk
		r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2 + sqrt(*(x+2) * *(x+2) + b2*b2),2));
		
		a2x = -G*M2/(r2*r2*r2) * *x;
		a2y = -G*M2/(r2*r2*r2) * *(x+1);
		a2z = -G*M2/(r2*r2*r2) * (a2 + sqrt(*(x+2) * *(x+2) + b2*b2))/sqrt(*(x+2) * *(x+2) + b2*b2) * *(x+2);
        
		//Log Halo
		r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]);
        
		a3x = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *x/r3;
		a3y = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *(x+1)/r3;
		a3z = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *(x+2)/(parameter[10]*parameter[10]*r3);
        
	} else if (gpot == 2) {
		//Log Halo from Koposov et al. (2010)
		r3 = *x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]; //R^2!
        
		a3x = -G*parameter[3]/parameter[11] *  *x/r3;
		a3y = -G*parameter[3]/parameter[11] *  *(x+1)/r3;
		a3z = -G*parameter[3]/parameter[11] *  *(x+2)/(parameter[10]*parameter[10]*r3);
        
		a1x = a1y = a1z = a2x = a2y = a2z = 0.0;
		
	} else if (gpot == 3) {
		//potential from Johnston/Law/Majewski/Helmi
		
		//Jaffe bulge
		r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));
		
		a1x = -G*M1_LMJ/((r1+b1_LMJ)*r1)**x/r1;
		a1y = -G*M1_LMJ/((r1+b1_LMJ)*r1)**(x+1)/r1;
		a1z = -G*M1_LMJ/((r1+b1_LMJ)*r1)**(x+2)/r1;
		
		//Miyamato disk
		r2 = sqrt(*x * *x + *(x+1) * *(x+1) + (a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ))*(a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ)));
		
		a2x = -G*M2_LMJ/(r2*r2*r2) * *x;
		a2y = -G*M2_LMJ/(r2*r2*r2) * *(x+1);
		a2z = -G*M2_LMJ/(r2*r2*r2) * (a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ))/sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ) * *(x+2);
		
		//NFW Halo
		r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]);
		
		a3x = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *x/r3;
		a3y = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+1)/r3;
		a3z = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+2)/(parameter[10]*parameter[10]*r3);
        
    } else if (gpot == 4) {
        //potential from Law&Majewski'10
        double fac, C1, C2, C3;
        
        //Hernquist bulge
        r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));
        
        a1x = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**x/r1;
        a1y = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+1)/r1;
        a1z = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+2)/r1;
        
        //Miyamato disk
        r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ),2));
        
        a2x = -G*M2_LMJ/(r2*r2*r2) * *x;
        a2y = -G*M2_LMJ/(r2*r2*r2) * *(x+1);
        a2z = -G*M2_LMJ/(r2*r2*r2) * (a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ))/sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ) * *(x+2);
        
        //Log Halo
        C1 = pow(1.0*cos(phi_lm)/q1_lm,2)+pow(1.0*sin(phi_lm)/q2_lm,2);
        C2 = pow(1.0*cos(phi_lm)/q2_lm,2)+pow(1.0*sin(phi_lm)/q1_lm,2);
        C3 = 2.*sin(phi_lm)*cos(phi_lm)*(1./q1_lm/q1_lm-1./q2_lm/q2_lm);
        
        r3 = C1* *x * *x + C2* *(x+1) * *(x+1) + C3* *x * *(x+1) + *(x+2)/qz_lm * *(x+2)/qz_lm;
        
        fac = 0.5*vh2_lm/(r3+c_lm*c_lm);
        a3x = -fac*(2.*C1* *x +C3 * *(x+1));
        a3y = -fac*(2.*C2* *(x+1)+C3* *x);
        a3z = -fac*(2.* *(x+2)/qz_lm/qz_lm);
        
    }

    
	*(a+0) = a1x + a2x + a3x;
	*(a+1) = a1y + a2y + a3y;
	*(a+2) = a1z + a2z + a3z;
    
}

/* ----------- force tail ----------- */
void getforcetail(double *x, double *v, double *a, double *xc, double *parameter){
	double r1, r2, r3, r4;
	double a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, a4x, a4y, a4z;
	double actualmass;
	actualmass = parameter[0]+parameter[9]*parameter[8];
	
	if (gpot == 1) {
		//Allen & Santillan (1991) potential w updated values from Irrgang et al. (2013)
		
		//Point mass
		r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2) + b1*b1);
		
		a1x = -G*M1/(r1*r1*r1)**x;
		a1y = -G*M1/(r1*r1*r1)**(x+1);
		a1z = -G*M1/(r1*r1*r1)**(x+2);
		
		//Miyamato disk
		r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2 + sqrt(*(x+2) * *(x+2) + b2*b2),2));
		
		a2x = -G*M2/(r2*r2*r2) * *x;
		a2y = -G*M2/(r2*r2*r2) * *(x+1);
		a2z = -G*M2/(r2*r2*r2) * (a2 + sqrt(*(x+2) * *(x+2) + b2*b2))/sqrt(*(x+2) * *(x+2) + b2*b2) * *(x+2);
		
		//Log Halo
		r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]);
		
		a3x = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *x/r3;
		a3y = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *(x+1)/r3;
		a3z = -G*parameter[3]/(parameter[11]*parameter[11] +parameter[11]*r3) * *(x+2)/(parameter[10]*parameter[10]*r3);
		
	} else if (gpot == 2) {
		//Log Halo from Koposov et al. (2010)
		r3 = *x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]; //R^2!
		
		a3x = -G*parameter[3]/parameter[11] *  *x/r3;
		a3y = -G*parameter[3]/parameter[11] *  *(x+1)/r3;
		a3z = -G*parameter[3]/parameter[11] *  *(x+2)/(parameter[10]*parameter[10]*r3);
		
		a1x = a1y = a1z = a2x = a2y = a2z = 0.0;
		
	} else if (gpot == 3) {
		//potential from Johnston/Law/Majewski/Helmi
		
		//Jaffe bulge
		r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));
		
		a1x = -G*M1_LMJ/((r1+b1_LMJ)*r1)**x/r1;
		a1y = -G*M1_LMJ/((r1+b1_LMJ)*r1)**(x+1)/r1;
		a1z = -G*M1_LMJ/((r1+b1_LMJ)*r1)**(x+2)/r1;
		
		//Miyamato disk
		r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ),2));
		
		a2x = -G*M2_LMJ/(r2*r2*r2) * *x;
		a2y = -G*M2_LMJ/(r2*r2*r2) * *(x+1);
		a2z = -G*M2_LMJ/(r2*r2*r2) * (a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ))/sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ) * *(x+2);
		
		//NFW Halo
		r3 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2)/parameter[10] * *(x+2)/parameter[10]);
		
		a3x = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *x/r3;
		a3y = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+1)/r3;
		a3z = -G*parameter[3]/r3 * (log(1.0 + r3/parameter[11])/r3 - 1.0/(parameter[11]+r3)) * *(x+2)/(parameter[10]*parameter[10]*r3);
		
    } else if (gpot == 4) {
        //potential from Law&Majewski'10
        double fac, C1, C2, C3;

        //Hernquist bulge
        r1 = sqrt(*x * *x + *(x+1) * *(x+1) + *(x+2) * *(x+2));
        
        a1x = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**x/r1;
        a1y = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+1)/r1;
        a1z = -G*M1_LMJ/((r1+b1_LMJ)*(r1+b1_LMJ))**(x+2)/r1;
        
        //Miyamato disk
        r2 = sqrt(*x * *x + *(x+1) * *(x+1) + pow(a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ),2));
        
        a2x = -G*M2_LMJ/(r2*r2*r2) * *x;
        a2y = -G*M2_LMJ/(r2*r2*r2) * *(x+1);
        a2z = -G*M2_LMJ/(r2*r2*r2) * (a2_LMJ + sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ))/sqrt(*(x+2) * *(x+2) + b2_LMJ*b2_LMJ) * *(x+2);
        
        //Log Halo
        C1 = pow(1.0*cos(phi_lm)/q1_lm,2)+pow(1.0*sin(phi_lm)/q2_lm,2);
        C2 = pow(1.0*cos(phi_lm)/q2_lm,2)+pow(1.0*sin(phi_lm)/q1_lm,2);
        C3 = 2.*sin(phi_lm)*cos(phi_lm)*(1./q1_lm/q1_lm-1./q2_lm/q2_lm);
        
        r3 = C1* *x * *x + C2* *(x+1) * *(x+1) + C3* *x * *(x+1) + *(x+2)/qz_lm * *(x+2)/qz_lm;
        
        fac = 0.5*vh2_lm/(r3+c_lm*c_lm);
        a3x = -fac*(2.*C1* *x +C3 * *(x+1));
        a3y = -fac*(2.*C2* *(x+1)+C3* *x);
        a3z = -fac*(2.* *(x+2)/qz_lm/qz_lm);
        
    }
	
	//Cluster
	r4 = sqrt((*x-*xc) * (*x-*xc) + (*(x+1)-*(xc+1)) * (*(x+1)-*(xc+1)) + (*(x+2)-*(xc+2)) * (*(x+2)-*(xc+2)) + R4*R4);
	
	a4x = -G*actualmass/(r4*r4*r4)*(*x-*xc);
	a4y = -G*actualmass/(r4*r4*r4)*(*(x+1)-*(xc+1));
	a4z = -G*actualmass/(r4*r4*r4)*(*(x+2)-*(xc+2));
	
	*(a+0) = a1x + a2x + a3x + a4x;
	*(a+1) = a1y + a2y + a3y + a4y;
	*(a+2) = a1z + a2z + a3z + a4z;
}

/* ---------- advancement ---------- */
void do_step(double dt, double *x, double *v, double *parameter) {
	double hh, acc0[3], acc1[3], acc2[3], acc3[3],xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	int k;
    
	hh = dt*0.5;
	
	getforce(x,v,acc0, parameter);
	for (k=0;k<3;k++) {                /* first half-step */
		xt1[k] = *(x+k)+hh**(v+k);
		vt1[k] = *(v+k)+hh**(acc0+k);
	}
	
	getforce(&xt1[0], &vt1[0], acc1, parameter);
	for (k=0;k<3;k++) {                /* second half-step */
		xt2[k] = *(x+k)+hh*vt1[k];
		vt2[k] = *(v+k)+hh**(acc1+k);
	}
	
	getforce(&xt2[0], &vt2[0], acc2, parameter);
	for (k=0;k<3;k++) {                /* thrid half-step with results of second half-step */
		xt3[k] = *(x+k)+dt*vt2[k];
		vt3[k] = *(v+k)+dt**(acc2+k);
	}
	
	getforce(&xt3[0], &vt3[0], acc3, parameter);
	for (k=0;k<3;k++) {                /* Runge-Kutta formula */
		*(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
		*(v+k) += dt/6.0*(*(acc0+k)+2.0*(*(acc1+k)+*(acc2+k))+*(acc3+k));
	}
    
}

/* ---------- advancement tail ---------- */
void do_steptail(double dt, double *x, double *v, double *xc, double *vc, double *parameter) {
	double hh, acc0[3], acc1[3], acc2[3], acc3[3] ,xt1[3],xt2[3],xt3[3],vt1[3],vt2[3],vt3[3];
	double ac0[3],ac1[3],ac2[3],ac3[3],xct1[3],xct2[3],xct3[3],vct1[3],vct2[3],vct3[3];
	int k;
	
	hh = dt*0.5;
	
	getforcetail(x,v,acc0,xc,parameter);
	getforce(xc,vc,ac0,parameter);
	for (k=0;k<3;k++) {
		xt1[k] = *(x+k)+hh**(v+k);
		vt1[k] = *(v+k)+hh**(acc0+k);
		xct1[k] = *(xc+k)+hh**(vc+k);
		vct1[k] = *(vc+k)+hh**(ac0+k);
	}
	
	getforcetail(&xt1[0], &vt1[0], acc1, xct1, parameter);
	getforce(&xct1[0], &vct1[0], ac1, parameter);
	for (k=0;k<3;k++) {
		xt2[k] = *(x+k)+hh*vt1[k];
		vt2[k] = *(v+k)+hh**(acc1+k);
		xct2[k] = *(xc+k)+hh*vct1[k];
		vct2[k] = *(vc+k)+hh**(ac1+k);
	}
	
	getforcetail(&xt2[0], &vt2[0], acc2, xct2, parameter);
	getforce(&xct2[0], &vct2[0], ac2, parameter);
	for (k=0;k<3;k++) {            
		xt3[k] = *(x+k)+dt*vt2[k];
		vt3[k] = *(v+k)+dt**(acc2+k);
		xct3[k] = *(xc+k)+dt*vct2[k];
		vct3[k] = *(vc+k)+dt**(ac2+k);
	}
	
	getforcetail(&xt3[0], &vt3[0], acc3, xct3, parameter);
	getforce(&xct3[0], &vct3[0], ac3, parameter);
	for (k=0;k<3;k++) {               
		*(x+k) += dt/6.0*(*(v+k)+2.0*(vt1[k]+vt2[k])+vt3[k]);
		*(v+k) += dt/6.0*(*(acc0+k)+2.0*(*(acc1+k)+*(acc2+k))+*(acc3+k));
		*(xc+k) += dt/6.0*(*(vc+k)+2.0*(vct1[k]+vct2[k])+vct3[k]);
		*(vc+k) += dt/6.0*(*(ac0+k)+2.0*(*(ac1+k)+*(ac2+k))+*(ac3+k));
	}
    
}

/* ---------- gaussian distribution ---------- */
double get_gauss(void){
	double random[2],p,q;
	do {
		random[0] = 2.0*drand48()-1.0;
		random[1] = 2.0*drand48()-1.0;
		q = random[0]*random[0]+random[1]*random[1];
	} while (q>1.0);
	
	p = sqrt(-2.0*log(q)/q);
	return random[0]*p;
	
}

/* ---------- coordinate conversion ---------- */
void convert(double *xtemp, double *vtemp, double *dsuntemp, double *vrsuntemp, double *vrtemp, double *ltemp, double *btemp, double *lcosbtemp, double *RAtemp, double *DECtemp, double *mu_alphatemp, double *mu_alphacosdeltatemp, double *mu_deltatemp, double *mutemp, double *PAtemp, int coordtype_coco, int vcoordtype, int radiococo, double vLSRtemp, double rgalsun){
	
	double x,y,z;        //galactic coordinates [kpc]
	double xsun = -rgalsun/1000.0;//galactocentric distance of sun [kpc]
	double dsun_coco;			//heliocentric radial distance [kpc]
	double dxy;             //heliocentric distance in xy-plane [kpc]
	double vx,vy,vz;        //cluster velocity [km/s]
	double vrx,vry,vrz;     //cluster radial velocity in 3d coordinates [km/s]
	double vtx,vty,vtz;     //cluster tangential velocity in 3d coordinates [km/s]
	double T[3][3], A[3][3], B[3][3];
	double TI[3][3];
	double RArad, DECrad;
	double brad, lrad, lcosbrad;
	double mu_alphacosdelta_coco, mu_delta_coco, mu, PArad;
	double vr_coco, vrsun;
	double RAENP = 0.0, DECENP = PI/2.0;  //equatorial coordinates of equatorial north pole
	double xENP, yENP, zENP, dxyENP; //cartesian vector pointing to the equatorial north pole
	double bENP, lENP; //galactic coordinates of ENP
	double FAK;
	double xdelta, ydelta, zdelta;
	double nx, ny, nz, dvt;
	double vrLSR, vrGSR;
	
	
	//transformation matrix equ. -> gal. from Johnson & Soderblom (1987)
/*	double t,d,a;
    double detT;
	t = PAGNP/360.0*2.0*PI;
	d = deltaGNP/360.0*2.0*PI;
	a = alphaGNP/360.0*2.0*PI;
	
	T[0][0] = -cos(t)*sin(d)*cos(a)-sin(t)*sin(a);
	T[0][1] = -cos(t)*sin(d)*sin(a)+sin(t)*cos(a);
	T[0][2] = cos(t)*cos(d);
	
	T[1][0] = -sin(t)*sin(d)*cos(a)+cos(t)*sin(a);
	T[1][1] = -sin(t)*sin(d)*sin(a)-cos(t)*cos(a);
	T[1][2] = sin(t)*cos(d);
	
	T[2][0] = cos(d)*cos(a);
	T[2][1] = cos(d)*sin(a);
	T[2][2] = sin(d);

	//invert matrix T in the most general way
	detT = T[0][0]*T[1][1]*T[2][2] + T[1][0]*T[2][1]*T[0][2] + T[2][0]*T[0][1]*T[1][2] - T[0][0]*T[1][2]*T[2][1] - T[1][0]*T[2][2]*T[0][1] - T[2][0]*T[0][2]*T[1][1];

	TI[0][0] = (T[1][1]*T[2][2]-T[1][2]*T[2][1])/detT;
	TI[1][0] = (T[1][2]*T[2][0]-T[1][0]*T[2][2])/detT;
	TI[2][0] = (T[1][0]*T[2][1]-T[1][1]*T[2][0])/detT;
	
	TI[0][1] = (T[0][2]*T[2][1]-T[0][1]*T[2][2])/detT;
	TI[1][1] = (T[0][0]*T[2][2]-T[2][0]*T[0][2])/detT;
	TI[2][1] = (T[0][1]*T[2][0]-T[0][0]*T[2][1])/detT;
	
	TI[0][2] = (T[0][1]*T[1][2]-T[0][2]*T[1][1])/detT;
	TI[1][2] = (T[0][2]*T[1][0]-T[0][0]*T[1][2])/detT;
	TI[2][2] = (T[0][0]*T[1][1]-T[0][1]*T[1][0])/detT;
*/

    //use result right away, careful when changing epochs though
    T[0][0] = -0.0548765333890783;
    T[0][1] = -0.8734366584039325;
    T[0][2] = -0.4838356847519307;
    T[1][0] = 0.4941106597543148;
    T[1][1] = -0.4448303583524283;
    T[1][2] = 0.7469809958795511;
    T[2][0] = -0.8676653859641706;
    T[2][1] = -0.1980766418440660;
    T[2][2] = 0.4559851115501738;

    TI[0][0] = -0.0548765333890783;
    TI[0][1] = 0.4941106597543147;
    TI[0][2] = -0.8676653859641704;
    TI[1][0] = -0.8734366584039324;
    TI[1][1] = -0.4448303583524283;
    TI[1][2] = -0.1980766418440659;
    TI[2][0] = -0.4838356847519305;
    TI[2][1] = 0.7469809958795510;
    TI[2][2] = 0.4559851115501737;
    
    
	//convert to kpc
	x = xtemp[0]/1000.0;
	y = xtemp[1]/1000.0;
	z = xtemp[2]/1000.0;
	
	dsun_coco = *dsuntemp/1000.0;
	
	vx = vtemp[0];
	vy = vtemp[1];
	vz = vtemp[2];
	
	vr_coco = *vrtemp;
	vrsun = *vrsuntemp;
	
	//convert to radians
	DECrad = *DECtemp/360.0*2.0*PI;
	RArad = *RAtemp/360.0*2.0*PI;
	PArad = *PAtemp/360.0*2.0*PI;
	
	
	//get the galactic coordinates first
	if (coordtype_coco == 1) {
		if (radiococo) printf("\nConverting equatorial to galactic coordinates using the transformation matrix:\n");
		if (radiococo) printf("%f\t%f\t%f\n",T[0][0],T[0][1],T[0][2]);
		if (radiococo) printf("%f\t%f\t%f\n",T[1][0],T[1][1],T[1][2]);
		if (radiococo) printf("%f\t%f\t%f\n",T[2][0],T[2][1],T[2][2]);
		
		brad = asin(T[2][0]*cos(DECrad)*cos(RArad) + T[2][1]*cos(DECrad)*sin(RArad) + T[2][2]*sin(DECrad));
		if (asin((T[1][0]*cos(DECrad)*cos(RArad) + T[1][1]*cos(DECrad)*sin(RArad) + T[1][2]*sin(DECrad))/cos(brad))>=0.0) {
			lrad = acos((T[0][0]*cos(DECrad)*cos(RArad) + T[0][1]*cos(DECrad)*sin(RArad) + T[0][2]*sin(DECrad))/cos(brad));
		} else {
			lrad = 2.0*PI-acos((T[0][0]*cos(DECrad)*cos(RArad) + T[0][1]*cos(DECrad)*sin(RArad) + T[0][2]*sin(DECrad))/cos(brad));
		}
		lcosbrad = lrad*cos(brad);
	} else if (coordtype_coco == 2) {
		brad = *btemp/360.0*2.0*PI;
		if (*ltemp) {
			lrad = *ltemp/360.0*2.0*PI;
			lcosbrad = lrad*cos(brad);
		} else if (*lcosbtemp) {
			lcosbrad = *lcosbtemp/360.0*2.0*PI;
			lrad = lcosbrad/cos(brad);
		}
	} else if (coordtype_coco == 3) {
		if (y >= 0.0)
			lrad = acos((x-xsun)/sqrt(pow(x-xsun,2)+y*y));
		else
			lrad = 2.0*Pi-acos((x-xsun)/sqrt(pow(x-xsun,2)+y*y));
		brad =  atan(z/sqrt(pow(x-xsun,2)+y*y));
		lcosbrad = lrad*cos(brad);
	}
	
	
	//get 3d position of cluster [kpc] from galactic coordinates
	if (coordtype_coco < 3) {
		z = sin(brad)*dsun_coco;
		dxy = sqrt(dsun_coco*dsun_coco-z*z);
		x = cos(lrad)*dxy + xsun;
		y = sin(lrad)*dxy;
	} else if (coordtype_coco == 3) {
		dsun_coco = sqrt(pow(x-xsun,2)+y*y+z*z);
	}
	
	
	//finally get the equatorial coordinates from galactic coordinates
	if (coordtype_coco > 1) {
		if (radiococo) printf("\nConverting galactic to equatorial coordinates using the transformation matrix:\n");
		if (radiococo) printf("%f\t%f\t%f\n",TI[0][0],TI[0][1],TI[0][2]);
		if (radiococo) printf("%f\t%f\t%f\n",TI[1][0],TI[1][1],TI[1][2]);
		if (radiococo) printf("%f\t%f\t%f\n",TI[2][0],TI[2][1],TI[2][2]);
		
		if (radiococo) {
			//unit matrix B = T * TI
			B[0][0] = T[0][0]*TI[0][0] + T[0][1]*TI[1][0] + T[0][2]*TI[2][0];
			B[0][1] = T[0][0]*TI[0][1] + T[0][1]*TI[1][1] + T[0][2]*TI[2][1];
			B[0][2] = T[0][0]*TI[0][2] + T[0][1]*TI[1][2] + T[0][2]*TI[2][2];
			
			B[1][0] = T[1][0]*TI[0][0] + T[1][1]*TI[1][0] + T[1][2]*TI[2][0];
			B[1][1] = T[1][0]*TI[0][1] + T[1][1]*TI[1][1] + T[1][2]*TI[2][1];
			B[1][2] = T[1][0]*TI[0][2] + T[1][1]*TI[1][2] + T[1][2]*TI[2][2];
			
			B[2][0] = T[2][0]*TI[0][0] + T[2][1]*TI[1][0] + T[2][2]*TI[2][0];
			B[2][1] = T[2][0]*TI[0][1] + T[2][1]*TI[1][1] + T[2][2]*TI[2][1];
			B[2][2] = T[2][0]*TI[0][2] + T[2][1]*TI[1][2] + T[2][2]*TI[2][2];
			
			printf("\nCalculating T*T^{-1} = 1 for consistency check:\n");
			printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
			printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
			printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		}
		
		DECrad = asin(TI[2][0]*cos(brad)*cos(lrad)+TI[2][1]*cos(brad)*sin(lrad)+TI[2][2]*sin(brad));
		if (asin((TI[1][0]*cos(brad)*cos(lrad) + TI[1][1]*cos(brad)*sin(lrad) + TI[1][2]*sin(brad))/cos(DECrad))>=0.0) {
			RArad = acos((TI[0][0]*cos(brad)*cos(lrad) + TI[0][1]*cos(brad)*sin(lrad) + TI[0][2]*sin(brad))/cos(DECrad));
		} else {
			RArad = 2.0*PI-acos((TI[0][0]*cos(brad)*cos(lrad) + TI[0][1]*cos(brad)*sin(lrad) + TI[0][2]*sin(brad))/cos(DECrad));
		}
	}
	
	
	
	//get tangential velocity in [km/s] from different sets of velocity-measurement types
	
	//get coordinates of equatorial north pole on great circle
	bENP = asin(T[2][0]*cos(DECENP)*cos(RAENP) + T[2][1]*cos(DECENP)*sin(RAENP) + T[2][2]*sin(DECENP));
	if (asin((T[1][0]*cos(DECENP)*cos(RAENP) + T[1][1]*cos(DECENP)*sin(RAENP) + T[1][2]*sin(DECENP))/cos(bENP))>=0.0) {
		lENP = acos((T[0][0]*cos(DECENP)*cos(RAENP) + T[0][1]*cos(DECENP)*sin(RAENP) + T[0][2]*sin(DECENP))/cos(bENP));
	} else {
		lENP = 2.0*PI-acos((T[0][0]*cos(DECENP)*cos(RAENP) + T[0][1]*cos(DECENP)*sin(RAENP) + T[0][2]*sin(DECENP))/cos(bENP));
	}
	if (radiococo) printf("\nCoordinates of equatorial north pole:\n");
	if (radiococo) printf("bENP = %f\tlENP = %f\n", bENP, lENP);
	zENP = sin(bENP)*dsun_coco;
	dxyENP = sqrt(dsun_coco*dsun_coco-zENP*zENP);
	xENP = cos(lENP)*dxyENP + xsun;
	yENP = sin(lENP)*dxyENP;
	if (radiococo) printf("xENP = %f\tyENP = %f\tzENP = %f\n", xENP, yENP, zENP);
	
	
	if (vcoordtype == 1) {
		
		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun_coco*vrsun;
		vry = y/dsun_coco*vrsun;
		vrz = z/dsun_coco*vrsun;
		if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));
		
		//convert to km/s
		mu = *mutemp*dsun_coco*4.74057;
		
		//compute proper motion components
		mu_alphacosdelta_coco = mu*sin(PArad);
		mu_delta_coco = mu*cos(PArad);
		
		A[0][0] = cos(RArad)*cos(DECrad);
		A[0][1] = -sin(RArad);
		A[0][2] = -cos(RArad)*sin(DECrad);
		
		A[1][0] = sin(RArad)*cos(DECrad);
		A[1][1] = cos(RArad);
		A[1][2] = -sin(RArad)*sin(DECrad);
		
		A[2][0] = sin(DECrad);
		A[2][1] = 0.0;
		A[2][2] = cos(DECrad);
		
		//printf("%f\t%f\t%f\n",A[0][0],A[0][1],A[0][2]);
		//printf("%f\t%f\t%f\n",A[1][0],A[1][1],A[1][2]);
		//printf("%f\t%f\t%f\n",A[2][0],A[2][1],A[2][2]);
		
		//B = T * A
		B[0][0] = T[0][0]*A[0][0] + T[0][1]*A[1][0] + T[0][2]*A[2][0];
		B[0][1] = T[0][0]*A[0][1] + T[0][1]*A[1][1] + T[0][2]*A[2][1];
		B[0][2] = T[0][0]*A[0][2] + T[0][1]*A[1][2] + T[0][2]*A[2][2];
		
		B[1][0] = T[1][0]*A[0][0] + T[1][1]*A[1][0] + T[1][2]*A[2][0];
		B[1][1] = T[1][0]*A[0][1] + T[1][1]*A[1][1] + T[1][2]*A[2][1];
		B[1][2] = T[1][0]*A[0][2] + T[1][1]*A[1][2] + T[1][2]*A[2][2];
		
		B[2][0] = T[2][0]*A[0][0] + T[2][1]*A[1][0] + T[2][2]*A[2][0];
		B[2][1] = T[2][0]*A[0][1] + T[2][1]*A[1][1] + T[2][2]*A[2][1];
		B[2][2] = T[2][0]*A[0][2] + T[2][1]*A[1][2] + T[2][2]*A[2][2];
		
		//printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
		//printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
		//printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		
		vx = vrsun*B[0][0] + mu_alphacosdelta_coco*B[0][1] + mu_delta_coco*B[0][2] +vxsun;
		vy = vrsun*B[1][0] + mu_alphacosdelta_coco*B[1][1] + mu_delta_coco*B[1][2] +vysun+vLSRtemp;
		vz = vrsun*B[2][0] + mu_alphacosdelta_coco*B[2][1] + mu_delta_coco*B[2][2] +vzsun;
		
		
		if (radiococo) printf("\nCartesian velocity:\n");
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));
		
	} else if (vcoordtype == 2) {
		
		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun_coco*vrsun;
		vry = y/dsun_coco*vrsun;
		vrz = z/dsun_coco*vrsun;
		if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));
		
		if (*mu_alphatemp) *mu_alphacosdeltatemp = *mu_alphatemp*cos(DECrad);
		else if (*mu_alphacosdeltatemp) *mu_alphatemp = *mu_alphacosdeltatemp/cos(DECrad);
		
		//convert to km/s
		mu_alphacosdelta_coco = *mu_alphacosdeltatemp*dsun_coco*4.74057;
		mu_delta_coco = *mu_deltatemp*dsun_coco*4.74057;
		mu = sqrt(mu_alphacosdelta_coco*mu_alphacosdelta_coco+mu_delta_coco*mu_delta_coco);
		
		A[0][0] = cos(RArad)*cos(DECrad);
		A[0][1] = -sin(RArad);
		A[0][2] = -cos(RArad)*sin(DECrad);
		
		A[1][0] = sin(RArad)*cos(DECrad);
		A[1][1] = cos(RArad);
		A[1][2] = -sin(RArad)*sin(DECrad);
		
		A[2][0] = sin(DECrad);
		A[2][1] = 0.0;
		A[2][2] = cos(DECrad);
		
		//printf("%f\t%f\t%f\n",A[0][0],A[0][1],A[0][2]);
		//printf("%f\t%f\t%f\n",A[1][0],A[1][1],A[1][2]);
		//printf("%f\t%f\t%f\n",A[2][0],A[2][1],A[2][2]);
		
		//B = T * A
		B[0][0] = T[0][0]*A[0][0] + T[0][1]*A[1][0] + T[0][2]*A[2][0];
		B[0][1] = T[0][0]*A[0][1] + T[0][1]*A[1][1] + T[0][2]*A[2][1];
		B[0][2] = T[0][0]*A[0][2] + T[0][1]*A[1][2] + T[0][2]*A[2][2];
		
		B[1][0] = T[1][0]*A[0][0] + T[1][1]*A[1][0] + T[1][2]*A[2][0];
		B[1][1] = T[1][0]*A[0][1] + T[1][1]*A[1][1] + T[1][2]*A[2][1];
		B[1][2] = T[1][0]*A[0][2] + T[1][1]*A[1][2] + T[1][2]*A[2][2];
		
		B[2][0] = T[2][0]*A[0][0] + T[2][1]*A[1][0] + T[2][2]*A[2][0];
		B[2][1] = T[2][0]*A[0][1] + T[2][1]*A[1][1] + T[2][2]*A[2][1];
		B[2][2] = T[2][0]*A[0][2] + T[2][1]*A[1][2] + T[2][2]*A[2][2];
		
		//printf("%f\t%f\t%f\n",B[0][0],B[0][1],B[0][2]);
		//printf("%f\t%f\t%f\n",B[1][0],B[1][1],B[1][2]);
		//printf("%f\t%f\t%f\n",B[2][0],B[2][1],B[2][2]);
		
		vx = vrsun*B[0][0] + mu_alphacosdelta_coco*B[0][1] + mu_delta_coco*B[0][2] +vxsun;
		vy = vrsun*B[1][0] + mu_alphacosdelta_coco*B[1][1] + mu_delta_coco*B[1][2] +vysun+vLSRtemp;
		vz = vrsun*B[2][0] + mu_alphacosdelta_coco*B[2][1] + mu_delta_coco*B[2][2] +vzsun;
		
		if (radiococo) printf("\nCartesian velocity:\n");
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));
		
		//get position angle of proper motion
		
		//heliocentric transverse velocity
		vtx = vx-vxsun-vrx;
		vty = vy-vysun-vLSRtemp-vry;
		vtz = vz-vzsun-vrz;
		if (radiococo) printf("\nTransverse velocity:\n");
		if (radiococo) printf("vtx = %f\tvty = %f\tvtz = %f\tvt = %f [km/s] (heliocentric)\n", vtx, vty, vtz, sqrt(vtx*vtx+vty*vty+vtz*vtz));
		
		//get tangential vector pointing to ENP
		FAK = -((xENP-xsun)*(x-xsun)+yENP*y+zENP*z)/(pow(x-xsun,2)+y*y+z*z);
		xdelta = FAK*(x-xsun)+(xENP-xsun);
		ydelta = FAK*y+yENP;
		zdelta = FAK*z+zENP;
		
		//determine distance (pos or neg) of Xobject + Vt from plane connecting ENP, Xobject and observer for position angle
		nx = y*zENP-z*yENP;
		ny = z*(xENP-xsun)-(x-xsun)*zENP;
		nz = (x-xsun)*yENP-y*(xENP-xsun);
		dvt = nx*(x+vtx)+ny*(y+vty)+nz*(z+vtz)-nx*xsun;
		
		//get position angle of proper motion with respect to tangential vector pointing to ENP
		if (dvt <= 0)
			PArad = acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		else
			PArad = 2.0*PI-acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		
		if (radiococo) printf("\nProper motion and position angle:\n");
		if (radiococo) printf("mu = %f\tPA = %f\n", mu, PArad);
		
	} else if (vcoordtype == 3) {
		
		if (radiococo) printf("\nCartesian velocity:\n");
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s]\n",vx,vy,vz, sqrt(vx*vx+vy*vy+vz*vz));
		if (radiococo) printf("vx = %.3f\tvy = %.3f\tvz = %.3f\tv = %.3f [km/s] (heliocentric)\n",vx-vxsun,vy-vysun-vLSRtemp,vz-vzsun, sqrt(pow(vx-vxsun,2)+pow(vy-vysun-vLSRtemp,2)+pow(vz-vzsun,2)));
		
		//heliocentric radial velocity
		vrsun = ((vx-vxsun)*(x-xsun)+(vy-vysun-vLSRtemp)*y+(vz-vzsun)*z)/sqrt(pow(x-xsun,2)+y*y+z*z);
		
		//get radial velocity in 3d coordinates [km/s]
		vrx = (x - xsun)/dsun_coco*vrsun;
		vry = y/dsun_coco*vrsun;
		vrz = z/dsun_coco*vrsun;
		if (radiococo) printf("\nHeliocentric radial velocity in cartesian coordinates:\n");
		if (radiococo) printf("vrx = %.3f\tvry = %.3f\tvrz = %.3f\tvr = %.3f [km/s] (heliocentric)\n",vrx,vry,vrz,sqrt(vrx*vrx+vry*vry+vrz*vrz));
		
		//get position angle of proper motion
		
		//heliocentric transverse velocity
		vtx = vx-vxsun-vrx;
		vty = vy-vysun-vLSRtemp-vry;
		vtz = vz-vzsun-vrz;
		if (radiococo) printf("\nTransverse velocity:\n");
		if (radiococo) printf("vtx = %f\tvty = %f\tvtz = %f\tvt = %f [km/s] (heliocentric)\n", vtx, vty, vtz, sqrt(vtx*vtx+vty*vty+vtz*vtz));
		
		//get tangential vector pointing to ENP
		FAK = -((xENP-xsun)*(x-xsun)+yENP*y+zENP*z)/(pow(x-xsun,2)+y*y+z*z);
		xdelta = FAK*(x-xsun)+(xENP-xsun);
		ydelta = FAK*y+yENP;
		zdelta = FAK*z+zENP;
		
		//determine distance (pos or neg) of Xobject + Vt from plane connecting ENP, Xobject and observer for position angle
		nx = y*zENP-z*yENP;
		ny = z*(xENP-xsun)-(x-xsun)*zENP;
		nz = (x-xsun)*yENP-y*(xENP-xsun);
		dvt = nx*(x+vtx)+ny*(y+vty)+nz*(z+vtz)-nx*xsun;
		
		//get position angle of proper motion with respect to tangential vector pointing to ENP
		if (dvt <= 0)
			PArad = acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		else
			PArad = 2.0*PI-acos((xdelta*vtx+ydelta*vty+zdelta*vtz)/(sqrt(vtx*vtx+vty*vty+vtz*vtz)*sqrt(xdelta*xdelta+ydelta*ydelta+zdelta*zdelta)));
		
		if (radiococo) printf("\nProper motion and position angle:\n");
		if (radiococo) printf("mu = %f\tPA = %f\n", mu, PArad);
		
		mu = sqrt(vtx*vtx+vty*vty+vtz*vtz);
		mu_delta_coco = mu*cos(PArad);
		mu_alphacosdelta_coco = mu*sin(PArad);
		
	}
	
	
	if (radiococo) printf("\nProper motion:\n");
	if (radiococo) printf("mu_alphacosdelta  = %f\tmu_delta = %f\tmu = %f [km/s]\t PA = %f\n", mu_alphacosdelta_coco, mu_delta_coco, mu, PArad);
	
	vr_coco = (vx*(x-xsun)+vy*y+vz*z)/sqrt(pow(x-xsun,2)+y*y+z*z);
	if (radiococo) printf("\nRadial velocity:\n");
	if (radiococo) printf("vr = %.3f\tvr = %.3f (heliocentric) [km/s]\n", vr_coco, vrsun);
	
	//consistency check with formula for GSR radial velocity from script of Steven Majewski
	if (radiococo) vrLSR = vrsun + (vxsun*cos(brad)*cos(lrad)+vysun*cos(brad)*sin(lrad)+vzsun*sin(brad));
	if (radiococo) vrGSR = vrLSR + vLSRtemp*cos(brad)*sin(lrad);
	if (radiococo) printf("\nConsistency check with formula for Galactic standard of rest (GSR) radial velocity (should be equal to vr):\n");
	if (radiococo) printf("vr_LSR = %f\tvr_GSR = %.3f [km/s]\n", vrLSR, vrGSR);
	
	
	
	//convert back to input units and write to output
	*xtemp = 1000.0*x;
	*(xtemp+1) = 1000.0*y;
	*(xtemp+2) = 1000.0*z;
	
	*dsuntemp = 1000.0*dsun_coco;
	
	*vtemp = vx;
	*(vtemp+1) = vy;
	*(vtemp+2) = vz;
	
	*vrsuntemp = vrsun;
	*vrtemp = vr_coco;
	
	*DECtemp = DECrad*180.0/PI;
	*RAtemp = RArad*180.0/PI;
	
	*btemp = brad*180.0/PI;
	*ltemp = lrad*180.0/PI;
	*lcosbtemp = *ltemp*cos(brad);
	
	*mutemp = mu/(dsun_coco*4.74057);
	*PAtemp = PArad*180.0/PI;
	*mu_deltatemp = mu_delta_coco/(dsun_coco*4.74057);
	*mu_alphacosdeltatemp = mu_alphacosdelta_coco/(dsun_coco*4.74057);
	*mu_alphatemp = *mu_alphacosdeltatemp/cos(DECrad);
	
}

