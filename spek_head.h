#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <complex.h>							
#include <fftw3.h>  

#include "eispack.h"     // Lineare Algebra, EW u EV berechnen 
#include "normal.h"      // Normalverteilte Zufallszahlen

#include "spek_funk.h"   // Funktions-Deklarationen


// Array
#define NTAU   101   //  siehe "lineshapeArray( )" !!
#define DTAU   50.0  // [fs] 
#define DTAU_0 0.001  
#define NGAM   21       
#define DGAM   0.05  
#define DGAM_0 0.10 // gamma werte [0.5; 1.0] 

//#define NTAU   276 //101//1001//301//501//  siehe "lineshapeArray( )" !!
//#define DTAU   10.//6.// [fs] 
//#define DTAU_0 0.001  //130.// 
//#define NGAM   81     //20   
//#define DGAM   0.0125 //0.025 
//#define DGAM_0 0.0001  //0.1    //0.50 // gamma werte [0.5; 1.0] 

// Groebere Rasterung verfaelscht die Temperaturabhaengigkeit ein bisschen!!


// FFT-Parameter
#define NT 32768//65536// 2^14=16384// anzahl element FFT (musst be 2^n)
                // 16384 => 177, 32768 => 354, 65536 => 708  intervals
  
#define DYNTHEO 1   // Dynamische Theorie: 0 = Markov, 1 = nonMarkov  
#define SPECDENS 1  // 1 = B777 spectral density
#define SCALE_COUP 1.75

#define COUPLINGS 1  // 1 = point dipole approx, 2 = read in couplings

#define SPECIES 2     // 1 = Aestuarii, 2 = Tepidum
#define SITE_3 12200 //200//180  // 12200 if SPECIES == 1, 12180 if SPEC == 2 

#define BCL_8 0.35  // Anteil BCL 8, BCL7+BCL8 = 1
 
#define N 2000  // N mal wuerfeln bei MonteCarlo

#define NTT 10//0000 //262144//65536 // Number of Train/Test-Sets  

#define TAU_0 1000.//2750.     // [fs]   Pure Dephasing, ca. 5.3 cm-1 = 1000 fs

#define LIFETIMEBROADENING 1    //1 Ein, 0 Aus

#define NNATM 4 // numb of N-Atom-Coord per pigment

//Grenzen f�r Spektrumsberechnung
//Nicht "zu knapp" waehlen, wenn ib in dynArray positiv wird, gibts falsche speicherzugriffe...

#define SPEKMIN 11700 //13000  //14600//[cm-1] 
#define SPEKMAX 13500 //17800  //15800// 
#define SPEKDIFF (SPEKMAX-SPEKMIN)

#define NEXC 200           // pigmExcDistrib()


// Float-Makros -----------------------------------------------------------------------------------------------

#define TEMP 4.//298.//77.      // [K]   Achtung: als float x.y schreiben!!  Tep 6K, Aest 4K

#define HRFAK 0.615 //0.4 // Huang-Rhys-Skalierungsfaktor: so w�hlen, dass S1+S2 = Huang Rhys Faktor = 0.8
                      
#define BETA 0. //7. //0.    //Winkel zwischen qy und NB-ND (N-Atome im Pigment), nur PDIP(Punktdipol) !!!

// beta > 0 Drehung in mathem. Drehsinn, von 13-1-keto-Gruppe weg
// beat < 0 Drehung in Richtung 13-1-keto-Gruppe

//Vakuum Dipolstaerke fuer QY Uebergang [Einheit Debye] von 1 Pigment (�m aus denen dann �M berechnet wird)

#define DIPST sqrt(37.1) //3.6//4.0 // [DEBYE]   // Chl_a 4.0, Chl_b 3.6, BChl_a sqrt(37.1)
 
//-------------Site Energies----------------------------

// Werte aus Adolphs & Renger 2006 fuer Tepidum
// 12445. 12520. 12205. 12335. 12490. 12640. 12450.

#define SITE_MIN 12200
#define SITE_MAX 12800

#define SITE_8 12700

// Breite der Gauss-Vert der site-Energies   // [1/cm]
#define FWHM_0 100. //111.// 90.0 //100. //0.01  // [1/cm]
#define FWHM_1 100. //157.//167.2 //100. //0.01
#define FWHM_2 100. //195.//182.0 //100. //0.01
#define FWHM_3 100. //126.//122.2 //100. //0.01
#define FWHM_4 100. //111.// 89.7 //100. //0.01
#define FWHM_5 100. //172.//142.2 //100. //0.01
#define FWHM_6 100. //147.//118.2 //100. //0.01 
#define FWHM_7 100. //281.//227.2 //100. //0.01


// ------------ Burning ----------------------------------

#define E_BURN 12150.//15047. // [1/cm], Chl_b: 15047=664.6, 15239=656.2, 15293=653.9 // BurnFreq 
                      //         Chl_a: 14663=682.0, 14925=670.0, 15035=665.1   

#define WBZL 7.5// 2.5// 7.5 //10.//5.     // [1/cm] BurnWidth in ZL 

 
#define DIELEK 0.8  //Dielktr. Abschirm. Faktor, fuer PointDip-Calc und TransCharge (empty cavity faktor f)
                    //0.615 empty cav faktor fuer eps=2.56, 0.796 aus mflx mit eps=2 

#define NINDEX 1.414213562 // Brechungsindex, n = sqrt(epsilon), und epsilon=2

#define DIPFAK 1.0  //Skalieren der Dipolst�rke unabh vom Dielektrikum

#define E_LAMB HRFAK*102 //[1/cm] Reorganisationsenergie E_lambda (4.17)



//--------------Fuer SpecDens J(w) aus Renger Marcus 2002 (2.16)

#define S1 HRFAK*0.8  // SpecDens Renger Marcus (2.16)
#define S2 HRFAK*0.5

#define W1 1.048294157e-4 // [fs-1]   // = 0.5565 1/cm * WZ2WFS => 1/fs Renger Marcus (2.16) 
#define W2 3.646240546e-4             // = 1.9357 1/cm * WZ2WFS => 1/fs Renger Marcus (2.16) 



#define S31 HRFAK*0.8 //HRFAK*0.9   //SpecDens Renger Marcus modified to 3 Terms
#define S32 HRFAK*0.3 //HRFAK*0.4
#define S33 HRFAK*0.2 //0.0

#define W31 1.04e-4 //1.04e-4 //1.048294157e-4 // [fs-1]   // = 0.5565 1/cm * WZ2WFS => 1/fs Renger Marcus (2.16) 
#define W32 3.64e-4 //3.04e-4 //3.646240546e-4             // = 1.9357 1/cm * WZ2WFS => 1/fs Renger Marcus (2.16) 
#define W33 7.64e-4 //0.15e-4


#define W1CM 0.5565         // [cm-1] dasselbe wie oben
#define W2CM 1.9357         // [cm-1]

#define PI    3.141592654
#define TWOPI 6.283185307
#define PISQ  9.869604401  // PI^2

#define EXCLENG (1.0*(SPEKMAX-SPEKMIN))/NEXC


// FFT-Parameter
#define DELTA 0.2 //1.0 //5.0 // t-steps FFT [fs]  // Freq-Aufloesung reziprok!  // Max. 5.0 !! 
#define WSTEP (TWOPI/(DELTA*NT))  // [fs-1] Schrittweite in Frequenz-Domaene, egal ob NT oder (NT-1)



// Makros:

#define SQ(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define ABS(x) sqrt( ((x) * (x)) )



// Conversion Factors -----------------------------------------------------------------------------

#define EVWZ   8.06554                // (meV => 1/cm)
#define WZ2WFS 1.88369796e-4          // E(1/cm)*WZ2WFS = w(1/fs), WZ2WFS=2*Pi*c*100/e+15
#define ALPHFAC 4.340277778e-12  // Conversion [Debye^2*s/m]*4.34e-12 = e^2*Ang*s, 4.34e-12=10^-10/4.8^2



// Basic Constants

#define HBAR    0.6582119282           // eV*fs   (hbar/e)
#define HBDKB   7638.23301             // hquer/kB in fs*K
#define HBARC   1.973270524e-5         // hquer*c in meV*m 
#define CVAK    2.99792458e+8          // vacuum speed of light m/s

#define ECHARGE 1.60217733e-19
#define EPS_0   8.854187818e-12
#define ANG     1.0e-10
#define MILLI   0.001
#define DEBFAC  4.8         // Dipolst�rke: e*Ang = 4.8 Debye

#define FWHMSIG 2.354820045 // 2*sqrt(2*ln(2)), fwhm=2*sqrt(2*ln(2))*sig (Gauss)  




