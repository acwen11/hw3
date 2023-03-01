#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <math.h>

using namespace std;

int cases;
double xinit,yinit,zinit,vxinit,vyinit,vzinit;

void rk4(double* xvec, double step) ;
void derivs(double* xvec, double* dxvec) ;

double phisfact=0.01;
double rtoverr0 = 0.03 * 2.0/3.0;

double potential(double xx, double yy, double zz, double* potderivs) {

  //Below, make sure the potential and its gradient are consistent!  The gradient is in "potderivs", and we
  //return the potential itself via the function

  double r = sqrt(xx*xx+yy*yy+zz*zz);
  if(cases==1) {
    //
    potderivs[0] = -1*xx/pow(r,3);
    potderivs[1] = -1*yy/pow(r,3);
    potderivs[2] = -1*zz/pow(r,3);
    return -1.0/r;
  } else if (cases==2) {
    potderivs[0] = -1*xx;
    potderivs[1] = -2*yy;
    potderivs[2] = -4*zz;
    return 0.5*(xx*xx+2*yy*yy+4*zz*zz);
  } else if (cases==3) {
    potderivs[0] = -1*xx/(xx*xx+yy*yy+zz*zz/pow(0.9,2));
    potderivs[1] = -1*yy/(xx*xx+yy*yy+zz*zz/pow(0.9,2));
    potderivs[2] = -1*zz/pow(0.9,2)/(xx*xx+yy*yy+zz*zz/pow(0.9,2));
    return 0.5*log(xx*xx+yy*yy+zz*zz/pow(0.9,2));
  } else if (cases==4) {
    potderivs[0] = -1*xx/(0.14*0.14+xx*xx+yy*yy/0.9/0.9);
    potderivs[1] = -1*yy/0.9/0.9/(0.14*0.14+xx*xx+yy*yy/0.9/0.9);
    potderivs[2] = 0.0;
    return 0.5*log(0.14*0.14+xx*xx+yy*yy/0.9/0.9);
  } else if (cases==5) {
    potderivs[0] = -1*(xx-0.2*pow(xx,3));
    potderivs[1] = -0.8*(yy-0.2*pow(yy,3));
    potderivs[2] = -0.6*(zz-0.2*pow(zz,3));
    return 0.5*(xx*xx*(1.0-0.1*xx*xx)+0.8*yy*yy*(1.0-0.1*yy*yy)+0.6*zz*zz*(1.0-0.1*zz*zz));
  } else if (cases==6) {
    potderivs[0] = -1*xx/(xx*xx+yy*yy/0.9+zz*zz/pow(0.9,2));
    potderivs[1] = -1*yy/0.9/(xx*xx+0.9*yy*yy/0.9+zz*zz/pow(0.9,2));
    potderivs[2] = -1*zz/pow(0.9,2)/(xx*xx+yy*yy/0.9+zz*zz/pow(0.9,2));
    return 0.5*log(xx*xx+yy*yy/0.9+zz*zz/pow(0.9,2));
  } else if (cases==7) {
    potderivs[0] = -1*xx/pow(r,3)-phisfact*(xx/r-rtoverr0*(xx*(xx*xx+yy*yy+4*zz*zz)/pow(r,3)));
    potderivs[1] = -1*yy/pow(r,3)-phisfact*(yy/r-rtoverr0*(yy*(xx*xx+yy*yy+4*zz*zz)/pow(r,3)));
    potderivs[2] = -1*zz/pow(r,3)-phisfact*(zz/r+rtoverr0*(zz*(5*xx*xx+5*yy*yy+2*zz*zz)/pow(r,3)));
    return -1.0/r+phisfact*(r+rtoverr0/r*(2*zz*zz-xx*xx-yy*yy));
  }
  return 0;
}

// Compile the code and call it with ./orbit ncase Lz L
// where ncase is the potential model used in "potential" above
// Lz is the specific azimuthal angular momentum
// and L is the relative angular momentum


int main(int argc, char** argv) {

  cases = atoi(argv[1]);
  //double lzinit = atof(argv[2]);
  //double linit = atof(argv[3]);


  ////The initial position lies in the x-z plane, as determiined by the angular momentum values
  //xinit = lzinit/linit;
  //zinit = sqrt(1.0-xinit*xinit);

  //// The inital velocity is in the y-direction, and is set by the angular momentum as well
  //double eccentricity = sqrt(1-linit*linit);
  //vyinit = sqrt(1-eccentricity);

  /* Allen's Low IQ Version */
  /* ./orbit ncase R z E Lz */


  double potderivs[3];

  xinit = atof(argv[2]);
  zinit = atof(argv[3]);
  double E = atof(argv[4]);
  double lzinit = atof(argv[5]);

  vzinit = 0.0;
  yinit = 0.0;
  vyinit = lzinit / xinit;
  vxinit = sqrt(2 * (E - potential(xinit, yinit, zinit, potderivs)) - vyinit * vyinit);

  cout << "Initial data" << " " << xinit << " " << yinit << " " << zinit << "\n" << endl;
  cout << vxinit << " " << vyinit << " " << vzinit << " " << potential(xinit, yinit, zinit, potderivs) <<  "\n" << endl;
  double tf = 700; //5000; 
  double step0=1.0e-3;
  int np = (int) (tf/step0);

  double step;

  double* t = new double[np];
  double* x = new double[np];
  double* y = new double[np];
  double* z = new double[np];
  double* vx = new double[np];
  double* vy = new double[np];
  double* vz = new double[np];

  
  int i,j;

  //initial data
  t[0]=0.0;
  x[0]=xinit;
  y[0]=yinit;
  z[0]=zinit;
  vx[0]=vxinit;
  vy[0]=vyinit;
  vz[0]=vzinit;

  double vec[6];
  vec[0]=xinit;
  vec[1]=yinit;
  vec[2]=zinit;
  vec[3]=vxinit;
  vec[4]=vyinit;
  vec[5]=vzinit;
  
  for (i=1; i<np; i++) {

    //The routine uses a variable step size along with an RK4 timestepper
    double rr2=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
    step = step0*sqrt(rr2);
    
    rk4(vec,step);
    x[i]=vec[0];
    y[i]=vec[1];
    z[i]=vec[2];
    vx[i]=vec[3];
    vy[i]=vec[4];
    vz[i]=vec[5];
    
    t[i]=t[i-1]+step;
  }

  
  ofstream outfile,outfile2;
  outfile.open("orbit.dat");
  //outfile2.open("orbit2.dat");
  for(j=0; j<i; j++) {
    if(j%100==0) {
      double E = 0.5*(vx[j]*vx[j]+vy[j]*vy[j]+vz[j]*vz[j])+potential(x[j],y[j],z[j],potderivs);
      double LX = y[j]*vz[j]-z[j]*vy[j];
      double LY = z[j]*vx[j]-x[j]*vz[j];
      double LZ = x[j]*vy[j]-y[j]*vx[j];
      double LL = sqrt(LX*LX+LY*LY+LZ*LZ);
      double r = sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]);
      double rcyl=sqrt(x[j]*x[j]+y[j]*y[j]);
      double vr = (x[j]*vx[j]+y[j]*vy[j]+z[j]*vz[j])/r;
      double vrcyl = (x[j]*vx[j]+y[j]*vy[j])/rcyl;

      //This is the "eccentricity vector" for the orbit
      double evecx = (vy[j]*LZ-vz[j]*LY)-x[j]/r;
      double evecy = (vz[j]*LX-vx[j]*LZ)-y[j]/r;
      double evecz = (vx[j]*LY-vy[j]*LX)-z[j]/r;
      double eccen = sqrt(evecx*evecx+evecy*evecy+evecz*evecz);

	//The longitude of the ascending node has coords (-LY, LX,0)
      double omega = acos((-1*LY*evecx+LX*evecy)/eccen/sqrt(LX*LX+LY*LY));
	
	
      outfile<<t[j]<<" "<<x[j]<<" "<<y[j]<<" "<<z[j]<<
	" "<<vx[j]<<" "<<vy[j]<<" "<<vz[j]<<" "<<E<<" "<<LX<<" "<<LY<<" "<<LZ<<" "<<r<<
	" "<<LL <<" "<<omega<<" "<<eccen<<" "<<evecx<<" "<<evecy<<" "<<evecz<<endl;
      //if(z[j]>0 && z[j-1]<0 && j>=2) outfile2<<rcyl<<" "<<vrcyl<<endl;
    }
  }
  outfile.close();
  //outfile2.close();
  return 0;

}

void rk4(double* xvec, double step) {

  double h=step/2.0;    
  double dxvec[6],k1[6],k2[6],k3[6],k4[6],t1[6],t2[6],t3[6];
  
  int i;

  derivs(xvec,dxvec);
  for (i=0; i<6; i++) {
    k1[i]=step*dxvec[i];
    t1[i] = xvec[i]+0.5*k1[i];
  }

  derivs(t1,dxvec);
  for (i=0; i<6; i++) {
    k2[i]=step*dxvec[i];
    t2[i] = xvec[i]+0.5*k2[i];
  }
  derivs(t2,dxvec);	
  for (i=0; i<6; i++) {
    k3[i]=step*dxvec[i];
    t3[i] = xvec[i]+k3[i];
  }

  derivs(t3,dxvec);
  for (i=0; i<6; i++) {
    k4[i] = step*dxvec[i];
  }
  for (i=0; i<6; i++) xvec[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
}

void derivs(double* xvec, double* dxvec) {
  dxvec[0]=xvec[3];
  dxvec[1]=xvec[4];
  dxvec[2]=xvec[5];

  double xx = xvec[0];
  double yy = xvec[1];
  double zz = xvec[2];

  double potderivs[3];
  double potval=potential(xx,yy,zz,potderivs);

  dxvec[3]=potderivs[0];
  dxvec[4]=potderivs[1];
  dxvec[5]=potderivs[2];
}
