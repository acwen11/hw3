#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdlib.h>
#include <math.h>

using namespace std;

double xinit,yinit,zinit,vxinit,vyinit,vzinit;

void rk4(double* xvec, double step, double kappa) ;
void derivs(double* xvec, double* dxvec, double kappa) ;

double potential(double xx, double yy, double zz, double kappa, double* potderivs) {

  //Below, make sure the potential and its gradient are consistent!  The gradient is in "potderivs", and we
  //return the potential itself via the function

  double r = sqrt(xx*xx+yy*yy+zz*zz);
  potderivs[0] = -1*xx/pow(r,3)+3*kappa*xx/pow(r,5);
  potderivs[1] = -1*yy/pow(r,3)+3*kappa*yy/pow(r,5);
  potderivs[2] = -1*zz/pow(r,3)+3*kappa*zz/pow(r,5);
  return -1.0/r+kappa/pow(r,3);
}

// Compile the code and call it with ./orbit ncase Lz L
// where ncase is the potential model used in "potential" above
// Lz is the specific azimuthal angular momentum
// and L is the relative angular momentum


int main(int argc, char** argv) {

 
  double kappa = atof(argv[1]);
  
  //The initial position lies in the x-z plane, as determiined by the angular momentum values
  //xinit = lzinit/linit;
  //zinit = sqrt(1.0-xinit*xinit);

  //change the initial radius
  double rinit = atof(argv[2]);
  //cout<<"rinit?:"<<rinit<<endl;

  xinit=rinit;
  yinit=0.0;
  zinit=0.0;

  double potd[3];
  double potval =potential(xinit,0.0,0.0,kappa,potd);
  //cout<<"potval:"<<potval<<endl;
    

  
  // The inital velocity is in the y-direction, and is set by the angular momentum as well
  //  double eccentricity = sqrt(1-linit*linit);

  vyinit=1.1;
  //cout<<"vyinit:"<<vyinit<<endl;

  vxinit=0.0;
  vzinit=0.0;

  double tf = 100;
  double step0=1.0e-4;
  int np = (int) (tf/step0);

  double step;

  double* t = new double[np];
  double* x = new double[np];
  double* y = new double[np];
  double* z = new double[np];
  double* vx = new double[np];
  double* vy = new double[np];
  double* vz = new double[np];

  double potderivs[3];
  
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
    
    rk4(vec,step,kappa);
    x[i]=vec[0];
    y[i]=vec[1];
    z[i]=vec[2];
    vx[i]=vec[3];
    vy[i]=vec[4];
    vz[i]=vec[5];
    
    t[i]=t[i-1]+step;

    //    cout<<"orbit"<<x[i]<<" "<<y[i]<<endl;
  }


  double theta=0.0;
  int norbitphi=0;
  int norbitr=0;
  for (j=1; j<i-1; j++) {

    double thetanew=atan2(y[j],x[j]);
    if(thetanew<0)thetanew+=2*M_PI;
    if(y[j]>0 && y[j-1]<0) {
      norbitphi+=1;
      double omegap = (2*M_PI*norbitphi+thetanew)/t[j];
      //      cout<<"phi:"<<j<<" "<<t[j]<<" "<<norbitphi<<" "<<thetanew<<" "<<omegap<<" "<<x[j]<<" "<<y[j]<<endl;
      //cout<<"phi:"<<j<<" "<<omegap<<endl;
    }
	
    double rr0=x[j]*x[j]+y[j]*y[j]+z[j]*z[j];
    double rrm1=x[j-1]*x[j-1]+y[j-1]*y[j-1]+z[j-1]*z[j-1];
    double rrp1=x[j+1]*x[j+1]+y[j+1]*y[j+1]+z[j+1]*z[j+1];
    if((rr0<rrm1) && (rr0<rrp1)) {
	norbitr+=1;
	//cout<<"rad:"<<j<<" "<<t[j]<<" "<<2*M_PI*norbitr/t[j]<<endl;
	//cout<<"rad:"<<j<<" "<<2*M_PI*norbitr/t[j]<<endl;
    }
  }
  
  ofstream outfile;
  outfile.open("orbit.dat");
  for(j=0; j<i; j++) {
    if(j%100==0) {
      double E = 0.5*(vx[j]*vx[j]+vy[j]*vy[j]+vz[j]*vz[j])+potential(x[j],y[j],z[j],kappa,potderivs);
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
	" "<<vx[j]<<" "<<vy[j]<<" "<<vz[j]<<" "<<E<<" "<<LX<<" "<<LY<<" "<<LZ<<" "<<r<<" "<<rcyl<<
	" "<<LL <<" "<<omega<<" "<<eccen<<" "<<evecx<<" "<<evecy<<" "<<evecz<<endl;
      //     if(z[j]>0 && z[j-1]<0 && j>=2)outfile2<<rcyl<<" "<<vrcyl<<endl;
    }
  }
  outfile.close();
  return 0;

}

void rk4(double* xvec, double step, double kappa) {

  double h=step/2.0;    
  double dxvec[6],k1[6],k2[6],k3[6],k4[6],t1[6],t2[6],t3[6];
  
  int i;

  derivs(xvec,dxvec,kappa);
  for (i=0; i<6; i++) {
    k1[i]=step*dxvec[i];
    t1[i] = xvec[i]+0.5*k1[i];
  }

  derivs(t1,dxvec,kappa);
  for (i=0; i<6; i++) {
    k2[i]=step*dxvec[i];
    t2[i] = xvec[i]+0.5*k2[i];
  }
  derivs(t2,dxvec,kappa);	
  for (i=0; i<6; i++) {
    k3[i]=step*dxvec[i];
    t3[i] = xvec[i]+k3[i];
  }

  derivs(t3,dxvec,kappa);
  for (i=0; i<6; i++) {
    k4[i] = step*dxvec[i];
  }
  for (i=0; i<6; i++) xvec[i] += (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
}

void derivs(double* xvec, double* dxvec, double kappa) {
  dxvec[0]=xvec[3];
  dxvec[1]=xvec[4];
  dxvec[2]=xvec[5];

  double xx = xvec[0];
  double yy = xvec[1];
  double zz = xvec[2];

  double potderivs[3];
  double potval=potential(xx,yy,zz,kappa,potderivs);

  dxvec[3]=potderivs[0];
  dxvec[4]=potderivs[1];
  dxvec[5]=potderivs[2];
}
