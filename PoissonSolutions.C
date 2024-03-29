// $Id$
//==============================================================================
//!
//! \file PoissonSolutions.C
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Analytic solutions for Poisson problems.
//!
//==============================================================================

#include "PoissonSolutions.h"
#include "Vec3.h"

#include <cmath>


Vec3 Square2D::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  Vec3 temp;
  temp.x = -pi*sin(pi*x)*(2.0-y);
  temp.y = -cos(pi*x);

  return temp;
}

double Square2DHeat::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  return -pi*pi*cos(pi*x)*(2.0-y);
}


Vec3 LshapePoisson::evaluate (const Vec3& X) const
{
/*
  double x = X.x;
  double y = X.y;
  double R = hypot(x,y);
  double pi = M_PI;
  double Rp = pow(R,-1.0/3.0);
  double th = x > 0.0 ? asin(y/R) : pi-asin(y/R);
  double frac = 2.0/3.0;

  Vec3 temp;
  temp.x =  frac*Rp*sin(th/3.0);
  temp.y = -frac*Rp*cos(th/3.0);

  return temp;
*/
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  double r2 = x*x+y*y;
  double theta = atan2(y,x);
  if(theta<=0) theta += 2*pi;
  if (r2 < 1e-16) {  // truncate the singularity to avoid NaN values
    r2 = 1e-16; 
    theta = 2*pi;
  }

  Vec3 temp;
  temp.x = (2.0/3) * (cos(2.0/3*theta + pi/6)*x + sin(2.0/3*theta + pi/6)*y) / pow(r2, 2.0/3);
  temp.y = (2.0/3) * (cos(2.0/3*theta + pi/6)*y - sin(2.0/3*theta + pi/6)*x) / pow(r2, 2.0/3);

  return temp;

}


Vec3 SquareSinus::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  Vec3 temp;
  temp.x = 2*pi*cos(2*pi*x)*sin(2*pi*y);
  temp.y = 2*pi*sin(2*pi*x)*cos(2*pi*y);

  return temp;
}

double SquareSinusSource::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;

  return -2*2*2*pi*pi*sin(2*pi*x)*sin(2*pi*y);
}


Vec3 PoissonInteriorLayer::evaluate (const Vec3& X) const
{
  double x = X.x - 1.25;
  double y = X.y + 0.25;
  double pi = M_PI;
  Vec3 result;

  double root = hypot(x,y);
  double u =  SLOPE*(root-pi/3.0);
  result.x = -SLOPE*x / root / (1.0+u*u);
  result.y = -SLOPE*y / root / (1.0+u*u);
  return result;
}


/*!
  \class PoissonInteriorLayerSol

  Picked up from http://hpfem.org/hermes2d/doc/src/benchmarks.html
*/
double PoissonInteriorLayerSol::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double pi = M_PI;
  return atan(SLOPE*(hypot(x-1.25,y+0.25)-pi/3.0));
}

double PoissonInteriorLayerSource::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;

  /******       maple did this for me (nabla of the u above)       ************/
  return SLOPE * pow(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1), -0.3e1 / 0.2e1) * pow(0.2e1 * x - 0.250e1, 0.2e1) / (0.1e1 + SLOPE * SLOPE * pow(sqrt(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) - 0.3141592654e1 / 0.3e1, 0.2e1)) / 0.4e1 - 0.2e1 * SLOPE * pow(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1), -0.1e1 / 0.2e1) / (0.1e1 + SLOPE * SLOPE * pow(sqrt(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) - 0.3141592654e1 / 0.3e1, 0.2e1)) + pow(SLOPE, 0.3e1) / (pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) * pow(0.2e1 * x - 0.250e1, 0.2e1) * pow(0.1e1 + SLOPE * SLOPE * pow(sqrt(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) - 0.3141592654e1 / 0.3e1, 0.2e1), -0.2e1) * (sqrt(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) - 0.3141592654e1 / 0.3e1) / 0.2e1 + SLOPE * pow(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1), -0.3e1 / 0.2e1) * pow(0.2e1 * y + 0.50e0, 0.2e1) / (0.1e1 + SLOPE * SLOPE * pow(sqrt(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) - 0.3141592654e1 / 0.3e1, 0.2e1)) / 0.4e1 + pow(SLOPE, 0.3e1) / (pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) * pow(0.2e1 * y + 0.50e0, 0.2e1) * pow(0.1e1 + SLOPE * SLOPE * pow(sqrt(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) - 0.3141592654e1 / 0.3e1, 0.2e1), -0.2e1) * (sqrt(pow(x - 0.125e1, 0.2e1) + pow(y + 0.25e0, 0.2e1)) - 0.3141592654e1 / 0.3e1) / 0.2e1;
  /****************************************************************************/
}

Vec3 PoissonWaterfall::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double z = X.z;
  double r = sqrt(x*x + y*y + z*z);
  double coef = 1.0 / r / epsilon / (1+pow((r-0.5)/epsilon,2));
  Vec3 result;

  result.x = -coef*x;
  result.y = -coef*y;
  result.z = -coef*z;
  return result;
}


/*!
  \class PoissonWaterfallSol

  3D waterfall. Should be a sharp (as defined by epsilon) edge at r=0.5
*/
double PoissonWaterfallSol::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double z = X.z;
  double r = sqrt(x*x + y*y + z*z);
  return atan((r-0.5)/epsilon);
}

double PoissonWaterfallSource::evaluate (const Vec3& X) const
{
  double x  = X.x;
  double y  = X.y;
  double z  = X.z;
  double r2 = x*x + y*y + z*z;
  double r  = sqrt(r2);
  double e  = epsilon;

  return 8*e * (-4*e*e+2*r-1) / r / pow(-4*e*e -4*r2 + 4*r - 1, 2);

}


Vec3 PoissonCube::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double z = X.z;
  double pi = M_PI;

  Vec3 temp;
  temp.x = -pi*cos(x*pi)*sin(y*pi)*sin(z*pi);
  temp.y = -pi*sin(x*pi)*cos(y*pi)*sin(z*pi);
  temp.z = -pi*sin(x*pi)*sin(y*pi)*cos(z*pi);

  return temp;
}

double PoissonCubeSource::evaluate (const Vec3& X) const
{
  double x = X.x;
  double y = X.y;
  double z = X.z;
  double pi = M_PI;

  return 3.0*pi*pi*sin(x*pi)*sin(y*pi)*sin(z*pi);
}


Vec3 PoissonLine::evaluate (const Vec3& X) const
{
  double x = X.x;
  double pi = M_PI;

  Vec3 temp;
  temp.x = -(pi/L)*cos(x*pi/L);

  return temp;
}

double PoissonLineSource::evaluate (const Vec3& X) const
{
  double x = X.x;
  double pi = M_PI;

  return (pi*pi)/(L*L)*sin(pi*x/L);
}
