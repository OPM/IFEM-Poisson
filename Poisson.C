// $Id$
//==============================================================================
//!
//! \file Poisson.C
//!
//! \date Apr 16 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Integrand implementations for Poisson problems.
//!
//==============================================================================

#include "Poisson.h"
#include "FiniteElement.h"
#include "Utilities.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "AnaSol.h"
#include "VTF.h"


Poisson::Poisson (unsigned short int n)
{
  npv = 1; // One primary unknown per node (scalar equation)
  nsd = n;

  kappa = 1.0;

  tracFld = 0;
  fluxFld = 0;
  heatSrc = 0;
}


double Poisson::getHeat (const Vec3& X) const
{
  return heatSrc ? (*heatSrc)(X) : 0.0;
}


double Poisson::getFlux (const Vec3& X, const Vec3& n) const
{
  if (fluxFld)
    return (*fluxFld)(X);
  else if (tracFld)
    return (*tracFld)(X)*n;
  else
    return 0.0;
}


void Poisson::setMode (SIM::SolutionMode mode)
{
  m_mode = mode;

  if (mode == SIM::RECOVERY)
    primsol.resize(1);
  else
    primsol.clear();
}


LocalIntegral* Poisson::getLocalIntegral (size_t nen, size_t,
                                          bool neumann) const
{
  ElmMats* result = new ElmMats();
  switch (m_mode)
  {
    case SIM::STATIC:
      result->rhsOnly = neumann;
      result->withLHS = !neumann;
      result->resize(neumann ? 0 : 1, 1+galerkin.size());
      break;

    case SIM::VIBRATION:
      result->resize(1,0);
      break;

    case SIM::RHS_ONLY:
      result->resize(neumann ? 0 : 1, 1);

    case SIM::RECOVERY:
      result->rhsOnly = true;
      result->withLHS = false;
      break;

    default:
      ;
  }

  result->redim(nen);
  return result;
}


bool Poisson::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  if (!elMat.A.empty())
  {
    // Evaluate the constitutive matrix at this point
    Matrix C, CB;
    this->formCmatrix(C,X);

    // Integrate the coefficient matrix
    CB.multiply(C,fe.dNdX,false,true).multiply(fe.detJxW); // = C*dNdX^T*|J|*w
    elMat.A.front().multiply(fe.dNdX,CB,false,false,true); // EK += dNdX * CB
  }

  // Integrate heat source, if defined
  if (heatSrc && !elMat.b.empty())
    WeakOps::Source(elMat.b.front(), fe, (*heatSrc)(X));

  // Galerkin projections a(u^h,v^h) = a(Pu,v^h) = a(w,v^h)
  for (size_t a = 1; a <= galerkin.size(); a++)
  {
    Vec3 Gw = (*galerkin[a-1])(X) * fe.detJxW;
    fe.dNdX.multiply(Gw.vec(nsd),elMat.b[a],false,true); // b += dNdX * Gw
  }

  return true;
}


bool Poisson::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const
{
  if (!tracFld && !fluxFld)
  {
    std::cerr <<" *** Poisson::evalBou: No heat flux."<< std::endl;
    return false;
  }

  ElmMats& elMat = static_cast<ElmMats&>(elmInt);
  if (elMat.b.empty())
  {
    std::cerr <<" *** Poisson::evalBou: No load vector."<< std::endl;
    return false;
  }

  // Evaluate the Neumann value h = -q(X)*n
  double h = -this->getFlux(X,normal);

  // Store flux value for visualization
  if (fe.iGP < fluxVal.size() && fabs(h) > 1.0e-8)
  {
    fluxVal[fe.iGP].first = X;
    fluxVal[fe.iGP].second += h*normal;
  }

  // Integrate the Neumann value
  elMat.b.front().add(fe.N,h*fe.detJxW);

  return true;
}


void Poisson::initIntegration (size_t, size_t nBp)
{
  fluxVal.clear();
  fluxVal.resize(nBp,std::make_pair(Vec3(),Vec3()));
}


bool Poisson::writeGlvT (VTF* vtf, int iStep, int& geoBlk, int& nBlock) const
{
  if (fluxVal.empty())
    return true;
  else if (!vtf)
    return false;

  // Write boundary heat flux as discrete point vectors to the VTF-file
  return vtf->writeVectors(fluxVal,geoBlk,++nBlock,"Heat flux",iStep);
}


bool Poisson::formCmatrix (Matrix& C, const Vec3&, bool invers) const
{
  C.diag(invers && kappa != 0.0 ? 1.0/kappa : kappa, nsd);
  return true;
}


bool Poisson::evalSol (Vector& q, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const
{
  if (primsol.empty() || primsol.front().empty())
  {
    std::cerr <<" *** Poisson::evalSol: No primary solution."<< std::endl;
    return false;
  }

  Vector eV;
  int ierr = utl::gather(MNPC,1,primsol.front(),eV);
  if (ierr > 0)
  {
    std::cerr <<" *** Poisson::evalSol: Detected "<< ierr
	      <<" node numbers out of range."<< std::endl;
    return false;
  }

  return this->evalSol(q,eV,fe.dNdX,X);
}


bool Poisson::evalSol (Vector& q, const Vector& eV,
                       const Matrix& dNdX, const Vec3& X) const
{
  if (eV.size() != dNdX.rows())
  {
    std::cerr <<" *** Poisson::evalSol: Invalid solution vector."
              <<"\n     size(eV) = "<< eV.size() <<"   size(dNdX) = "
              << dNdX.rows() <<","<< dNdX.cols() << std::endl;
    return false;
  }

  // Evaluate the constitutive matrix at this point
  Matrix C, CB;
  this->formCmatrix(C,X);

  // Evaluate the heat flux vector
  CB.multiply(C,dNdX,false,true).multiply(eV,q); // q = C*dNdX^T*eV
  q *= -1.0;

  return true;
}


std::string Poisson::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "u";

  return prefix + std::string(" u");
}


std::string Poisson::getField2Name (size_t i, const char* prefix) const
{
  if (i >= nsd) return 0;

  static const char* s[3] = { "q_x","q_y","q_z" };
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


NormBase* Poisson::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new PoissonNorm(*const_cast<Poisson*>(this),
			   asol->getScalarSecSol());
  else
    return new PoissonNorm(*const_cast<Poisson*>(this));
}


void Poisson::clearGalerkinProjections ()
{
  for (size_t i = 0; i < galerkin.size(); i++)
    delete galerkin[i];

  galerkin.clear();
}


PoissonNorm::PoissonNorm (Poisson& p, VecFunc* a) : NormBase(p), anasol(a)
{
  nrcmp = myProblem.getNoFields(2);
}


bool PoissonNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
			   const Vec3& X) const
{
  Poisson& problem = static_cast<Poisson&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the inverse constitutive matrix at this point
  Matrix Cinv;
  problem.formCmatrix(Cinv,X,true);

  // Evaluate the finite element heat flux field
  Vector sigmah, sigma, error;
  if (!problem.evalSol(sigmah,pnorm.vec.front(),fe.dNdX,X))
    return false;

  size_t ip = 0;
  // Integrate the energy norm a(u^h,u^h)
  pnorm[ip++] += sigmah.dot(Cinv*sigmah)*fe.detJxW;
  // Evaluate the temperature field
  double u = pnorm.vec.front().dot(fe.N);
  // Integrate the external energy (h,u^h)
  pnorm[ip++] += problem.getHeat(X)*u*fe.detJxW;

  if (anasol)
  {
    // Evaluate the analytical heat flux
    sigma.fill((*anasol)(X).ptr(),nrcmp);
    // Integrate the energy norm a(u,u)
    pnorm[ip++] += sigma.dot(Cinv*sigma)*fe.detJxW;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    error = sigma - sigmah;
    pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
  }

  size_t i, j;
  for (i = 0; i < pnorm.psol.size(); i++)
    if (!pnorm.psol[i].empty())
    {
      // Evaluate projected heat flux field
      Vector sigmar(nrcmp);
      for (j = 0; j < nrcmp; j++)
	sigmar[j] = pnorm.psol[i].dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += sigmar.dot(Cinv*sigmar)*fe.detJxW;
      // Integrate the estimated error in energy norm a(u^r-u^h,u^r-u^h)
      error = sigmar - sigmah;
      pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;

      if (anasol)
      {
	// Integrate the error in the projected solution a(u-u^r,u-u^r)
	error = sigma - sigmar;
	pnorm[ip++] += error.dot(Cinv*error)*fe.detJxW;
	ip++; // Make room for the local effectivity index here
      }
    }

  return true;
}


bool PoissonNorm::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
			   const Vec3& X, const Vec3& normal) const
{
  Poisson& problem = static_cast<Poisson&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface heat flux
  double h = problem.getFlux(X,normal);
  // Evaluate the temperature field
  double u = pnorm.vec.front().dot(fe.N);

  // Integrate the external energy (h,u^h)
  pnorm[1] += h*u*fe.detJxW;
  return true;
}


bool PoissonNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as sqrt(a(e^r,e^r)/a(e,e))
  // with e^r = u^r - u^h  and  e = u - u^h
  for (size_t ip = 7; ip < pnorm.size(); ip += 4)
    pnorm[ip] = sqrt(pnorm[ip-2] / pnorm[3]);

  return true;
}


void PoissonNorm::addBoundaryTerms (Vectors& gNorm, double energy) const
{
  gNorm.front()[1] += energy;
}


size_t PoissonNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else
    return anasol ? 4 : 2;
}


std::string PoissonNorm::getName (size_t i, size_t j, const char* prefix) const
{
  if (i == 0 || j == 0 || j > 4)
    return this->NormBase::getName(i,j,prefix);

  static const char* s[8] = {
    "a(u^h,u^h)^0.5",
    "(h,u^h)^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h",
    "a(u^r,u^r)^0.5",
    "a(e,e)^0.5, e=u^r-u^h",
    "a(e,e)^0.5, e=u-u^r",
    "effectivity index"
  };

  size_t k = i > 1 ? j+3 : j-1;

  if (!prefix)
    return s[k];

  return prefix + std::string(" ") + s[k];
}


bool PoissonNorm::hasElementContributions (size_t i, size_t j) const
{
  return i > 1 || j != 2;
}
