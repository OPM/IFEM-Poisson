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
#include "ElmMats.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "ExprFunctions.h"
#include "Utilities.h"
#include "VTF.h"
#include "tinyxml.h"


Poisson::Poisson (unsigned short int n) : IntegrandBase(n)
{
  kappa = 1.0;

  tracFld = nullptr;
  fluxFld = nullptr;
  heatSrc = nullptr;
  reacInt = nullptr;
  extEner = false;

  normIntegrandType = ELEMENT_CORNERS;
}


bool Poisson::parse (const TiXmlElement* elem)
{
  if (!elem) return false;

  if (elem->FirstChildElement("residual"))
    normIntegrandType |= SECOND_DERIVATIVES;
  else if (strcasecmp(elem->Value(),"galerkin"))
    return false;
  else if (elem->FirstChild() && elem->FirstChild()->Value())
    galerkin.push_back(new VecFuncExpr(elem->FirstChild()->Value()));

  return true;
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

  primsol.resize(mode >= SIM::RHS_ONLY ? 1 : 0);
}


GlobalIntegral& Poisson::getGlobalInt (GlobalIntegral* gq) const
{
  if (m_mode == SIM::RHS_ONLY && reacInt)
    return *reacInt;

  return this->IntegrandBase::getGlobalInt(gq);
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
    // Conductivity scaled by integration point weight at this point
    double cw = kappa*fe.detJxW;
    // Integrate the coefficient matrix // EK += kappa * dNdX * dNdX^T * |J|*w
    elMat.A.front().multiply(fe.dNdX,fe.dNdX,false,true,true,cw);
  }

  if (!elMat.b.empty())
  {
    // Integrate heat source, if defined
    if (heatSrc)
      elMat.b.front().add(fe.N,(*heatSrc)(X)*fe.detJxW); // EV += N*h(x)*|J|*w

    if (m_mode == SIM::RHS_ONLY && !elmInt.vec.empty())
    {
      // Integrate the internal forces based on current solution
      Vector q;
      if (!this->evalSol(q,elmInt.vec.front(),fe.dNdX,X))
        return false;
      if (!fe.dNdX.multiply(q,elMat.b.front(),fe.detJxW,1.0)) // b += dNdX * q
        return false;
    }
  }

  // Galerkin projections a(u^h,v^h) = a(Pu,v^h) = a(w,v^h)
  bool ok = true;
  for (size_t a = 1; a <= galerkin.size() && a < elMat.b.size() && ok; a++)
  {
    Vec3 Gw = (*galerkin[a-1])(X) * fe.detJxW;
    ok = fe.dNdX.multiply(Gw.vec(nsd),elMat.b[a],false,true); // b += dNdX * Gw
  }

  return ok;
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
  if (h == 0.0) return true; // Skip for homogeneous Neumann

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
  static int partID = 1+geoBlk;
  return vtf->writeVectors(fluxVal,partID,geoBlk,++nBlock,"Heat flux",iStep);
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
  // Evaluate the heat flux vector, q = -kappa*du/dX = -kappa*dNdX^T*eV
  if (!dNdX.multiply(eV,q,true))
  {
    std::cerr <<" *** Poisson::evalSol: Invalid solution vector."<< std::endl;
    return false;
  }

  q *= -kappa;
  return true;
}


std::string Poisson::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "u";

  return prefix + std::string(" u");
}


std::string Poisson::getField2Name (size_t i, const char* prefix) const
{
  if (i >= nsd) return "";

  static const char* s[3] = { "q_x","q_y","q_z" };
  if (!prefix) return s[i];

  return prefix + std::string(" ") + s[i];
}


/*!
  \note The Integrand object is allocated dynamically and has to be deleted
  manually when leaving the scope of the pointer variable receiving the
  returned pointer value.
*/

NormBase* Poisson::getNormIntegrand (AnaSol* asol) const
{
  if (asol)
    return new PoissonNorm(*const_cast<Poisson*>(this), normIntegrandType,
                           asol->getScalarSecSol());
  else
    return new PoissonNorm(*const_cast<Poisson*>(this), normIntegrandType);
}


void Poisson::clearGalerkinProjections ()
{
  for (VecFunc* f : galerkin) delete f;
  galerkin.clear();
}


PoissonNorm::PoissonNorm (Poisson& p, int itype, VecFunc* a) :
  NormBase(p), anasol(a), integrandType(itype)
{
  nrcmp = myProblem.getNoFields(2);
  projBou = true;
}


bool PoissonNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
			   const Vec3& X) const
{
  Poisson& problem = static_cast<Poisson&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the finite element heat flux field
  Vector sigmah, sigma, error;
  if (!problem.evalSol(sigmah,pnorm.vec.front(),fe.dNdX,X))
    return false;

  // Evaluate the temperature field
  double u = pnorm.vec.front().dot(fe.N);
  // Evaluate the heat source field
  double h = problem.getHeat(X);

  // Evaluate the inverse conductivity scaled by the integration point weight
  double cwInv = fe.detJxW / problem.getMaterial();

  // Integrate the energy norm a(u^h,u^h)
  pnorm[0] += sigmah.dot(sigmah)*cwInv;
  // Integrate the external energy (h,u^h)
  if (problem.extEner)
    pnorm[1] += h*u*fe.detJxW;

  size_t ip = 2;
  if (anasol)
  {
    // Evaluate the analytical heat flux
    sigma.fill((*anasol)(X).ptr(),nrcmp);
    // Integrate the energy norm a(u,u)
    pnorm[ip++] += sigma.dot(sigma)*cwInv;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    error = sigma - sigmah;
    pnorm[ip++] += error.dot(error)*cwInv;
  }

  for (const Vector& psol : pnorm.psol)
    if (!psol.empty())
    {
      // Evaluate projected heat flux field
      Vector sigmar(nrcmp);
      for (size_t j = 0; j < nrcmp; j++)
        sigmar[j] = psol.dot(fe.N,j,nrcmp);

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += sigmar.dot(sigmar)*cwInv;
      // Integrate the estimated error in energy norm a(u^r-u^h,u^r-u^h)
      error = sigmar - sigmah;
      pnorm[ip++] += error.dot(error)*cwInv;

      // Evaluate the projected heat flux gradient.
      // Notice that the matrix multiplication method used here treats
      // the element vector (psol) as a matrix whose number of columns
      // equals the number of rows in the matrix fe.dNdX.
      Matrix dSigmadX;
      if (!dSigmadX.multiplyMat(psol,fe.dNdX)) // dSigmadX = psol*dNdX
        return false;

      // Evaluate the interior residual of the projected solution
      double Res = h - dSigmadX.trace();
      // Integrate the residual error in the projected solution
      pnorm[ip++] += fe.h*fe.h*Res*Res*fe.detJxW;

      if (anasol)
      {
        // Integrate the error in the projected solution a(u-u^r,u-u^r)
        error = sigma - sigmar;
        pnorm[ip++] += error.dot(error)*cwInv;
        ip += 2; // Make room for the local effectivity indices here
      }
    }
    else if (integrandType & SECOND_DERIVATIVES)
    {
      // Integrate the residual error in the FE solution
      double Res = h;
      for (size_t j = 1; j <= fe.N.size(); j++)
        Res += fe.d2NdX2.trace(j)*pnorm.vec.front()(j);

      pnorm[ip+1] += fe.h*fe.h*Res*Res*fe.detJxW;
      ip += anasol ? 6 : 3; // Dummy entries in order to get norm in right place
    }

  return true;
}


bool PoissonNorm::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
			   const Vec3& X, const Vec3& normal) const
{
  Poisson& problem = static_cast<Poisson&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface heat flux
  double h = -problem.getFlux(X,normal);
  // Evaluate the temperature field
  double u = pnorm.vec.front().dot(fe.N);

  // Integrate the external energy (h,u^h)
  if (problem.extEner)
    pnorm[1] += h*u*fe.detJxW;

  size_t ip = anasol ? 6 : 4;
  for (const Vector& psol : pnorm.psol)
    if (!psol.empty())
    {
      // Evaluate projected heat flux field
      Vec3 sigmar;
      for (size_t j = 0; j < nrcmp; j++)
        sigmar[j] = psol.dot(fe.N,j,nrcmp);

      // Evaluate the boundary jump term
      double Jump = h + sigmar*normal;
      // Integrate the residual error in the projected solution
      pnorm[ip] += 0.5*fe.h*Jump*Jump*fe.detJxW;
      ip += anasol ? 6 : 3;
    }
    else if (integrandType & SECOND_DERIVATIVES)
      ip += anasol ? 6 : 3; // TODO: Add residual jump terms?

  return true;
}


bool PoissonNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as a(e^r,e^r)/a(e,e)
  // with e^r = u^r - u^h  and  e = u - u^h
  for (size_t ip = 9; ip < pnorm.size(); ip += 6)
  {
    pnorm[ip-1] = pnorm[ip-4] / pnorm[3];
    pnorm[ip] = (pnorm[ip-4]+pnorm[ip-3]) / pnorm[3];
  }

  return true;
}


size_t PoissonNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else if (group > 1)
    return anasol ? 6 : 3;
  else
    return anasol ? 4 : 2;
}


std::string PoissonNorm::getName (size_t i, size_t j, const char* prefix) const
{
  if (i == 0 || j == 0 || j > 6)
    return this->NormBase::getName(i,j,prefix);

  static const char* s[10] = {
    "a(u^h,u^h)^0.5",
    "(h,u^h)^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h",
    "a(u^r,u^r)^0.5",
    "a(e,e)^0.5, e=u^r-u^h",
    "res(u^r)^0.5",
    "a(e,e)^0.5, e=u-u^r",
    "effectivity index^*",
    "effectivity index^RES"
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
