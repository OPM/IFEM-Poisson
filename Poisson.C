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

#include "Fields.h"
#include "Poisson.h"

#include "AnaSol.h"
#include "ElmMats.h"
#include "ElmNorm.h"
#include "ExprFunctions.h"
#include "Function.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "Vec3Oper.h"
#include "VTF.h"

// TODO: Are all these really necessary??
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ext/alloc_traits.h>
#include <iostream>
#include <strings.h>
#include "tinyxml2.h"
#include <utility>


Poisson::Poisson (unsigned short int n) : IntegrandBase(n)
{
  kappaC  = 1.0;
  fluxFld = heatSrc = kappaF = nullptr;
  tracFld = nullptr;
  reacInt = nullptr;
  dualRHS = nullptr;
  aSol    = nullptr;
  extEner = false;

  normIntegrandType = ELEMENT_CORNERS;
}


bool Poisson::parse (const tinyxml2::XMLElement* elem)
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

  if (mode == SIM::STATIC || mode == SIM::NORMS)
    primsol.resize(1+dualFld.size());
  else if (mode >= SIM::RHS_ONLY)
    primsol.resize(1);
  else
    primsol.clear();

#ifdef INT_DEBUG
  std::cout <<"Poisson::setMode: "<< m_mode
            <<"  size(primsol) = "<< primsol.size() << std::endl;
#endif
}


void Poisson::addExtrFunction (FunctionBase* extr)
{
  if (dynamic_cast<RealFunc*>(extr))
  {
    dualFld.push_back(extr);
    if (!dualRHS)
      dualRHS = extr;
  }
  else
  {
    dualFld.clear();
    dualRHS = nullptr;
  }
}


bool Poisson::initElement (const std::vector<int>& MNPC,
                           const FiniteElement&, const Vec3& XC,
                           size_t, LocalIntegral& elmInt)
{
  // In case of global Lagrange multiplier, the last element node in MNPC
  // needs to be omitted when extracting the element solution vectors
  if (!this->initElement2(MNPC,elmInt.vec,setIntegratedSol))
    return false;

  for (size_t i = 1; i < elmInt.vec.size(); i++)
  {
    FunctionBase* extrFunc = nullptr;
    if (i == 1 && m_mode == SIM::STATIC)
      extrFunc = dualRHS;
    else if (i <= dualFld.size() && m_mode >= SIM::RECOVERY)
      extrFunc = dualFld[i-1];

    if (extrFunc && !extrFunc->inDomain(XC))
      elmInt.vec[i].clear(); // Erase extraction field values if outside domain
#ifdef INT_DEBUG
    else
      std::cout <<"Poisson::initElement: Point "<< XC
                <<" is inside domain "<< i << std::endl;
#endif
  }

  return true;
}


void Poisson::setSecondaryInt (GlobalIntegral* gq)
{
  delete reacInt;
  reacInt = gq;
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
      result->resize(neumann ? 0 : 1,
                     neumann ? 1 : (dualRHS ? 2 : 1+galerkin.size()));
      break;

    case SIM::VIBRATION:
      result->resize(2,0);
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

  result->redim(setIntegratedSol ? nen+1 : nen);
  return result;
}


bool Poisson::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  // Conductivity scaled by integration point weight at this point
  double cw = this->getMaterial(X)*fe.detJxW;

  if (!elMat.A.empty()) {
    Matrix Ab;
    Matrix& A = setIntegratedSol ? Ab : elMat.A.front();
    // Integrate the coefficient matrix // EK += kappa * dNdX * dNdX^T * |J|*w
    A.multiply(fe.dNdX,fe.dNdX,false,true,!setIntegratedSol,cw);
    if (this->setIntegratedSol) {
      elMat.A.front().addBlock(Ab, 1.0, 1, 1);
      const size_t nrow = elMat.A.front().rows();
      for (size_t col = 1; col <= fe.N.size(); ++col)
        elMat.A.front()(nrow, col) = elMat.A.front()(col, nrow) += fe.detJxW*fe.N(col);
    }
    if (m_mode == SIM::VIBRATION)
      elMat.A[1].outer_product(fe.N, fe.N, true, fe.detJxW);
  }

  // Lambda function for integration of the internal force vector
  auto&& evalIntForce = [cw,fe](Vector& S, const Vector& eV)
  {
    RealArray tmp;
    // S -= dNdX * (dNdX^t * eV) * cw
    return fe.dNdX.multiply(eV,tmp,true) && fe.dNdX.multiply(tmp,S,-cw,1.0);
  };

  if (!elMat.b.empty())
  {
    // Integrate heat source, if defined
    if (heatSrc)
      elMat.b.front().add(fe.N,(*heatSrc)(X)*fe.detJxW); // EV += N*h(x)*|J|*w

    if (m_mode == SIM::RHS_ONLY && !elmInt.vec.empty())
      // Integrate the internal forces based on current solution
      if (!evalIntForce(elMat.b.front(),elmInt.vec.front()))
        return false;
  }

  if (dualRHS && elmInt.vec.size() > 1)
  {
    if (elmInt.vec[1].empty())
      return true; // the extraction function is zero in this element

    // Integrate the dual load vector
    return evalIntForce(elMat.b[1],elmInt.vec[1]);
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
  return vtf->writeVectors(fluxVal,geoBlk,++nBlock,"Heat flux",iStep);
}


/*!
  This method is overridden to account for global %Lagrange multipliers.
*/

bool Poisson::evalSol (Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const
{
  return this->evalSol1(s, fe, X, MNPC, MNPC.size() > fe.N.size());
}


bool Poisson::evalSol2 (Vector& q, const Vectors& eV,
                        const FiniteElement& fe, const Vec3& X) const
{
  // Evaluate the heat flux vector, q = -kappa*du/dX = -kappa*dNdX^T*eV
  if (eV.empty() || !fe.dNdX.multiply(eV.front(),q,true))
  {
    std::cerr <<" *** Poisson::evalSol: Invalid solution vector."<< std::endl;
    return false;
  }

  const double kappa = this->getMaterial(X);
  q *= -kappa;
  const Vec3 ex = aSol ? (*aSol->getScalarSecSol())(X) : Vec3();
  for (size_t i = 0; i < nsd; ++i)
    q.push_back(ex[i]);
  return true;
}


Vector* Poisson::getExtractionField (size_t ifield)
{
  return dualRHS && ifield < primsol.size() ? &primsol[ifield] : nullptr;
}


std::string Poisson::getField1Name (size_t, const char* prefix) const
{
  if (!prefix) return "u";

  return prefix + std::string(" u");
}


std::string Poisson::getField2Name (size_t i, const char* prefix) const
{
  if (i >= this->getNoFields(2)) return "";

  if (nsd == 2 && i >= 2)
    ++i;

  static const char* s[6] = { "q_x","q_y","q_z",
                             "(q_ex)_x", "(q_ex)_y", "(q_ex)_z" };

  return prefix ? prefix + std::string(" ") + s[i] : s[i];
}


bool Poisson::suppressOutput (size_t j, ASM::ResultClass resultClass) const
{
  if (resultClass == ASM::SECONDARY)
    return j >= nsd;

  return false;
}


double Poisson::getMaterial (const Vec3& X) const
{
  return kappaF ? (*kappaF)(X) : kappaC;
}


/*!
  \note The Integrand object is allocated dynamically and has to be deleted
  manually when leaving the scope of the pointer variable receiving the
  returned pointer value.
*/

NormBase* Poisson::getNormIntegrand (AnaSol* asol) const
{
  return new PoissonNorm(*const_cast<Poisson*>(this), normIntegrandType,
                         asol ? asol->getScalarSecSol() : nullptr);
}


void Poisson::clearGalerkinProjections ()
{
  for (VecFunc* f : galerkin)
    delete f;
  galerkin.clear();
}


PoissonNorm::PoissonNorm (Poisson& p, int integrandType, VecFunc* a)
  : NormBase(p), anasol(a), integrdType(integrandType)
{
  nrcmp = myProblem.getNoFields(2);
  projBou = true;
}


PoissonNorm::~PoissonNorm() = default;


bool PoissonNorm::evalInt (LocalIntegral& elmInt, const FiniteElement& fe,
                           const Vec3& X) const
{
  const Poisson& problem = static_cast<const Poisson&>(myProblem);
  const size_t nsol = std::min(1+problem.numExtrFunction(), elmInt.vec.size());

  // Evaluate the finite element and dual solution gradient fields
  RealArray tmp;
  Vec3Vec gradh;
  gradh.reserve(nsol);
  for (size_t i = 0; i < nsol; i++)
    if (elmInt.vec[i].empty())
      gradh.push_back(Vec3());
    else if (fe.dNdX.multiply(elmInt.vec[i],tmp,true))
      gradh.push_back(tmp);
    else
      std::cerr <<" *** PoissonNorm::evalInt: Invalid solution vector "
                << i+1 << std::endl;
  if (gradh.size() < nsol)
    return false;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);
  const Vector& vh = elmInt.vec.front();
  const Vec3& epsh = gradh.front();

  // Evaluate the temperature field
  double u = vh.dot(fe.N);
  // Evaluate the heat source field
  double h = problem.getHeat(X);
  // Evaluate the heat conductivity
  double kappa = problem.getMaterial(X);
  // Evaluate the inverse conductivity scaled by the integration point weight
  double cwInv = fe.detJxW / kappa;

  // Integrate the energy norm a(u^h,u^h)
  pnorm[0] += kappa*epsh*epsh*fe.detJxW;
  // Integrate the external energy (h,u^h)
  if (problem.extEner)
    pnorm[1] += h*u*fe.detJxW;

  Vec3 sigma;
  size_t ip = 2;
  if (anasol)
  {
    // Evaluate the analytical heat flux
    sigma = (*anasol)(X);
    // Integrate the energy norm a(u,u)
    pnorm[ip++] += sigma*sigma*cwInv;
    // Integrate the error in energy norm a(u-u^h,u-u^h)
    Vec3 error = sigma + kappa*epsh;
    pnorm[ip++] += error*error*cwInv;
  }

  // Integrate the volume
  pnorm[ip++] += fe.detJxW;

  for (size_t j = 1; j < nsol; j++)
  {
    // Evaluate the variational-consistent postprocessing quantity, a(u^h,w)
    pnorm[ip++] -= kappa*epsh*gradh[j]*fe.detJxW;
    if (anasol) // Evaluate the corresponding exact quantity, a(u,w)
      pnorm[ip++] += sigma*gradh[j]*fe.detJxW;
  }

#if INT_DEBUG > 3
  std::cout <<"\nPoissonNorm::evalInt(X = "<< X <<")";
  if (anasol) std::cout <<"\nsigma ="<< sigma;
  std::cout <<"epsh ="<< epsh;
  for (size_t i = 1; i < nsol; i++)
    std::cout <<"epsz("<< i <<") ="<< gradh[i] <<"a(u^h,w"<< i
              <<"): "<< pnorm[ip+2*(i-nsol)] << std::endl;
#endif

  size_t pi = 0;
  const size_t nsd = problem.getNoSpaceDim();
  for (const Vector& psol : pnorm.psol)
  {
    if (!psol.empty() || (projFields.size() > pi && projFields[pi]))
    {
      // Evaluate projected heat flux field and projected heat flux gradient
      Vec3 sigmar, error, qex;
      Matrix dSigmadX;
      if (!projFields.empty() && projFields[pi])
      {
        Vector projSol(nrcmp);
        projFields[pi]->valueFE(fe, projSol);
        Matrix dS;
        projFields[pi]->gradFE(fe, dS);
        dSigmadX.resize(nsd, nsd);
        dS.extractBlock(dSigmadX, 1, 1);
        sigmar = Vec3(projSol.ptr(), nsd);
        qex = Vec3(projSol.ptr() + nsd, nsd);
      }
      else
      {
        for (size_t j = 0; j < nsd; j++)
          sigmar[j] = psol.dot(fe.N,j,nrcmp);
        for (size_t j = 0; j < nsd; j++)
          qex[j] = psol.dot(fe.N,j+nsd,nrcmp);
        // Evaluate the projected heat flux gradient.
        // Notice that the matrix multiplication method used here treats
        // the element vector (psol) as a matrix whose number of columns
        // equals the number of rows in the matrix fe.dNdX.
        Matrix dS;
        if (!dS.multiplyMat(psol,fe.dNdX)) // dSigmadX = psol*dNdX
          return false;
        dSigmadX.resize(nsd, nsd);
        dS.extractBlock(dSigmadX, 1, 1);
      }

      // Integrate the energy norm a(u^r,u^r)
      pnorm[ip++] += sigmar*sigmar*cwInv;
      // Integrate the estimated error in energy norm a(u^r-u^h,u^r-u^h)
      error = sigmar + kappa*epsh;
      pnorm[ip++] += error*error*cwInv;

      // Evaluate the interior residual of the projected solution
      double Res = h - dSigmadX.trace();
      // Integrate the residual error in the projected solution
      pnorm[ip++] += fe.h*fe.h*Res*Res*fe.detJxW;

      if (anasol)
      {
        // Integrate the error in the projected solution a(u-u^r,u-u^r)
        error = sigma - sigmar;
        pnorm[ip++] += error*error*cwInv;
        // Integrate the error in projected solution norm a(u_ex^r-u^h,u_ex^r-u^h)
        error = sigma - qex;
        pnorm[ip++] += error*error*cwInv;
        ip += 2; // Make room for the local effectivity indices here
      }
    }
    else if (integrdType & SECOND_DERIVATIVES)
    {
      // Integrate the residual error in the FE solution
      double Res = h;
      for (size_t j = 1; j <= fe.N.size(); j++)
        Res += fe.d2NdX2.trace(j)*vh(j);

      pnorm[ip+1] += fe.h*fe.h*Res*Res*fe.detJxW;
      ip += this->getNoFields(pi+2); // Dummy entries in order to get norm in right place
    }
    ++pi;
  }

  return true;
}


bool PoissonNorm::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                           const Vec3& X, const Vec3& normal) const
{
  const Poisson& problem = static_cast<const Poisson&>(myProblem);
  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate the surface heat flux
  double h = -problem.getFlux(X,normal);
  // Evaluate the temperature field
  double u = elmInt.vec.front().dot(fe.N);

  // Integrate the external energy (h,u^h)
  if (problem.extEner)
    pnorm[1] += h*u*fe.detJxW;

  const size_t nsd = problem.getNoSpaceDim();
  size_t ip = this->getNoFields(1) + 2;
  size_t pi = 0;
  for (const Vector& psol : pnorm.psol)
  {
    if (!psol.empty() || (projFields.size() > pi && projFields[pi]))
    {
      // Evaluate projected heat flux field
      Vec3 sigmar;
      if (!projFields.empty() && projFields[pi])
      {
        Vector projSol(nrcmp);
        projFields[pi]->valueFE(fe,projSol);
        sigmar = Vec3(projSol.ptr(), nsd);
      }
      else
        for (size_t j = 0; j < nsd; j++)
          sigmar[j] = psol.dot(fe.N,j,nrcmp);

      // Evaluate the boundary jump term
      double Jump = h + sigmar*normal;
      // Integrate the residual error in the projected solution
      pnorm[ip] += 0.5*fe.h*Jump*Jump*fe.detJxW;
      ip += problem.getNoFields(pi+2);
    }
    else if (integrdType & SECOND_DERIVATIVES)
      ip += problem.getNoFields(pi+2); // TODO: Add residual jump terms?
    ++pi;
  }

  return true;
}


bool PoissonNorm::finalizeElement (LocalIntegral& elmInt)
{
  if (!anasol) return true;

  ElmNorm& pnorm = static_cast<ElmNorm&>(elmInt);

  // Evaluate local effectivity indices as a(e^r,e^r)/a(e,e)
  // with e^r = u^r - u^h  and  e = u - u^h
  for (size_t ip = this->getNoFields(1), g = 2; ip < pnorm.size(); ip += this->getNoFields(g++))
  {
    pnorm[ip+5] =  pnorm[ip] / pnorm[3];
    pnorm[ip+6] = (pnorm[ip]+pnorm[ip+1]) / pnorm[3];
  }

  return true;
}


size_t PoissonNorm::getNoFields (int group) const
{
  if (group < 1)
    return this->NormBase::getNoFields();
  else if (group > 1)
    return anasol ? 7 : 3;

  size_t nExt = static_cast<Poisson&>(myProblem).numExtrFunction();
  return anasol ? 5 + 2*nExt : 3 + nExt;
}


std::string PoissonNorm::getName (size_t i, size_t j, const char* prefix) const
{
  const size_t nx = i > 1 ? 0 : static_cast<Poisson&>(myProblem).numExtrFunction() * 2;

  static const auto primary = std::array{
    "a(u^h,u^h)^0.5",
    "(h,u^h)^0.5",
    "a(u,u)^0.5",
    "a(e,e)^0.5, e=u-u^h",
    "volume",
  };

  static const auto extract = std::array{
    "a(u^h,w)",
    "a(u,w)",
  };

  const auto proj = std::array{
    "a(u^r,u^r)^0.5",
    (integrdType & SECOND_DERIVATIVES ? "|u^h|_res" : "a(e,e)^0.5, e=u^r-u^h"),
    "res(u^r)^0.5",
    "a(e,e)^0.5, e=u-u^r",
    "a(e,e)^0.5, e=u-u_ex^r",
    "effectivity index^*",
    "effectivity index^RES",
  };

  static const auto dual = std::array{
    "a(z^h,z^h)^0.5",
    "(E(u)*E(z))^0.5, E(v)=a(e,e), e=v^r-v^h"
  };


  if (i == 0 || j == 0 || j > primary.size() + extract.size() + nx)
    return this->NormBase::getName(i,j,prefix);

  if (i >= 2 && j > proj.size())
    return this->NormBase::getName(i,j,prefix);

  if (i == 1)
  {
    if (!anasol && j > 2 && j < 4)
      j += 2;

    if (j > primary.size() && nx > 0)
    {
      j -= primary.size();
      char comp[32];
      if (j < 3)
        return extract[j-1];
      else if (!anasol)
        sprintf(comp,"a(u^h,w%zu)",j);
      else if (j%2)
        sprintf(comp,"a(u^h,w%zu)",j/2+1);
      else
        sprintf(comp,"a(u,w%zu)",j/2);
      return comp;
    }
    else if (prefix && !strncmp(prefix,"Dual",4))
    {
      if (j <= dual.size())
        return dual[j-1];
      else
        return "";
    }
  }

  std::string s(prefix ? prefix + std::string(" ") : "");
  s += (i == 1 ? primary[j-1] : proj[j-1]);
  return s;
}


bool PoissonNorm::hasElementContributions (size_t i, size_t j) const
{
  if (integrdType & SECOND_DERIVATIVES && i == 2)
    return j == 2 || j == 7;

  return i > 1 || j != 2;
}


int PoissonNorm::getIntegrandType () const
{
  return integrdType | myProblem.getIntegrandType();
}


void PoissonNorm::setProjectedFields(Fields* field, size_t idx)
{
  if (idx >= projFields.size())
    projFields.resize(idx+1);
  projFields[idx].reset(field);
}


Poisson::Robin::Robin (unsigned short int n, const Poisson& itg) :
  IntegrandBase(n), integrand(itg)
{
  alpha = nullptr;
  g = nullptr;
}


bool Poisson::Robin::evalBou (LocalIntegral& elmInt, const FiniteElement& fe,
                              const Vec3& X, const Vec3& normal) const
{
  ElmMats& elMat = static_cast<ElmMats&>(elmInt);

  Vec3 ax = alpha ? (*alpha)(X) : Vec3(1.0, 1.0, 1.0);
  if (g) ax.y = (*g)(X);

  elMat.A.front().outer_product(fe.N, fe.N, true, fe.detJxW * ax.x); // mass
  elMat.b.front().add(fe.N, fe.detJxW * ax.y); // source

  return true;
}
