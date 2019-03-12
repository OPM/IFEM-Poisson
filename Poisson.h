// $Id$
//==============================================================================
//!
//! \file Poisson.h
//!
//! \date Apr 16 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Integrand implementations for Poisson problems.
//!
//==============================================================================

#ifndef _POISSON_H
#define _POISSON_H

#include "IntegrandBase.h"
#include "GlobalIntegral.h"
#include "Vec3.h"

class RealFunc;
class VecFunc;


/*!
  \brief Class representing the integrand of the Poisson problem.

  \details This class supports constant isotropic conductivity only.

  See the document doc/Integrands/Poisson.pdf for the theoretical foundation
  of the integrand implemented in the evalInt() and evalBou() methods.
*/

class Poisson : public IntegrandBase
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  explicit Poisson(unsigned short int n = 3);
  //! \brief The destructor deletes the functions to be Galerkin-projected.
  virtual ~Poisson() { this->clearGalerkinProjections(); }

  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const TiXmlElement* elem);

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { tracFld = tf; }
  //! \brief Defines the heat flux field to use in Neumann boundary conditions.
  void setTraction(RealFunc* ff) { fluxFld = ff; }
  //! \brief Defines the heat source field.
  void setSource(RealFunc* src) { heatSrc = src; }

  //! \brief Defines the conductivity.
  void setMaterial(double K) { kappa = K; }
  //! \brief Returns the conductivity.
  double getMaterial() const { return kappa; }

  //! \brief Returns the number of Galerkin projections.
  size_t getNoGalerkin() const { return galerkin.size(); }
  //! \brief Clears up the Galerkin projections.
  void clearGalerkinProjections();

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  using IntegrandBase::initIntegration;
  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  //! \brief Defines the global integral for calculating reaction forces only.
  void setReactionIntegral(GlobalIntegral* gq) { delete reacInt; reacInt = gq; }
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using IntegrandBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  using IntegrandBase::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s, const FiniteElement& fe,
                       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE solution values at current point
  //! \param[in] eV Element solution vector
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol(Vector& s, const Vector& eV,
               const Matrix& dNdX, const Vec3& X) const;

  //! \brief Evaluates the boundary heat flux (if any) at specified point.
  double getFlux(const Vec3& X, const Vec3& n) const;
  //! \brief Evaluates the heat source (if any) at specified point.
  double getHeat(const Vec3& X) const;

  //! \brief Writes the heat flux vector for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the heat flux vectors
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& geoBlk, int& nBlock) const;

  //! \brief Returns whether there are any heat flux values to write to VTF.
  virtual bool hasTractionValues() const { return !fluxVal.empty(); }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \param[in] asol Pointer to analytical solution (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = nullptr) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const { return fld > 1 ? nsd : 1; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual std::string getField1Name(size_t, const char* prefix) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix) const;

  //! \brief Defines the properties of the resulting linear system.
  virtual LinAlg::LinearSystemType getLinearSystemType() const
  {
    return LinAlg::SPD;
  }

private:
  // Physical properties (constant)
  double kappa; //!< Conductivity

  VecFunc*  tracFld; //!< Pointer to boundary traction field
  RealFunc* fluxFld; //!< Pointer to boundary normal flux field
  RealFunc* heatSrc; //!< Pointer to interior heat source

  GlobalIntegral* reacInt; //!< Reaction-forces-only integral

  int normIntegrandType; //!< Integrand type for norm class

  mutable std::vector<Vec3Pair> fluxVal; //!< Heat flux point values

  std::vector<VecFunc*> galerkin; //!< Functions to be Galerkin-projected

public:
  char extEner; //!< If \e true, external energy is to be computed
};


/*!
  \brief Class representing the integrand of Poisson energy norms.
*/

class PoissonNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] integrandType Integrand type flag
  //! \param[in] a The analytical heat flux (optional)
  PoissonNorm(Poisson& p, int integrandType, VecFunc* a = nullptr);
  //! \brief Empty destructor.
  virtual ~PoissonNorm() {}

  using NormBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using NormBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X, const Vec3& normal) const;

  using NormBase::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to compute effectivity indices.
  //! \param elmInt The local integral object to receive the contributions
  virtual bool finalizeElement(LocalIntegral& elmInt);

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual std::string getName(size_t i, size_t j, const char* prefix) const;

  //! \brief Returns whether a norm quantity stores element contributions.
  virtual bool hasElementContributions(size_t i, size_t j) const;

  //! \brief Defines which FE quantities are needed by the integrand.
  virtual int getIntegrandType() const { return integrandType; }

private:
  VecFunc* anasol; //!< Analytical heat flux
  int integrandType; //!< Integrand type
};

#endif
