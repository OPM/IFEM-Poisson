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

#include "GlobalIntegral.h"
#include "Integrand.h"
#include "IntegrandBase.h"
#include "MatVec.h"
#include "LinAlgenums.h"
#include "SIMenums.h"
#include "Vec3.h"

#include <string>
#include <vector>


class AnaSol;
class FiniteElement;
class FunctionBase;
class LocalIntegral;
class RealFunc;
class VecFunc;
class VTF;


/*!
  \brief Class representing the integrand of the Poisson problem.

  \details This class supports constant isotropic conductivity only.

  See the document doc/Integrands/Poisson.pdf for the theoretical foundation
  of the integrand implemented in the evalInt() and evalBou() methods.
*/

class Poisson : public IntegrandBase
{
public:
  //! \brief Class representing the Robin boundary conditions.
  class Robin : public IntegrandBase
  {
  public:
    //! \brief The constructor forwards to the parent class constructor.
    //! \param[in] n Number of spatial dimensions
    //! \param[in] itg Main integrand instance
    Robin(unsigned short int n, const Poisson& itg);
    //! \brief Empty destructor.
    virtual ~Robin() {}

    //! \brief Returns that this integrand has no interior contributions.
    bool hasInteriorTerms() const override { return false; }

    using IntegrandBase::getLocalIntegral;
    //! \brief Returns a local integral contribution object for given element.
    //! \param[in] nen Number of nodes on element
    //! \param[in] iEl Element number
    LocalIntegral* getLocalIntegral(size_t nen, size_t iEl, bool) const override
    {
      return integrand.getLocalIntegral(nen, iEl, false);
    }

    using IntegrandBase::evalBou;
    //! \brief Evaluates the integrand at a boundary point.
    //! \param elmInt The local integral object to receive the contributions
    //! \param[in] fe Finite element data of current integration point
    //! \param[in] X Cartesian coordinates of current integration point
    //! \param[in] normal Boundary normal vector at current integration point
    bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
                 const Vec3& X, const Vec3& normal) const override;

    //! \brief Sets the coefficient function \a alpha to \a f.
    void setAlpha(const VecFunc* f) { alpha = f; g = nullptr; }
    //! \brief Sets the flux function \a g to \a f.
    void setFlux(const RealFunc* f) { alpha = nullptr; g = f; }

  protected:
    const VecFunc* alpha = nullptr; //!< Coefficient - alpha(1) * u + du/dn = alpha(2)
    const RealFunc* g = nullptr; //!< Coefficient - u + du/dn = g
    const Poisson& integrand; //!< Main integrand instance
  };

  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  explicit Poisson(unsigned short int n = 3);
  //! \brief The destructor deletes the functions to be Galerkin-projected.
  virtual ~Poisson() { this->clearGalerkinProjections(); }

  //! \brief Parses a data section from an XML-element.
  bool parse(const TiXmlElement* elem) override;

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { tracFld = tf; }
  //! \brief Defines the heat flux field to use in Neumann boundary conditions.
  void setTraction(RealFunc* ff) { fluxFld = ff; }
  //! \brief Defines the heat source field.
  void setSource(RealFunc* src) { heatSrc = src; }
  //! \brief Defines the extraction function of the dual problem.
  void setDualRHS(FunctionBase* df) { dualRHS = df; }
  //! \brief Defines an extraction function for VCP.
  void addExtrFunction(FunctionBase* exf);
  //! \brief Returns the number of extraction functions.
  size_t numExtrFunction() const { return dualFld.size(); }

  //! \brief Defines the conductivity.
  void setMaterial(double K) { kappa = K; }
  //! \brief Defines the conductivity.
  void setMaterial(RealFunc* K) { kappaF = K; }
  //! \brief Evaluates the conductivity at specified point.
  double getMaterial(const Vec3& = Vec3()) const;

  //! \brief Returns the number of Galerkin projections.
  size_t getNoGalerkin() const { return galerkin.size(); }
  //! \brief Clears up the Galerkin projections.
  void clearGalerkinProjections();

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  void setMode(SIM::SolutionMode mode) override;

  using IntegrandBase::initIntegration;
  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  void initIntegration(size_t nGp, size_t nBp) override;

  using IntegrandBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] X0 Cartesian coordinates of the element center
  //! \param elmInt Local integral for element
  bool initElement(const std::vector<int>& MNPC, const FiniteElement&,
                   const Vec3& X0, size_t, LocalIntegral& elmInt) override;

  //! \brief Defines the global integral for calculating reaction forces only.
  void setReactionIntegral(GlobalIntegral* gq) { delete reacInt; reacInt = gq; }
  //! \brief Returns the system quantity to be integrated by \a *this.
  GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const override;

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BCs
  LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                  bool neumann) const override;

  using IntegrandBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  using IntegrandBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE solution values at current point
  //! \param[in] eV Element solution vectors
  //! \param[in] fe Finite element data at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol2(Vector& s, const Vectors& eV,
                const FiniteElement& fe, const Vec3& X) const override;

  //! \brief Evaluates the boundary heat flux (if any) at specified point.
  double getFlux(const Vec3& X, const Vec3& n) const;
  //! \brief Evaluates the heat source (if any) at specified point.
  double getHeat(const Vec3& X) const;

  //! \brief Writes the heat flux vector for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the heat flux vectors
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  bool writeGlvT(VTF* vtf, int iStep, int& geoBlk, int& nBlock) const override;

  //! \brief Returns whether there are any heat flux values to write to VTF.
  bool hasTractionValues() const override { return !fluxVal.empty(); }

  //! \brief Returns the patch-wise extraction function field, if any.
  //! \param[in] ifield 1-based index of the field to return
  Vector* getExtractionField(size_t ifield) override;

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \param[in] asol Pointer to analytical solution (optional)
  NormBase* getNormIntegrand(AnaSol* asol) const override;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  size_t getNoFields(int fld) const override { return fld > 1 ? nsd : 1; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  std::string getField1Name(size_t, const char* prefix) const override;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  std::string getField2Name(size_t i, const char* prefix) const override;

  //! \brief Defines the properties of the resulting linear system.
  LinAlg::LinearSystemType getLinearSystemType() const override
  {
    return LinAlg::SPD;
  }

  //! \brief Defines which FE quantities are needed by the integrand.
  int getIntegrandType() const override
  {
    return dualFld.empty() ? STANDARD : ELEMENT_CENTER;
  }

private:
  // Physical properties
  double kappa; //!< Conductivity (constant)
  const RealFunc* kappaF = nullptr; //!< Conductivity as a function

  VecFunc*  tracFld; //!< Pointer to boundary traction field
  RealFunc* fluxFld; //!< Pointer to boundary normal flux field
  RealFunc* heatSrc; //!< Pointer to interior heat source

  FunctionBase*              dualRHS; //!< Extraction function for dual RHS
  std::vector<FunctionBase*> dualFld; //!< Extraction functions for VCP

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
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X) const override;

  using NormBase::evalBou;
  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3& X, const Vec3& normal) const override;

  using NormBase::finalizeElement;
  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to compute effectivity indices.
  //! \param elmInt The local integral object to receive the contributions
  bool finalizeElement(LocalIntegral& elmInt) override;

  //! \brief Returns whether this norm has explicit boundary contributions.
  bool hasBoundaryTerms() const override { return true; }

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  size_t getNoFields(int group) const override;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  std::string getName(size_t i, size_t j, const char* prefix) const override;

  //! \brief Returns whether a norm quantity stores element contributions.
  bool hasElementContributions(size_t i, size_t j) const override;

  //! \brief Defines which FE quantities are needed by the integrand.
  int getIntegrandType() const override;

private:
  VecFunc* anasol; //!< Analytical heat flux
  int integrdType; //!< Integrand type flag
};

#endif
