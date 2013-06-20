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
#include "Vec3.h"


/*!
  \brief Class representing the integrand of the Poisson problem.
  \details This class supports constant isotropic constitutive properties only.
  Properties with spatial variation has to be implemented as sub-classes.
*/

class Poisson : public IntegrandBase
{
public:
  //! \brief The default constructor initializes all pointers to zero.
  //! \param[in] n Number of spatial dimensions
  Poisson(unsigned short int n = 3);
  //! \brief Nearly empty destructor.
  virtual ~Poisson() { clearGalerkinProjections(); }

  //! \brief Defines the traction field to use in Neumann boundary conditions.
  void setTraction(VecFunc* tf) { tracFld = tf; }
  //! \brief Defines the heat flux field to use in Neumann boundary conditions.
  void setTraction(RealFunc* ff) { fluxFld = ff; }
  //! \brief Defines the heat source field.
  void setSource(RealFunc* src) { heatSrc = src; }

  //! \brief Defines the conductivity (constitutive property).
  void setMaterial(double K) { kappa = K; }

  //! \brief Add a Galerkin projection
  void addGalerkin(VecFunc* g) { galerkin.push_back(g); }

  //! \brief Returns the number of Galerkin projections
  size_t getNoGalerkin() const { return galerkin.size(); }

  //! \brief Clear up Galerkin projections
  void clearGalerkinProjections();

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  //! \param[in] nBp Total number of boundary integration points
  virtual void initIntegration(size_t nGp, size_t nBp);

  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  //! \param[in] neumann Whether or not we are assembling Neumann BC's
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t,
                                          bool neumann) const;

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s Array of solution field values at current point
  //! \param[in] N Basis function values at current point
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] MNPC Nodal point correspondance for the basis function values
  virtual bool evalSol(Vector& s,
		       const Vector& N, const Matrix& dNdX,
		       const Vec3& X, const std::vector<int>& MNPC) const;

  //! \brief Evaluates the finite element (FE) solution at an integration point.
  //! \param[out] s The FE solution values at current point
  //! \param[in] eV Element solution vector
  //! \param[in] dNdX Basis function gradients at current point
  //! \param[in] X Cartesian coordinates of current point
  bool evalSol(Vector& s, const Vector& eV,
               const Matrix& dNdX, const Vec3& X) const;

  //! \brief Evaluates the analytical solution at an integration point.
  //! \param[out] s The analytical solution values at current point
  //! \param[in] asol The analytical solution field (heat flux)
  //! \param[in] X Cartesian coordinates of current point
  virtual bool evalSol(Vector& s, const VecFunc& asol, const Vec3& X) const;

  //! \brief Evaluates the boundary heat flux (if any) at specified point.
  double getFlux(const Vec3& X, const Vec3& n) const;
  //! \brief Evaluates the heat source (if any) at specified point.
  double getHeat(const Vec3& X) const;

  //! \brief Writes the heat flux vector for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the heat flux vectors
  //! \param[in] iStep Load/time step identifier
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& nBlock) const;

  //! \brief Returns whether there are any heat flux values to write to VTF.
  virtual bool hasTractionValues() const { return !fluxVal.empty(); }

  //! \brief Returns a pointer to an Integrand for solution norm evaluation.
  //! \note The Integrand object is allocated dynamically and has to be deleted
  //! manually when leaving the scope of the pointer variable receiving the
  //! returned pointer value.
  //! \param[in] asol Pointer to analytical solution (optional)
  virtual NormBase* getNormIntegrand(AnaSol* asol = 0) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld = 2) const { return fld > 1 ? nsd : 1; }
  //! \brief Returns the name of the primary solution field.
  //! \param[in] prefix Name prefix
  virtual const char* getField1Name(size_t, const char* prefix = 0) const;
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual const char* getField2Name(size_t i, const char* prefix = 0) const;

  //! \brief Sets up the constitutive matrix at current point.
  //! \param[out] C \f$ nsd\times nsd\f$-matrix
  //! \param[in] X Cartesian coordinates of current point
  //! \param[in] invers If \e true, set up the inverse matrix instead
  virtual bool formCmatrix(Matrix& C, const Vec3& X, bool invers = false) const;

private:
  // Physical properties (constant)
  double kappa; //!< Conductivity

protected:
  VecFunc*  tracFld; //!< Pointer to boundary traction field
  RealFunc* fluxFld; //!< Pointer to boundary normal flux field
  RealFunc* heatSrc; //!< Pointer to interior heat source

  mutable std::vector<Vec3Pair> fluxVal; //!< Heat flux point values

  std::vector<VecFunc*> galerkin; //!< Pointer to functions to be galerkin projected

  unsigned short int nsd; //!< Number of space dimensions (1, 2 or, 3)
};


/*!
  \brief Class representing the integrand of Poisson energy norms.
*/

class PoissonNorm : public NormBase
{
public:
  //! \brief The only constructor initializes its data members.
  //! \param[in] p The Poisson problem to evaluate norms for
  //! \param[in] a The analytical heat flux (optional)
  PoissonNorm(Poisson& p, VecFunc* a = 0);
  //! \brief Empty destructor.
  virtual ~PoissonNorm() {}

  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
		       const Vec3& X) const;

  //! \brief Evaluates the integrand at a boundary point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite Element quantities
  //! \param[in] X Cartesian coordinates of current integration point
  //! \param[in] normal Boundary normal vector at current integration point
  virtual bool evalBou(LocalIntegral& elmInt,  const FiniteElement& fe,
		       const Vec3& X, const Vec3& normal) const;

  //! \brief Finalizes the element norms after the numerical integration.
  //! \details This method is used to compute effectivity indices.
  //! \param elmInt The local integral object to receive the contributions
  virtual bool finalizeElement(LocalIntegral& elmInt, const TimeDomain&,size_t);

  //! \brief Returns whether this norm has explicit boundary contributions.
  virtual bool hasBoundaryTerms() const { return true; }

  //! \brief Adds external energy terms to relevant norms.
  //! \param gNorm Global norm quantities
  //! \param[in] energy Global external energy
  virtual void addBoundaryTerms(Vectors& gNorm, double energy) const;

  //! \brief Returns the number of norm groups or size of a specified group.
  //! \param[in] group The norm group to return the size of
  //! (if zero, return the number of groups)
  virtual size_t getNoFields(int group = 0) const;

  //! \brief Returns the name of a norm quantity.
  //! \param[in] i The norm group (one-based index)
  //! \param[in] j The norm number (one-based index)
  //! \param[in] prefix Common prefix for all norm names
  virtual const char* getName(size_t i, size_t j, const char* prefix) const;

  //! \brief Returns whether a norm quantity stores element contributions.
  virtual bool hasElementContributions(size_t i, size_t j) const;

private:
  VecFunc* anasol; //!< Analytical heat flux
};

#endif
