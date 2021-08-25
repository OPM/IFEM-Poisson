// $Id$
//==============================================================================
//!
//! \file SIMPoisson.h
//!
//! \date May 25 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#ifndef _SIM_POISSON_H
#define _SIM_POISSON_H

#include "SIMMultiPatchModelGen.h"
#include "SIMbase.h"
#include "Poisson.h"
#include "TextureProperties.h"


class DataExporter;
class TimeStep;


/*!
  \brief Driver class for NURBS-based FEM analysis of Poisson problems.
*/

template<class Dim> class SIMPoisson : public SIMMultiPatchModelGen<Dim>
{
public:
  //! \brief Default constructor.
  explicit SIMPoisson(bool checkRHS = false, bool ds = false);

  //! \brief The destructor zero out the integrand pointer (deleted by parent).
  virtual ~SIMPoisson();

  //! \brief Initializes the property containers of the model.
  void clearProperties() override;

  //! \brief Returns the number of right-hand-side vectors.
  size_t getNoRHS() const override;

  //! \brief Returns whether a dual solution is available or not.
  bool haveDualSol() const override;

  //! \brief Registers data fields for output.
  void registerFields(DataExporter& exporter);

  //! \brief Sets the solution vector for output.
  void setSol(const Vector* sol) { solution = sol; }

  //! \brief Toggles writing the load vector to VTF.
  void setVizRHS(bool viz) { vizRHS = viz; }

  //! \brief Sets the ASCII file name prefix.
  void setASCIIfile(const char* filename);

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  //!
  //! \details This method is not used in adaptive simulations.
  //! It also writes out the boundary tractions (if any) and the Dirichlet
  //! boundary conditions, as this data is regarded part of the model
  //! and not as simulation results.
  bool saveModel(char* fileName, int& geoBlk, int& nBlock);

  //! \brief Assembles and solves the linear system.
  bool solveStep(TimeStep&);

  //! \brief Saves solution-dependent quantities to file for postprocessing.
  bool saveStep(TimeStep&, int& nBlock);

  using SIMMultiPatchModelGen<Dim>::solveSystem;
  //! \brief Solves the assembled linear system of equations for a given load.
  //! \param[out] solution Global primary solution vector
  //! \param[in] printSol Print solution if its size is less than \a printSol
  //! \param[out] rCond Reciprocal condition number
  //! \param[in] compName Solution name to be used in norm output
  //! \param[in] newLHS If \e false, reuse the LHS-matrix from previous call
  //! \param[in] idxRHS Index to the right-hand-side vector to solve for
  //!
  //! This overloaded version also computes the reaction forces along a given
  //! boundary. This requires an additional assembly loop calculating the
  //! internal forces only, since we only are doing a linear solve here.
  bool solveSystem(Vector& solution, int printSol, double* rCond,
                   const char* compName, bool newLHS, size_t idxRHS) override;

  //! \brief Returns current reaction force vector.
  const Vector* getReactionForces() const override
  {
    return myReact.empty() ? nullptr : &myReact;
  }

  //! \brief Prints a norm group to the log stream.
  //! \param[in] gNorm The global norm values
  //! \param[in] fNorm Global reference norm values
  //! \param[in] name Name of norm group
  void printNormGroup(const Vector& gNorm, const Vector& fNorm,
                      const std::string& name) const override;

  //! \brief Returns the name of this simulator.
  //! \details This method is typically reimplemented in sub-classes that are
  //! parts of a partitioned solution method and are used to identify the basis
  //! for the result fields associated with each simulator in the HDF5 output.
  std::string getName() const override { return "Poisson"; }

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to resolve inhomogeneous boundary
  //! condition fields in case they are derived from the analytical solution.
  void preprocessA() override;

  //! \brief Performs some pre-processing tasks on the FE model.
  bool preprocessB() override;

  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  bool parse(char* keyWord, std::istream& is) override;

  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  bool parse(const TiXmlElement* elem) override;

private:
  //! \brief Parses a dimension-specific data section from an input stream.
  //! \details This function allows for specialization of the template while
  //! still reusing as much code as possible. Only for dimension-specific code.
  bool parseDimSpecific(char* keyWord, std::istream& is);
  //! \brief Parses a dimension-specific data section from the an XML element.
  //! \details This function allows for specialization of the template while
  //! still reusing as much code as possible. Only for dimension-specific code.
  bool parseDimSpecific(const TiXmlElement* child);

protected:
  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  bool initMaterial(size_t propInd) override;

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  bool initNeumann(size_t propInd) override;

  //! \brief Returns norm index of the integrated volume.
  size_t getVolumeIndex() const override
  {
    return this->haveAnaSol() ? 5 : 3;
  }

  //! \brief Reverts the square-root operation on the volume and VCP quantities.
  bool postProcessNorms(Vectors& gNorm, Matrix* eNorm) override
  {
    return this->revertSqrt(gNorm,eNorm);
  }

private:
  Poisson   prob;     //!< Data and methods for the Poisson problem
  Poisson::Robin robinBC; //!< Integrand for Robin conditions

  //! \brief Struct with either a constant or a function value for kappa.
  struct Kappa {
    double constant; //!< Constant value
    std::shared_ptr<RealFunc> func; //!< Function value
  };

  std::vector<Kappa> mVec; //!< Kappa properties
  int       aCode[2]; //!< Analytical BC code (used by destructor)
  TextureProperties tprops; //!< Texture property (for kappa)

  Vector    myLoad;   //!< External load vector (for VTF export)
  Vector    mySolVec; //!< Primary solution vector
  Vector    myReact;  //!< Nodal reaction forces
  Vectors   myProj;   //!< Projected solution vectors
  Matrix    myNorm;   //!< Element norms

  std::vector<Mode> modes; //!< Eigen modes

  const Vector* solution; //!< Pointer to primary solution vector

  bool dualS; //!< If \e true, also solve the dual problem

  bool        vizRHS;    //!< If \e true, store load vector to VTF
  std::string asciiFile; //!< ASCII output file prefix
};

#endif
