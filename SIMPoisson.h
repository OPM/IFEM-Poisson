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
#include "Poisson.h"
#include "ReactionsOnly.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "ASMbase.h"
#include "AnaSol.h"
#include "Utilities.h"
#include "DataExporter.h"
#include "IFEM.h"
#include "tinyxml.h"
#include <fstream>


/*!
  \brief Driver class for NURBS-based FEM analysis of Poisson problems.
*/

template<class Dim> class SIMPoisson : public SIMMultiPatchModelGen<Dim>
{
public:
  //! \brief Default constructor.
  explicit SIMPoisson(bool checkRHS = false)
    : SIMMultiPatchModelGen<Dim>(1,checkRHS),
      prob(Dim::dimension), solution(&mySolVec)
  {
    Dim::myProblem = &prob;
    aCode[0] = aCode[1] = 0;
    vizRHS = false;
  }

  //! \brief The destructor zero out the integrand pointer (deleted by parent).
  virtual ~SIMPoisson()
  {
    Dim::myProblem = nullptr; // Because it is not dynamically allocated
    Dim::myInts.clear();

    // To prevent the SIMbase destructor try to delete already deleted functions
    if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
    if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
  }

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties()
  {
    // To prevent SIMbase::clearProperties deleting the analytical solution
    if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
    if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
    aCode[0] = aCode[1] = 0;

    mVec.clear();
    prob.setSource(nullptr);
    prob.setTraction((RealFunc*)nullptr);
    prob.setTraction((VecFunc*)nullptr);
    prob.clearGalerkinProjections();
    this->Dim::clearProperties();
  }

  //! \brief Returns the number of right-hand-side vectors.
  virtual size_t getNoRHS() const { return 1 + prob.getNoGalerkin(); }

  //! \brief Registers data fields for output.
  void registerFields(DataExporter& exporter)
  {
    if (Dim::opt.eig > 0) {
      exporter.registerField("eigenmodes", "eigenmodes",
                             DataExporter::SIM, DataExporter::EIGENMODES);
      exporter.setFieldValue("eigenmodes", this, &modes);
    }
    else {
      int results = DataExporter::PRIMARY;
      if (!Dim::opt.pSolOnly)
        results |= DataExporter::SECONDARY;
      if (Dim::opt.saveNorms)
        results |= DataExporter::NORMS;

      exporter.registerField("u", "solution", DataExporter::SIM, results);
      exporter.setFieldValue("u", this, solution,
                             Dim::opt.project.empty() ? nullptr : &myProj,
                             (results & DataExporter::NORMS) ? &myNorm : nullptr);
    }
  }

  //! \brief Sets the solution vector for output.
  void setSol(const Vector* sol) { solution = sol; }

  //! \brief Toggles writing the load vector to VTF.
  void setVizRHS(bool viz) { vizRHS = viz; }

  //! \brief Sets the ASCII file name prefix.
  void setASCIIfile(const char* filename)
  {
    const char* end = strrchr(filename,'.');
    if (end)
      asciiFile.assign(filename,end);
    else
      asciiFile.assign(filename);
  }

  //! \brief Opens a new VTF-file and writes the model geometry to it.
  //! \param[in] fileName File name used to construct the VTF-file name from
  //! \param[out] geoBlk Running geometry block counter
  //! \param[out] nBlock Running result block counter
  //!
  //! \details This method is not used in adaptive simulations.
  //! It also writes out the boundary tractions (if any) and the Dirichlet
  //! boundary conditions, as this data is regarded part of the model
  //! and not as simulation results.
  bool saveModel(char* fileName, int& geoBlk, int& nBlock)
  {
    if (!asciiFile.empty())
    {
      // Write (possibly refined) model to g2-file
      std::ofstream osg(asciiFile+".g2");
      osg.precision(18);
      IFEM::cout <<"\nWriting updated g2-file "<< asciiFile+".g2" << std::endl;
      this->dumpGeometry(osg);
    }

    if (Dim::opt.format < 0)
      return true; // No VTF-output

    // Write geometry
    if (!this->writeGlvG(geoBlk,fileName))
      return false;

    // Write boundary tractions, if any
    if (!this->writeGlvT(1,geoBlk,nBlock))
      return false;

    // Write Dirichlet boundary conditions
    return this->writeGlvBC(nBlock);
  }

  //! \brief Assembles and solves the linear system.
  bool solveStep(TimeStep&)
  {
    if (!this->setMode(Dim::opt.eig == 0 ? SIM::STATIC : SIM::VIBRATION))
      return false;

    if (Dim::opt.eig == 0)
      this->initSystem(Dim::opt.solver,1,1);
    else
      this->initSystem(Dim::opt.solver,1,0);
    this->setQuadratureRule(Dim::opt.nGauss[0],true);

    if (!this->assembleSystem())
      return false;
    else if (vizRHS && Dim::opt.eig == 0)
      this->extractLoadVec(myLoad);

    if (Dim::opt.eig > 0) // Eigenvalue analysis (free vibration)
      return this->systemModes(modes);

    if (!this->solveSystem(mySolVec,1))
      return false;

    if (!this->setMode(SIM::RECOVERY))
      return false;

    if (!Dim::opt.project.empty())
    {
      // Project the secondary solution onto the splines basis
      size_t j = 0;
      for (auto& pit : Dim::opt.project)
        if (!this->project(myProj[j++],mySolVec,pit.first))
          return false;

      IFEM::cout << std::endl;
    }

    // Evaluate solution norms
    Vectors gNorm;
    this->setQuadratureRule(Dim::opt.nGauss[1]);
    if (!this->solutionNorms(mySolVec,myProj,myNorm,gNorm))
      return false;

    // Print global norm summary to console
    this->printNorms(gNorm);
    return true;
  }

  //! \brief Saves solution-dependent quantities to file for postprocessing.
  bool saveStep(TimeStep&, int& nBlock)
  {
    if (!asciiFile.empty())
    {
      // Write solution (control point values) to ASCII files
      std::ofstream osd(asciiFile+".sol");
      osd.precision(18);
      IFEM::cout <<"\nWriting primary solution (temperature field) to file "
                 << asciiFile+".sol" << std::endl;
      utl::LogStream log(osd);
      this->dumpPrimSol(mySolVec,log,false);
      std::ofstream oss(asciiFile+".sec");
      oss.precision(18);
      IFEM::cout <<"\nWriting all solutions (heat flux etc) to file "
                 << asciiFile+".sec" << std::endl;
      utl::LogStream log2(oss);
      this->dumpSolution(mySolVec,log2);
    }

    if (Dim::opt.format < 0)
      return true;

    // Write load vector to VTF-file
    if (!this->writeGlvV(myLoad,"Load vector",1,nBlock))
      return false;

    // Write solution fields to VTF-file
    if (!this->writeGlvS(mySolVec,1,nBlock))
      return false;

    // Write projected solution fields to VTF-file
    size_t i = 0;
    int iBlk = 100, iGrad = -1;
    std::string grdName;
    std::vector<std::string> prefix(Dim::opt.project.size());
    for (auto& pit : Dim::opt.project)
      if (i >= myProj.size())
        break;
      else if (!this->writeGlvP(myProj[i],1,nBlock,iBlk,pit.second.c_str()))
        return false;
      else
      {
        iBlk += 10;
        prefix[i++] = pit.second.c_str();
        if (iGrad < 0 || pit.first == SIMoptions::GLOBAL)
        {
          iGrad = i-1;
          grdName = pit.second + " q";
        }
      }

    // Write the projected solution gradient vector (heat flux) to VTF-file
    if (iGrad >= 0)
      if (!this->writeGlvV(myProj[iGrad],grdName.c_str(),1,nBlock,110,Dim::nsd))
        return false;

    // Write eigenmodes
    bool isFreq = Dim::opt.eig==3 || Dim::opt.eig==4 || Dim::opt.eig==6;
    for (i = 0; i < modes.size(); i++)
      if (!this->writeGlvM(modes[i],isFreq,nBlock))
        return false;

    // Write element norms
    if (!this->writeGlvN(myNorm,1,nBlock,prefix))
      return false;

    return this->writeGlvStep(1,0.0,1);
  }

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
  virtual bool solveSystem(Vector& solution, int printSol, double* rCond,
                           const char* compName, bool newLHS, size_t idxRHS)
  {
    if (!this->Dim::solveSystem(solution,printSol,rCond,compName,newLHS,idxRHS))
      return false;
    else if (idxRHS > 0 || !this->haveReactions() || prob.extEner != 'R')
      return true;

    // Assemble the reaction forces. Strictly, we only need to assemble those
    // elements that have nodes on the Dirichlet boundaries, but...
    prob.setReactionIntegral(new ReactionsOnly(myReact,Dim::mySam));
    AlgEqSystem* tmpEqSys = Dim::myEqSys;
    Dim::myEqSys = nullptr;
    bool ok = this->setMode(SIM::RHS_ONLY) && this->assembleSystem({solution});
    Dim::myEqSys = tmpEqSys;
    prob.setReactionIntegral(nullptr);

    return ok;
  }

  //! \brief Returns current reaction force vector.
  virtual const Vector* getReactionForces() const
  {
    return myReact.empty() ? nullptr : &myReact;
  }

  //! \brief Prints a norm group to the log stream.
  //! \param[in] gNorm The global norm values
  //! \param[in] fNorm Global reference norm values
  //! \param[in] name Name of norm group
  virtual void printNormGroup(const Vector& gNorm, const Vector& fNorm,
                              const std::string& name) const
  {
    double Rel = 100.0/(this->haveAnaSol() ? fNorm(3) : gNorm(1));
    const char* uRef = this->haveAnaSol() ? "|u|)  " : "|u^r|)";
    IFEM::cout <<"\n>>> Error estimates based on "<< name <<" <<<";
    if (name == "Pure residuals")
      IFEM::cout <<"\nResidual norm |u|_res = |f+nabla^2 u|: "<< gNorm(2);
    else
      IFEM::cout <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< gNorm(1)
                 <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< gNorm(2)
                 <<"\n- relative error (% of "<< uRef <<" : "<< gNorm(2)*Rel
                 <<"\nResidual error (r(u^r) + J(u^r))^0.5 : "<< gNorm(3)
                 <<"\n- relative error (% of "<< uRef <<" : "<< gNorm(3)*Rel;

    if (this->haveAnaSol())
    {
      if (gNorm.size() > 3 && name != "Pure residuals")
        IFEM::cout <<"\nExact error a(e,e)^0.5, e=u-u^r      : "<< gNorm(4)
                   <<"\n- relative error (% of |u|)   : "<< gNorm(4)*Rel;
      if (fNorm.size() > 3 && gNorm.size() > 3)
        IFEM::cout <<"\nEffectivity index, theta^*           : "
                   << gNorm(2)/fNorm(4)
                   <<"\nEffectivity index, theta^EX          : "
                   << (gNorm(2)+gNorm(4))/fNorm(4)
                   <<"\nEffectivity index, theta^RES         : "
                   << (gNorm(2)+gNorm(3))/fNorm(4);
    }
    IFEM::cout << std::endl;
  }

  //! \brief Returns the name of this simulator.
  //! \details This method is typically reimplemented in sub-classes that are
  //! parts of a partitioned solution method and are used to identify the basis
  //! for the result fields associated with each simulator in the HDF5 output.
  virtual std::string getName() const { return "Poisson"; }

protected:
  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to resolve inhomogeneous boundary
  //! condition fields in case they are derived from the analytical solution.
  virtual void preprocessA()
  {
    myProj.resize(Dim::opt.project.size());
    if (!Dim::mySol) return;

    // Define analytical boundary condition fields
    PropertyVec::iterator p;
    for (p = Dim::myProps.begin(); p != Dim::myProps.end(); ++p)
      if (p->pcode == Property::DIRICHLET_ANASOL)
      {
        if (!Dim::mySol->getScalarSol())
          p->pcode = Property::UNDEFINED;
        else if (aCode[0] == abs(p->pindx))
          p->pcode = Property::DIRICHLET_INHOM;
        else if (aCode[0] == 0)
        {
          aCode[0] = abs(p->pindx);
          Dim::myScalars[aCode[0]] = Dim::mySol->getScalarSol();
          p->pcode = Property::DIRICHLET_INHOM;
        }
        else
          p->pcode = Property::UNDEFINED;
      }
      else if (p->pcode == Property::NEUMANN_ANASOL)
      {
        if (!Dim::mySol->getScalarSecSol())
          p->pcode = Property::UNDEFINED;
        else if (aCode[1] == p->pindx)
          p->pcode = Property::NEUMANN;
        else if (aCode[1] == 0)
        {
          aCode[1] = p->pindx;
          Dim::myVectors[aCode[1]] = Dim::mySol->getScalarSecSol();
          p->pcode = Property::NEUMANN;
        }
        else
          p->pcode = Property::UNDEFINED;
      }
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  virtual bool preprocessB()
  {
    // Check if the model has constraints.
    // If not, we can calculate external energy also without reaction forces.
    if (this->getNoConstraints() == 0 && !prob.extEner)
      prob.extEner = 'y';
    return true;
  }

  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is)
  {
    if (this->parseDimSpecific(keyWord,is))
      return true;

    else if (!strncasecmp(keyWord,"ISOTROPIC",9))
    {
      char* cline = 0;
      int nmat = atoi(keyWord+10);
      std::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;
      for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
      {
        int    code  = atoi(strtok(cline," "));
        double kappa = atof(strtok(nullptr," "));
        if (code == 0)
          prob.setMaterial(kappa);
        else
          this->setPropertyType(code,Property::MATERIAL,mVec.size());
        mVec.push_back(kappa);
        std::cout <<"\tMaterial code "<< code <<": "<< kappa << std::endl;
      }
    }

    else
      return this->Dim::parse(keyWord,is);

    return true;
  }

  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (!strcasecmp(elem->Value(),"postprocessing"))
      prob.parse(elem->FirstChildElement("projection"));

    if (strcasecmp(elem->Value(),"poisson"))
      return this->Dim::parse(elem);

    bool result = true;
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (this->parseDimSpecific(child))
        continue;

      else if (!strcasecmp(child->Value(),"isotropic")) {
        int code = this->parseMaterialSet(child,mVec.size());
        double kappa = 1000.0;
        utl::getAttribute(child,"kappa",kappa);
        if (code == 0)
          prob.setMaterial(kappa);
        mVec.push_back(kappa);
        std::cout <<"\tMaterial code "<< code <<": "<< kappa << std::endl;
      }

      else if (!strcasecmp(child->Value(),"reactions"))
        prob.extEner = 'R';
      else if (!prob.parse(child))
        result &= this->Dim::parse(child);

    return result;
  }

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
  virtual bool initMaterial(size_t propInd)
  {
    if (propInd >= mVec.size()) return false;

    prob.setMaterial(mVec[propInd]);
    return true;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator sit = Dim::myScalars.find(propInd);
    typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);

    if (sit != Dim::myScalars.end())
      prob.setTraction(sit->second);
    else if (vit != Dim::myVectors.end())
      prob.setTraction(vit->second);
    else
      return false;

    return true;
  }

private:
  Poisson   prob;     //!< Data and methods for the Poisson problem
  RealArray mVec;     //!< Material data
  int       aCode[2]; //!< Analytical BC code (used by destructor)

  Vector    myLoad;   //!< External load vector (for VTF export)
  Vector    mySolVec; //!< Primary solution vector
  Vector    myReact;  //!< Nodal reaction forces
  Vectors   myProj;   //!< Projected solution vectors
  Matrix    myNorm;   //!< Element norms

  std::vector<Mode> modes; //!< Eigen modes

  const Vector* solution; //!< Pointer to primary solution vector

  bool        vizRHS;    //!< If \e true, store load vector to VTF
  std::string asciiFile; //!< ASCII output file prefix
};


typedef SIMPoisson<SIM1D> SIMPoisson1D; //!< 1D specific driver
typedef SIMPoisson<SIM2D> SIMPoisson2D; //!< 2D specific driver
typedef SIMPoisson<SIM3D> SIMPoisson3D; //!< 3D specific driver

//! \brief Template specialization - 1D specific input parsing.
template<> bool SIMPoisson1D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 1D specific input parsing.
template<> bool SIMPoisson1D::parseDimSpecific(const TiXmlElement* elem);

//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMPoisson2D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMPoisson2D::parseDimSpecific(const TiXmlElement* elem);

//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMPoisson3D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMPoisson3D::parseDimSpecific(const TiXmlElement* elem);

#endif
