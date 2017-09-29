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
#include "AnaSol.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "Utilities.h"
#include "tinyxml.h"
#include "TimeStep.h"
#include "Functions.h"
#include "ExprFunctions.h"
#include "DataExporter.h"
#include "IFEM.h"
#include <fstream>


/*!
  \brief Driver class for NURBS-based FEM analysis of Poisson problems.
*/

template<class Dim> class SIMPoisson : public SIMMultiPatchModelGen<Dim>
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  explicit SIMPoisson(bool checkRHS = false) :
    SIMMultiPatchModelGen<Dim>(1,checkRHS), prob(Dim::dimension)
  {
    Dim::myProblem = &prob;
    aCode[0] = aCode[1] = 0;
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

  //! \brief Returns number of right-hand-side vectors
  size_t getNoRHS() const { return 1 + prob.getNoGalerkin(); }

  //! \brief Register data fields for output
  void registerFields(DataExporter& exporter)
  {
    int results = DataExporter::PRIMARY | DataExporter::NORMS;

    if (!this->opt.pSolOnly && this->opt.project.empty())
      results |= DataExporter::SECONDARY;

    if (solution) {
      exporter.registerField("u", "solution", DataExporter::SIM, results);
      exporter.setFieldValue("u", this, solution);
    }

    if (!this->opt.project.empty() && projections) {
      size_t i = 0;
      for (auto& it : this->opt.project) {
        exporter.registerField(it.second.c_str(), "projected", DataExporter::SIM,
                                DataExporter::SECONDARY, it.second.c_str());
        exporter.setFieldValue(it.second.c_str(), this, &(*projections)[i++]);
      }
    }
  }

  //! \brief Set solution vector
  void setSol(const Vector* sol) { solution = const_cast<Vector*>(sol); }

  //! \brief Set project vectors
  void setProjections(const Vectors* sol) { projections = sol; }

  //! \brief No internal VTF handling.
  bool saveModel(char* infile, int& geoBlk, int&)
  {
    if (Dim::opt.format >= 0 && !this->writeGlvG(geoBlk, infile))
      return false;

    this->gBlk = geoBlk;

    return true;
  }

  //! \brief Assemble and solve linear system
  bool solveStep(TimeStep&)
  {
    this->setMode(this->opt.eig == 0 ? SIM::STATIC : SIM::VIBRATION);
    if (!this->assembleSystem())
      return false;

    if (this->opt.eig == 0) {
      if (!this->solveSystem(mySolVec,1))
        return false;

      // Project the FE stresses onto the splines basis
      this->setMode(SIM::RECOVERY);
      size_t i = 0;
      for (auto pit  = this->opt.project.begin();
                pit != this->opt.project.end(); ++i, ++pit) {
        Matrix ssol;
        if (!this->project(ssol,*solution,pit->first))
          return false;
        else
          myProj.push_back(ssol);
      }

      // Evaluate solution norms
      this->setQuadratureRule(this->opt.nGauss[1]);
      if (!this->solutionNorms(Vectors(1,*solution),*projections,eNorm,gNorm))
        return false;

      this->printSummary();
    } else {
      this->setMode(SIM::VIBRATION);
      if (!this->systemModes(modes))
        return false;
    }

    return true;
  }

  //! \brief No time stepping.
  bool advanceStep(TimeStep&) { return true; }

  //! \brief No internal VTF handling.
  bool saveStep(TimeStep& tp, int& nBlock)
  {
    if (this->opt.format < 0)
      return true;

    // Write boundary tractions, if any
    if (!this->writeGlvT(1,gBlk,nBlock))
      return false;

    // Write Dirichlet boundary conditions
    if (!this->writeGlvBC(nBlock))
      return false;

    // Write load vector to VTF-file
    if (vizRHS) {
      Vector load;
      this->extractLoadVec(load);
      if (!this->writeGlvV(load,"Load vector",1,nBlock))
        return false;
    }

    // Write solution fields to VTF-file
    if (!this->writeGlvS(*solution,1,nBlock))
      return false;

    // Write solution gradients to VTF-file
    if (!solution->empty()) {
      Matrix tmp;
      if (!this->project(tmp, *solution))
        return false;
      this->writeGlvV(tmp, "grad(u)", tp.step, nBlock, 110, this->nsd);
    }

    // Write projected solution fields to VTF-file
    size_t i = 0;
    int iBlk = 100;
    for (auto pit = this->opt.project.begin(); pit != this->opt.project.end(); pit++, i++, iBlk += 10)
      if (!this->writeGlvP((*projections)[i],1,nBlock,iBlk,pit->second.c_str()))
        return false;

    // Write eigenmodes
    bool isFreq = this->opt.eig==3 || this->opt.eig==4 || this->opt.eig==6;
    for (i = 0; i < modes.size(); i++)
      if (!this->writeGlvM(modes[i],isFreq,nBlock))
        return false;

    // Write element norms
    auto prefix = this->getNormPrefixes();
    if (!this->writeGlvN(eNorm,1,nBlock,prefix.data()))
      return false;

    this->writeGlvStep(1,0.0,1);

    if (!asciiFile.empty())
      this->writeASCII();

    return true;
  }

  //! \brief Print a summary of the solution to terminal
  void printSummary()
  {
    this->printNorms(gNorm);
    size_t j = 1;

    for (auto& pit : this->opt.project)
    {
      IFEM::cout <<"\n>>> Error estimates based on "<< pit.second <<" <<<";
      IFEM::cout <<"\nEnergy norm |u^r| = a(u^r,u^r)^0.5   : "<< gNorm[j](1);
      IFEM::cout <<"\nError norm a(e,e)^0.5, e=u^r-u^h     : "<< gNorm[j](2);
      IFEM::cout <<"\n- relative error (% of |u^r|) : "
                 << gNorm[j](2)/gNorm[j](1)*100.0;
      if (this->haveAnaSol() && j <= gNorm.size())
      {
        IFEM::cout <<"\nExact error a(e,e)^0.5, e=u-u^r      : "<< gNorm[j](3)
                   <<"\n- relative error (% of |u|)   : "
                   << gNorm[j](3)/gNorm[0](3)*100.0;
        IFEM::cout <<"\nEffectivity index             : "
                   << gNorm[j](2)/gNorm[0](4);
      }
      IFEM::cout << std::endl;
      ++j;
    }
  }

  //! \brief Set whether or not to dump results to ascii files
  void setDumpASCII(bool dump) { asciiFile.resize(dump?1:0); }

  //! \brief True to store laod vector to VTF file
  void setVizRHS(bool viz) { vizRHS = viz; }

  //! \brief Write results to ascii files
  void writeASCII()
  {
    // Write (refined) model to g2-file
    std::ofstream osg(asciiFile+".g2");
    osg.precision(18);
    IFEM::cout <<"\nWriting updated g2-file "<< asciiFile+".g2" << std::endl;
    this->dumpGeometry(osg);
    if (solution && !solution->empty())
    {
      // Write solution (control point values) to ASCII files
      std::ofstream osd(asciiFile+".sol");
      osd.precision(18);
      IFEM::cout <<"\nWriting primary solution (temperature field) to file "
                 << asciiFile+".sol" << std::endl;
      utl::LogStream log(osd);
      this->dumpPrimSol(*solution,log,false);
      std::ofstream oss(asciiFile+".sec");
      oss.precision(18);
      IFEM::cout <<"\nWriting all solutions (heat flux etc) to file "
                 << asciiFile+".sec" << std::endl;
      utl::LogStream log2(oss);
      this->dumpSolution(*solution,log2);
    }
  }

  //! \brief Read input file and store name of the parsed file.
  virtual bool read(const char* filename)
  {
    if (!asciiFile.empty()) {
      asciiFile = filename;
      asciiFile = asciiFile.substr(asciiFile.rfind('.')+1);
    }

    return this->SIMadmin::read(filename);
  }


protected:
  //! \brief Internal helper function to extract projection prefixes.
  std::vector<const char*> getNormPrefixes()
  {
    std::vector<const char*> result;
    for (auto& it : this->opt.project)
      result.push_back(it.second.c_str());

    return result;
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to resolve inhomogeneous boundary
  //! condition fields in case they are derived from the analytical solution.
  virtual void preprocessA()
  {
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
    if (strcasecmp(elem->Value(),"poisson"))
      return this->Dim::parse(elem);

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
      else if (!strcasecmp(child->Value(),"galerkin")) {
        if (child->FirstChild() && child->FirstChild()->Value())
          prob.addGalerkin(new VecFuncExpr(child->FirstChild()->Value()));
      }

      else
        this->Dim::parse(child);

    return true;
  }

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
  //! \brief Parses a dimension-specific data section from an input stream.
  //! \details This function allows for specialization of the template while
  //! still reusing as much code as possible. Only for dimension-specific code.
  bool parseDimSpecific(char* keyWord, std::istream& is);
  //! \brief Parses a dimension-specific data section from the an XML element.
  //! \details This function allows for specialization of the template while
  //! still reusing as much code as possible. Only for dimension-specific code.
  bool parseDimSpecific(const TiXmlElement* child);

  Poisson   prob;     //!< Data and methods for the Poisson problem
  RealArray mVec;     //!< Material data
  int       aCode[2]; //!< Analytical BC code (used by destructor)
  int gBlk = -1;      //!< Geometry block
  Matrix eNorm;       //!< Element norms
  Vectors gNorm;      //!< Global norms
  bool vizRHS = false;//!< True to store load vector to VTF file

  Vectors myProj; //!< Internal projection vectors
  const Vectors* projections = &myProj; //!< Pointer to projection vectors
  Vector mySolVec; //!< Internal solution vector
  const Vector* solution = &mySolVec; //!< Pointer to solution vector
  std::vector<Mode> modes; //!< Eigen modes
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
