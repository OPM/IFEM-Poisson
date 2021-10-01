// $Id$
//==============================================================================
//!
//! \file SIMPoisson.C
//!
//! \date May 25 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#include "SIMPoisson.h"
#include "PoissonSolutions.h"

#include "AnaSol.h"
#include "DataExporter.h"
#include "ExprFunctions.h"
#include "IFEM.h"
#include "ReactionsOnly.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "Utilities.h"

#ifdef HAS_LRSPLINE
#include "ASMu2D.h"
#include "Profiler.h"
#endif

#include <fstream>
#include "tinyxml.h"


template<class Dim>
SIMPoisson<Dim>::SIMPoisson (bool checkRHS, bool ds)
  : SIMMultiPatchModelGen<Dim>(1,checkRHS),
    prob(Dim::dimension),
    robinBC(Dim::dimension, prob),
    solution(&mySolVec)
{
  Dim::myProblem = &prob;
  aCode[0] = aCode[1] = 0;
  vizRHS = false;
  dualS = ds;
}


template<class Dim>
SIMPoisson<Dim>::~SIMPoisson ()
{
  Dim::myProblem = nullptr; // Because it is not dynamically allocated
  Dim::myInts.clear();

  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
  if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
}


template<class Dim>
void SIMPoisson<Dim>::clearProperties ()
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
  prob.addExtrFunction(nullptr);
  this->Dim::clearProperties();
}


template<class Dim>
size_t SIMPoisson<Dim>::getNoRHS () const
{
  return dualS ? 2 : 1 + prob.getNoGalerkin();
}


template<class Dim>
bool SIMPoisson<Dim>::haveDualSol () const
{
  return dualS && Dim::dualField;
}


template<class Dim>
void SIMPoisson<Dim>::registerFields (DataExporter& exporter)
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
                           Dim::opt.saveNorms ? &myNorm : nullptr);
  }
}


template<class Dim>
void SIMPoisson<Dim>::setASCIIfile (const char* filename)
{
  const char* end = strrchr(filename,'.');
  if (end)
    asciiFile.assign(filename,end);
  else
    asciiFile.assign(filename);
}


template<class Dim>
bool SIMPoisson<Dim>::saveModel (char* fileName, int& geoBlk, int& nBlock)
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


template<class Dim>
bool SIMPoisson<Dim>::solveStep (TimeStep&)
{
  if (!this->setMode(Dim::opt.eig == 0 ? SIM::STATIC : SIM::VIBRATION))
    return false;

  if (!this->initSystem(Dim::opt.solver, 1, Dim::opt.eig == 0 ? 1 : 0))
    return false;

  this->setQuadratureRule(Dim::opt.nGauss[0],true);
  if (!this->assembleSystem())
    return false;

  if (Dim::opt.eig > 0) // Eigenvalue analysis (free vibration)
    return this->systemModes(modes);
  else if (vizRHS)
    this->extractLoadVec(myLoad);

  if (!this->solveSystem(mySolVec,1))
    return false;

  if (!Dim::opt.project.empty())
  {
    this->setMode(SIM::RECOVERY);

    // Project the secondary solution onto the splines basis
    size_t j = 0;
    for (const SIMoptions::ProjectionMap::value_type& pit : Dim::opt.project)
      if (!this->project(myProj[j++],mySolVec,pit.first))
        return false;

    IFEM::cout << std::endl;
  }

  // Evaluate solution norms
  Vectors gNorm;
  this->setMode(SIM::NORMS);
  this->setQuadratureRule(Dim::opt.nGauss[1]);
  if (!this->solutionNorms(mySolVec,myProj,myNorm,gNorm))
    return false;

  // Print global norm summary to console
  this->printNorms(gNorm);

  if (this->hasResultPoints())
  {
    // Print point-wise result quantities
    this->setMode(SIM::RECOVERY);
    this->dumpResults(mySolVec,0.0,IFEM::cout,true);
    if (!myProj.empty())
      this->dumpVector(myProj.front(),nullptr,IFEM::cout);
  }

  return true;
}



template<class Dim>
bool SIMPoisson<Dim>::saveStep (TimeStep&, int& nBlock)
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

  size_t pos = 0;
  for (const Kappa& f : mVec)
    if (f.func && !this->writeGlvF(*f.func, ("kappa" + std::to_string(++pos)).c_str(), 1, nBlock))
      return false;

  // Write projected solution fields to VTF-file
  size_t i = 0;
  int iBlk = 100, iGrad = -1;
  std::string grdName;
  std::vector<std::string> prefix(Dim::opt.project.size());
  for (const SIMoptions::ProjectionMap::value_type& pit : Dim::opt.project)
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


template<class Dim>
bool SIMPoisson<Dim>::solveSystem (Vector& solution, int printSol,
                                   double* rCond, const char* compName,
                                   bool newLHS, size_t idxRHS)
{
  if (!this->Dim::solveSystem(solution,printSol,rCond,compName,newLHS,idxRHS))
    return false;
  else if (idxRHS > 0 || !this->haveReactions() || prob.extEner != 'R')
    return true;

  // Assemble the reaction forces. Strictly, we only need to assemble those
  // elements that have nodes on the Dirichlet boundaries, but...
  prob.setReactionIntegral(new ReactionsOnly(myReact,Dim::mySam,Dim::adm));
  AlgEqSystem* tmpEqSys = Dim::myEqSys;
  Dim::myEqSys = nullptr;
  bool ok = this->setMode(SIM::RHS_ONLY) && this->assembleSystem({solution});
  Dim::myEqSys = tmpEqSys;
  prob.setReactionIntegral(nullptr);

  return ok;
}


template<class Dim>
void SIMPoisson<Dim>::printNormGroup (const Vector& gNorm, const Vector& fNorm,
                                      const std::string& name) const
{
  IFEM::cout <<"\n>>> Error estimates based on "<< name <<" <<<";

  if (name == "Pure residuals") {
    IFEM::cout <<"\nResidual norm |u|_res = |f+nabla^2 u|: "<< gNorm(2);
    if (!this->haveAnaSol()) {
      IFEM::cout << std::endl;
      return;
    }
  }

  double Rel = 100.0/(this->haveAnaSol() ? fNorm(3) : gNorm(1));
  const char* uRef = this->haveAnaSol() ? "|u|)  " : "|u^r|)";
  if (name != "Pure residuals")
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


template<class Dim>
void SIMPoisson<Dim>::preprocessA ()
{
  if (Dim::dualField)
    prob.setDualRHS(Dim::dualField);

  myProj.resize(Dim::opt.project.size());
  if (!Dim::mySol) return;

  Dim::myInts.insert(std::make_pair(0,Dim::myProblem));

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
    } else if (p->pcode == Property::ROBIN)
      if (Dim::myInts.find(p->pindx) == Dim::myInts.end())
        Dim::myInts.insert(std::make_pair(p->pindx,&robinBC));
}


template<class Dim>
bool SIMPoisson<Dim>::preprocessB ()
{
  // Check if the model has constraints.
  // If not, we can calculate external energy also without reaction forces.
  if (this->getNoConstraints() == 0 && !prob.extEner)
    prob.extEner = 'y';
  return true;
}


template<class Dim>
bool SIMPoisson<Dim>::parse (char* keyWord, std::istream& is)
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
      mVec.push_back(Kappa{kappa, nullptr});
      std::cout <<"\tMaterial code "<< code <<": "<< kappa << std::endl;
    }
  }

  else
    return this->Dim::parse(keyWord,is);

  return true;
}


template<class Dim>
bool SIMPoisson<Dim>::parse (const TiXmlElement* elem)
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
      Kappa kappa{1000.0, nullptr};
      utl::getAttribute(child,"kappa",kappa.constant);
      if (code == 0)
        prob.setMaterial(kappa.constant);
      mVec.push_back(kappa);
      std::cout <<"\tMaterial code "<< code <<": "<< kappa.constant << std::endl;
    }
    else if (!strcasecmp(child->Value(),"propertymaterial")) {
      int code = this->parseMaterialSet(child,mVec.size());
      tprops.parse(child);
      if (tprops.hasProperty("kappa")) {
        Kappa kappa;
        kappa.func.reset(new PropertyFunc("kappa", tprops));
        mVec.push_back(kappa);
        if (code == 0)
          prob.setMaterial(mVec.back().func.get());
        std::cout <<"\tMaterial code "<< code <<": property"<< std::endl;
      }
    }

    else if (!strcasecmp(child->Value(),"reactions"))
      prob.extEner = 'R';
    else if (!strcasecmp(child->Value(),"dualfield"))
      prob.addExtrFunction(this->parseDualTag(child));
    else if (!prob.parse(child))
      result &= this->Dim::parse(child);

  return result;
}


template<class Dim>
bool SIMPoisson<Dim>::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) return false;

  if (mVec[propInd].func)
    prob.setMaterial(mVec[propInd].func.get());
  else
    prob.setMaterial(mVec[propInd].constant);

  return true;
}


template<class Dim>
bool SIMPoisson<Dim>::initNeumann (size_t propInd)
{
  typename Dim::SclFuncMap::const_iterator sit = Dim::myScalars.find(propInd);
  typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);

  if (sit != Dim::myScalars.end()) {
    prob.setTraction(sit->second);
    robinBC.setFlux(sit->second);
  } else if (vit != Dim::myVectors.end()) {
    prob.setTraction(vit->second);
    robinBC.setAlpha(vit->second);
  } else
    return false;

  return true;
}


//! \brief Template specialization - 1D specific input parsing.
template<>
bool SIMPoisson<SIM1D>::parseDimSpecific (char* keyWord, std::istream& is)
{
  char* cline;
  if (!strncasecmp(keyWord,"SOURCE",6))
  {
    int code = -1; // Reserve negative code(s) for the source term function
    while (myScalars.find(code) != myScalars.end()) --code;

    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"LINE",4))
    {
      double L = atof(strtok(nullptr," "));
      std::cout <<"\nHeat source function: Line L="<< L << std::endl;
      myScalars[code] = new PoissonLineSource(L);
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      cline = strtok(nullptr," ");
      std::cout <<"\nHeat source function: "<< cline << std::endl;
      myScalars[code] = new EvalFunction(cline);
    }
    else
    {
      std::cerr <<"  ** SIMPoisson1D::parse: Invalid source function "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }
    prob.setSource(myScalars[code]);
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    int code = -1;
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"LINE",4))
    {
      double L = atof(strtok(nullptr," "));
      std::cout <<"\nAnalytical solution: Line L="<< L << std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new PoissonLine(L));
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      if (!mySol)
        mySol = new AnaSol(is,lines);
    }
    else
    {
      std::cerr <<"  ** SIMPoisson1D::parse: Invalid analytical solution "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    if (code == -1)
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
    if (code > 0 && mySol->getScalarSecSol())
    {
      this->setPropertyType(code,Property::NEUMANN);
      myVectors[code] = mySol->getScalarSecSol();
      aCode[1] = code;
    }
  }
  else
    return false;

  return true;
}


//! \brief Template specialization - 1D specific input parsing.
template<>
bool SIMPoisson<SIM1D>::parseDimSpecific (const TiXmlElement* child)
{
  if (!strcasecmp(child->Value(),"source")) {
    int code = -1; // Reserve negative code(s) for the source term function
    while (myScalars.find(code) != myScalars.end()) --code;
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (type == "line") {
      double L = 0.0;
      utl::getAttribute(child,"L",L);
      std::cout <<"\tHeat source function: Line L="<< L << std::endl;
      myScalars[code] = new PoissonLineSource(L);
    }
    else if (type == "expression" && child->FirstChild()) {
      std::cout <<"\tHeat source function: "
                << child->FirstChild()->Value() << std::endl;
      myScalars[code] = new EvalFunction(child->FirstChild()->Value());
    }
    else
    {
      std::cerr <<"  ** SIMPoisson1D::parse: Invalid source function "
                << type <<" (ignored)"<< std::endl;
      return false;
    }
    prob.setSource(myScalars[code]);
  }

  else if (!strcasecmp(child->Value(),"anasol")) {
    int code = 0;
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (type == "line") {
      double L = 0.0;
      utl::getAttribute(child,"L",L);
      std::cout <<"\tAnalytical solution: Line L="<< L << std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new PoissonLine(L));
    }
    else if (type == "expression" || type == "fields") {
      type[0] = toupper(type[0]);
      std::cout <<"\tAnalytical solution: "<< type << std::endl;
      if (!mySol)
        mySol = new AnaSol(child);
    }
    else
      std::cerr <<"  ** SIMPoisson1D::parse: Invalid analytical solution "
                << type <<" (ignored)"<< std::endl;

    // Define the analytical boundary traction field
    if (code == 0 && utl::getAttribute(child,"code",code))
      if (code > 0 && mySol && mySol->getScalarSecSol())
      {
        this->setPropertyType(code,Property::NEUMANN);
        myVectors[code] = mySol->getScalarSecSol();
        aCode[1] = code;
      }
  }
  else
    return false;

  return true;
}


//! \brief Template specialization - 2D specific input parsing.
template<>
bool SIMPoisson<SIM2D>::parseDimSpecific (char* keyWord, std::istream& is)
{
  char* cline;

  // Test code, move to SIM2D maybe?)
  if (!strncasecmp(keyWord,"LRREFINE",8))
  {
#ifdef HAS_LRSPLINE
    ASMu2D* patch = dynamic_cast<ASMu2D*>(myModel.front());
    if (patch)
    {
      cline = strtok(keyWord+8," ");
      if (!strncasecmp(cline,"UNIFORM",7))
      {
        PROFILE("LR refinement");
        int nRef = atoi(strtok(nullptr," "));
        std::cout <<"\nLR refinement UNIFORM : "<< nRef << std::endl;
        patch->uniformRefine(nRef);
      }
      else if (!strncasecmp(cline,"CORNER",6))
      {
        PROFILE("LR refinement");
        int nRef = atoi(strtok(nullptr," "));
        std::cout <<"\nLR refinement CORNER : "<< nRef << std::endl;
        patch->cornerRefine(nRef);
      }
      else if (!strncasecmp(cline,"DIAGONAL",8))
      {
        PROFILE("LR refinement");
        int nRef = atoi(strtok(nullptr," "));
        std::cout <<"\nLR refinement DIAGONAL : "<< nRef << std::endl;
        patch->diagonalRefine(nRef);
      }
    }
#endif
  }

  else if (!strncasecmp(keyWord,"SOURCE",6))
  {
    int code = -1; // Reserve negative code(s) for the source term function
    while (myScalars.find(code) != myScalars.end()) --code;

    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"SQUARE",6))
    {
      double L = atof(strtok(nullptr," "));
      std::cout <<"\nHeat source function: Square L="<< L << std::endl;
      myScalars[code] = new Square2DHeat(L);
    }
    else if (!strncasecmp(cline,"SINUSSQUARE",11))
    {
      std::cout <<"\nHeat source function: SquareSinus"<< std::endl;
      myScalars[code] = new SquareSinusSource();
    }
    else if (!strncasecmp(cline,"INTERIORLAYER",13))
    {
      std::cout <<"\nHeat source function: InteriorLayer"<< std::endl;
      myScalars[code] = new PoissonInteriorLayerSource();
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      cline = strtok(nullptr," ");
      std::cout <<"\nHeat source function: "<< cline << std::endl;
      myScalars[code] = new EvalFunction(cline);
    }
    else
    {
      std::cerr <<"  ** SIMPoisson2D::parse: Invalid source function "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }
    prob.setSource(myScalars[code]);
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    int code = -1;
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"SQUARE",6))
    {
      double L = atof(strtok(nullptr," "));
      std::cout <<"\nAnalytical solution: Square L="<< L << std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new Square2D(L));
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      std::cout <<"\nAnalytical solution: Lshape"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new LshapePoisson());
    }
    else if (!strncasecmp(cline,"SINUSSQUARE",11))
    {
      std::cout <<"\nAnalytical solution: SquareSinus"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new SquareSinus());
    }
    else if (!strncasecmp(cline,"INTERIORLAYER",13))
    {
      std::cout <<"\nAnalytical solution: InteriorLayer"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(new PoissonInteriorLayerSol(),
                           new PoissonInteriorLayer());

      // Define the Dirichlet boundary condition from the analytical solution
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      if (code < 1)
      {
        std::cerr <<" *** SIMPoisson2D::parse: Specify code > 0 for the"
                  <<" inhomogenous DIRICHLET boundary on InteriorLayer\n";
        return false;
      }
      this->setPropertyType(code,Property::DIRICHLET_INHOM);
      myScalars[code] = mySol->getScalarSol();
      aCode[0] = code;
      code = 0; // Avoid definition of Neumann property
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      if (!mySol)
        mySol = new AnaSol(is,lines);
    }
    else
    {
      std::cerr <<"  ** SIMPoisson2D::parse: Invalid analytical solution "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    if (code == -1)
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
    if (code > 0 && mySol->getScalarSecSol())
    {
      this->setPropertyType(code,Property::NEUMANN);
      myVectors[code] = mySol->getScalarSecSol();
      aCode[1] = code;
    }
  }
  else
    return false;

  return true;
}


//! \brief Template specialization - 2D specific input parsing.
template<>
bool SIMPoisson<SIM2D>::parseDimSpecific (const TiXmlElement* child)
{
  if (!strcasecmp(child->Value(),"source")) {
    int code = -1; // Reserve negative code(s) for the source term function
    while (myScalars.find(code) != myScalars.end()) --code;
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (type == "square") {
      double L = 0.0;
      utl::getAttribute(child,"L",L);
      std::cout <<"\tHeat source function: Square L="<< L << std::endl;
      myScalars[code] = new Square2DHeat(L);
    }
    else if (type == "sinussquare") {
      std::cout <<"\tHeat source function: SquareSinus"<< std::endl;
      myScalars[code] = new SquareSinusSource();
    }
    else if (type == "interiorlayer") {
      double s = 60;
      utl::getAttribute(child,"slope",s);
      std::cout <<"\tHeat source function: InteriorLayer, slope="<< s
                << std::endl;
      myScalars[code] = new PoissonInteriorLayerSource(s);
    }
    else if (type == "expression" && child->FirstChild()) {
      std::cout <<"\tHeat source function: "
                << child->FirstChild()->Value() << std::endl;
      myScalars[code] = new EvalFunction(child->FirstChild()->Value());
    }
    else
    {
      std::cerr <<"  ** SIMPoisson2D::parse: Invalid source function "
                << type <<" (ignored)"<< std::endl;
      return true;
    }

    prob.setSource(myScalars[code]);
  }

  else if (!strcasecmp(child->Value(),"anasol")) {
    int code = 0;
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (type == "square") {
      double L = 0.0;
      utl::getAttribute(child,"L",L);
      std::cout <<"\tAnalytical solution: Square L="<< L << std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new Square2D(L));
    }
    else if (type == "lshape") {
      std::cout <<"\tAnalytical solution: Lshape"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new LshapePoisson());
    }
    else if (type == "sinussquare") {
      std::cout <<"\tAnalytical solution: SquareSinus"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new SquareSinus());
    }
    else if (type == "interiorlayer") {
      double s = 60;
      utl::getAttribute(child,"slope",s);
      std::cout <<"\tAnalytical solution: InteriorLayer, slope="<< s
                << std::endl;
      if (!mySol)
        mySol = new AnaSol(new PoissonInteriorLayerSol(s),
                           new PoissonInteriorLayer(s));

      // Define the Dirichlet boundary condition from the analytical solution
      utl::getAttribute(child,"code",code);
      if (code > 0)
      {
        this->setPropertyType(code,Property::DIRICHLET_INHOM);
        myScalars[code] = mySol->getScalarSol();
        aCode[0] = code;
      }
    }
    else if (type == "expression" || type == "fields") {
      type[0] = toupper(type[0]);
      std::cout <<"\tAnalytical solution: "<< type << std::endl;
      if (!mySol)
        mySol = new AnaSol(child);
    }
    else
      std::cerr <<"  ** SIMPoisson2D::parse: Invalid analytical solution "
                << type <<" (ignored)"<< std::endl;

    // Define the analytical boundary traction field
    if (code == 0 && utl::getAttribute(child,"code",code))
      if (code > 0 && mySol && mySol->getScalarSecSol())
      {
        this->setPropertyType(code,Property::NEUMANN);
        myVectors[code] = mySol->getScalarSecSol();
        aCode[1] = code;
      }
  }
  else
    return false;

  return true;
}


//! \brief Template specialization - 3D specific input parsing.
template<>
bool SIMPoisson<SIM3D>::parseDimSpecific (char* keyWord, std::istream& is)
{
  char* cline;
  if (!strncasecmp(keyWord,"SOURCE",6))
  {
    int code = -1; // Reserve negative code(s) for the source term function
    while (myScalars.find(code) != myScalars.end()) --code;

    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"CUBE",4))
    {
      std::cout <<"\nHeat source function: Cube"<< std::endl;
      myScalars[code] = new PoissonCubeSource();
    }
    else if (!strncasecmp(cline,"WATERFALL",9))
    {
      std::cout <<"\nHeat source function: Waterfall"<< std::endl;
      myScalars[code] = new PoissonWaterfallSource();
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      cline = strtok(nullptr," ");
      std::cout <<"\nHeat source function: "<< cline << std::endl;
      myScalars[code] = new EvalFunction(cline);
    }
    else
    {
      std::cerr <<"  ** SIMPoisson3D::parse: Invalid source function "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }
    prob.setSource(myScalars[code]);
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    int code = -1;
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"CUBE",4))
    {
      std::cout <<"\nAnalytical solution: Cube"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new PoissonCube());
    }
    else if (!strncasecmp(cline,"WATERFALL",9))
    {
      std::cout <<"\nAnalytical solution: Waterfall"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(new PoissonWaterfallSol(),
                           new PoissonWaterfall());

      // Define the Dirichlet boundary condition from the analytical solution
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      if (code < 1)
      {
        std::cerr <<" *** SIMPoisson2D::parse: Specify code > 0 for the"
                  <<" inhomogenous DIRICHLET boundary on Waterfall\n";
        return false;
      }
      this->setPropertyType(code,Property::DIRICHLET_INHOM);
      myScalars[code] = mySol->getScalarSol();
      aCode[0] = code;
      code = 0; // Avoid definition of Neumann property
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
      if (!mySol)
        mySol = new AnaSol(is,lines);
    }
    else
    {
      std::cerr <<"  ** SIMPoisson3D::parse: Invalid analytical solution "
                << cline <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    if (code == -1)
      code = (cline = strtok(nullptr," ")) ? atoi(cline) : 0;
    if (code > 0 && mySol->getScalarSecSol())
    {
      this->setPropertyType(code,Property::NEUMANN);
      myVectors[code] = mySol->getScalarSecSol();
      aCode[1] = code;
    }
  }
  else
    return false;

  return true;
}


//! \brief Template specialization - 3D specific input parsing.
template<>
bool SIMPoisson<SIM3D>::parseDimSpecific (const TiXmlElement* child)
{
  if (!strcasecmp(child->Value(),"source")) {
    int code = -1; // Reserve negative code(s) for the source term function
    while (myScalars.find(code) != myScalars.end()) --code;
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (type == "cube") {
      std::cout <<"\tHeat source function: Cube"<< std::endl;
      myScalars[code] = new PoissonCubeSource();
    }
    else if (type == "waterfall") {
      double eps = 0.002;
      utl::getAttribute(child,"epsilon",eps);
      std::cout <<"\tHeat source function: Waterfall, epsilon="<< eps
                << std::endl;
      myScalars[code] = new PoissonWaterfallSource(eps);
    }
    else if (type == "expression" && child->FirstChild()) {
      std::cout <<"\tHeat source function: "
                << child->FirstChild()->Value() << std::endl;
      myScalars[code] = new EvalFunction(child->FirstChild()->Value());
    }
    else
    {
      std::cerr <<"  ** SIMPoisson3D::parse: Invalid source function "
                << type <<" (ignored)"<< std::endl;
      return true;
    }
    prob.setSource(myScalars[code]);
  }

  else if (!strcasecmp(child->Value(),"anasol")) {
    int code = 0;
    std::string type;
    utl::getAttribute(child,"type",type,true);
    if (type == "cube") {
      std::cout <<"\tAnalytical solution: Cube"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(nullptr,new PoissonCube());
    }
    else if (type == "waterfall") {
      double eps = 0.002;
      utl::getAttribute(child,"epsilon",eps);
      std::cout <<"\tAnalytical solution: Waterfall, epsilon="<< eps
                << std::endl;
      if (!mySol)
        mySol = new AnaSol(new PoissonWaterfallSol(eps),
                           new PoissonWaterfall(eps));

      // Define the Dirichlet boundary condition from the analytical solution
      utl::getAttribute(child,"code",code);
      if (code > 0)
      {
        this->setPropertyType(code,Property::DIRICHLET_INHOM);
        myScalars[code] = mySol->getScalarSol();
        aCode[0] = code;
      }
    }
    else if (type == "expression" || type == "fields") {
      type[0] = toupper(type[0]);
      std::cout <<"\tAnalytical solution: "<< type << std::endl;
      if (!mySol)
        mySol = new AnaSol(child);
    }
    else
      std::cerr <<"  ** SIMPoisson3D::parse: Invalid analytical solution "
                << type <<" (ignored)"<< std::endl;

    // Define the analytical boundary traction field
    if (code == 0 && utl::getAttribute(child,"code",code))
      if (code > 0 && mySol && mySol->getScalarSecSol())
      {
        this->setPropertyType(code,Property::NEUMANN);
        myVectors[code] = mySol->getScalarSecSol();
        aCode[1] = code;
      }
  }
  else
    return false;

  return true;
}


template class SIMPoisson<SIM1D>;
template class SIMPoisson<SIM2D>;
template class SIMPoisson<SIM3D>;
