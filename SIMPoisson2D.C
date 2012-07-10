// $Id$
//==============================================================================
//!
//! \file SIMPoisson2D.C
//!
//! \date Apr 16 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for 2D NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#include "SIMPoisson2D.h"
#include "PoissonSolutions.h"
#include "Functions.h"
#include "Utilities.h"
#include "AnaSol.h"
#include "tinyxml.h"

#ifdef HAS_LRSPLINE
#include "LR/ASMu2D.h"
#include "Profiler.h"
#endif


SIMPoisson2D::~SIMPoisson2D ()
{
  myProblem = NULL; // Because it is not dynamically allocated

  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode[0] > 0) myScalars.erase(aCode[0]);
  if (aCode[1] > 0) myVectors.erase(aCode[1]);
}


void SIMPoisson2D::clearProperties ()
{
  // To prevent SIMbase::clearProperties deleting the analytical solution
  if (aCode[0] > 0) myScalars.erase(aCode[0]);
  if (aCode[1] > 0) myVectors.erase(aCode[1]);
  aCode[0] = aCode[1] = 0;

  mVec.clear();
  prob.setSource(NULL);
  prob.setTraction((RealFunc*)NULL);
  prob.setTraction((VecFunc*)NULL);
  this->SIMbase::clearProperties();
}


bool SIMPoisson2D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;

  if (!strncasecmp(keyWord,"ISOTROPIC",9))
  {
    int nmat = atoi(keyWord+10);
    std::cout <<"\nNumber of isotropic materials: "<< nmat << std::endl;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int    code  = atoi(strtok(cline," "));
      double kappa = atof(strtok(NULL," "));
      if (code == 0)
        prob.setMaterial(kappa);
      else
        this->setPropertyType(code,Property::MATERIAL,mVec.size());
      mVec.push_back(kappa);
      std::cout <<"\tMaterial code "<< code <<": "<< kappa << std::endl;
    }
  }

  // Test code, move to SIM2D maybe?)
  else if (!strncasecmp(keyWord,"LRREFINE",8))
  {
#ifdef HAS_LRSPLINE
    ASMu2D* patch = dynamic_cast<ASMu2D*>(myModel.front());
    if (patch)
    {
      cline = strtok(keyWord+8," ");
      if (!strncasecmp(cline,"UNIFORM",7))
      {
	PROFILE("LR refinement");
	int nRef = atoi(strtok(NULL," "));
	std::cout <<"\nLR refinement UNIFORM : "<< nRef << std::endl;
	patch->uniformRefine(nRef);
      }
      else if (!strncasecmp(cline,"CORNER",6))
      {
	PROFILE("LR refinement");
	int nRef = atoi(strtok(NULL," "));
	std::cout <<"\nLR refinement CORNER : "<< nRef << std::endl;
	patch->cornerRefine(nRef);
      }
      else if (!strncasecmp(cline,"DIAGONAL",8))
      {
	PROFILE("LR refinement");
	int nRef = atoi(strtok(NULL," "));
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
      double L = atof(strtok(NULL," "));
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
      cline = strtok(NULL," ");
      std::cout <<"\nHeat source function: " << cline << std::endl;
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
      double L = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Square L="<< L << std::endl;
      if (!mySol)
        mySol = new AnaSol(NULL,new Square2D(L));
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      std::cout <<"\nAnalytical solution: Lshape"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(NULL,new LshapePoisson());
    }
    else if (!strncasecmp(cline,"SINUSSQUARE",11))
    {
      std::cout <<"\nAnalytical solution: SquareSinus"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(NULL,new SquareSinus());
    }
    else if (!strncasecmp(cline,"INTERIORLAYER",13))
    {
      std::cout <<"\nAnalytical solution: InteriorLayer"<< std::endl;
      if (!mySol)
        mySol = new AnaSol(new PoissonInteriorLayerSol(),
                           new PoissonInteriorLayer());

      // Define the Dirichlet boundary condition from the analytical solution
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
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
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
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
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
    if (code > 0 && mySol->getScalarSecSol())
    {
      this->setPropertyType(code,Property::NEUMANN);
      myVectors[code] = mySol->getScalarSecSol();
      aCode[1] = code;
    }
  }

  else
    return this->SIM2D::parse(keyWord,is);

  return true;
}


bool SIMPoisson2D::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"poisson"))
    return this->SIM2D::parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"isotropic")) {
      int code = this->parseMaterialSet(child,mVec.size());
      double kappa = 1000.0;
      utl::getAttribute(child,"kappa",kappa);
      if (code == 0)
        prob.setMaterial(kappa);
      mVec.push_back(kappa);
      std::cout <<"\tMaterial code "<< code <<": "<< kappa << std::endl;
    }

    else if (!strcasecmp(child->Value(),"source")) {
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
        std::cout <<"\tHeat source function: InteriorLayer, slope = "<< s
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
        continue;
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
          mySol = new AnaSol(NULL,new Square2D(L));
      }
      else if (type == "lshape") {
        std::cout <<"\tAnalytical solution: Lshape"<< std::endl;
        if (!mySol)
          mySol = new AnaSol(NULL,new LshapePoisson());
      }
      else if (type == "sinussquare") {
        std::cout <<"\tAnalytical solution: SquareSinus"<< std::endl;
        if (!mySol)
          mySol = new AnaSol(NULL,new SquareSinus());
      }
      else if (type == "interiorlayer") {
        double s = 60;
        utl::getAttribute(child,"slope",s);
        std::cout <<"\tAnalytical solution: InteriorLayer, slope = "<< s
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
      else if (type == "expression") {
        std::cout <<"\tAnalytical solution: Expression"<< std::endl;
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

  return true;
}


bool SIMPoisson2D::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  if (mySol) // Define analytical boundary condition fields
    for (PropertyVec::iterator p = myProps.begin(); p != myProps.end(); p++)
      if (p->pcode == Property::DIRICHLET_ANASOL)
      {
        if (!mySol->getScalarSol())
          p->pcode = Property::UNDEFINED;
	else if (aCode[0] == abs(p->pindx))
          p->pcode = Property::DIRICHLET_INHOM;
	else if (aCode[0] == 0)
        {
          aCode[0] = abs(p->pindx);
          myScalars[aCode[0]] = mySol->getScalarSol();
          p->pcode = Property::DIRICHLET_INHOM;
        }
        else
          p->pcode = Property::UNDEFINED;
      }
      else if (p->pcode == Property::NEUMANN_ANASOL)
      {
        if (!mySol->getScalarSecSol())
          p->pcode = Property::UNDEFINED;
	else if (aCode[1] == p->pindx)
          p->pcode = Property::NEUMANN;
	else if (aCode[1] == 0)
        {
          aCode[1] = p->pindx;
          myVectors[aCode[1]] = mySol->getScalarSecSol();
          p->pcode = Property::NEUMANN;
        }
        else
          p->pcode = Property::UNDEFINED;
      }

  return this->SIM2D::preprocess(ignored,fixDup);
}


bool SIMPoisson2D::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) return false;

  prob.setMaterial(mVec[propInd]);
  return true;
}


bool SIMPoisson2D::initNeumann (size_t propInd)
{
  SclFuncMap::const_iterator sit = myScalars.find(propInd);
  VecFuncMap::const_iterator vit = myVectors.find(propInd);

  if (sit != myScalars.end())
    prob.setTraction(sit->second);
  else if (vit != myVectors.end())
    prob.setTraction(vit->second);
  else
    return false;

  return true;
}


std::ostream& SIMPoisson2D::printNorms (const Vectors& norms, std::ostream& os)
{
  if (norms.empty()) return os;

  NormBase* norm = this->getNormIntegrand();
  const Vector& gnorm = norms.front();

  os <<"Energy norm "<< norm->getName(1,1) <<": "<< gnorm(1)
     <<"\nExternal energy "<< norm->getName(1,2) <<": "<< gnorm(2);

  if (mySol)
    os <<"\nExact norm "<< norm->getName(1,3) <<": "<< gnorm(3)
       <<"\nExact error "<< norm->getName(1,4) <<": "<< gnorm(4)
       <<"\nExact relative error (%) : "<< 100.0*gnorm(4)/gnorm(3);

  delete norm;

  return os << std::endl;
}
