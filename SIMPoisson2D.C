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
  myProblem = 0;
  // To avoid that that SIMbase tries to delete already deleted functions
  if (myAFcode > 0) myVectors.erase(myAFcode);
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
      std::cout <<"\tMaterial code "<< code <<": "<< kappa << std::endl;
      if (code == 0)
	prob.setMaterial(kappa);
      else if (this->setPropertyType(code,Property::MATERIAL,mVec.size()))
	mVec.push_back(kappa);
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
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"SQUARE",6))
    {
      double L = atof(strtok(NULL," "));
      std::cout <<"\nHeat source function: Square L="<< L << std::endl;
      prob.setSource(new Square2DHeat(L));
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      cline = strtok(NULL," ");
      std::cout <<"\nHeat source function: " << cline << std::endl;
      prob.setSource(new EvalFunction(cline));
    }
    else
      std::cerr <<"  ** SIMPoisson2D::parse: Unknown source function "
		<< cline << std::endl;
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    int code = -1;
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"SQUARE",6))
    {
      double L = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Square L="<< L << std::endl;
      mySol = new AnaSol(NULL,new Square2D(L));
    }
    else if (!strncasecmp(cline,"LSHAPE",6))
    {
      std::cout <<"\nAnalytical solution: Lshape"<< std::endl;
      mySol = new AnaSol(NULL,new LshapePoisson());
    }
    else if (!strncasecmp(cline,"SINUSSQUARE",11))
    {
      std::cout <<"\nAnalytical solution: SquareSinus"
		<<"\nHeat source function: SquareSinusSource"<< std::endl;
      mySol = new AnaSol(NULL,new SquareSinus());
      prob.setSource(new SquareSinusSource());
    }
    else if (!strncasecmp(cline,"INTERIORLAYER",13))
    {
      std::cout <<"\nAnalytical solution: InteriorLayer"
		<<"\nHeat source function: InteriorLayerSource"<< std::endl;
      mySol = new AnaSol(new PoissonInteriorLayerSol(),
			 new PoissonInteriorLayer());
      prob.setSource(new PoissonInteriorLayerSource());

      // Define the Dirichlet boundary condition from the analytical solution
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      if (code < 1)
      {
        std::cerr <<" *** SIMPoisson2D::parse: Specify code > 0 for the"
		  <<" inhomogenous DIRICHLET boundary on InteriorLayer\n";
        return false;
      }
      this->setPropertyType(code,Property::DIRICHLET_INHOM);
      myScalars[code] = new PoissonInteriorLayerSol();
      code = 0; // Avoid definition of Neumann property
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::cout <<"\nAnalytical solution: Expression"<< std::endl;
      int lines = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
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
      myAFcode = code;
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
      int code = 0;
      double kappa = 1000.0;
      utl::getAttribute(child,"code",code);
      utl::getAttribute(child,"kappa",kappa);
      if (code == 0)
        prob.setMaterial(kappa);
      else if (setPropertyType(code,Property::MATERIAL,mVec.size()))
        mVec.push_back(kappa);
    }

    else if (!strcasecmp(child->Value(),"source")) {
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "square") {
        double L = 0.0;
        utl::getAttribute(child,"L",L);
        std::cout <<"\nHeat source function: Square L="<< L << std::endl;
        prob.setSource(new Square2DHeat(L));
      }
      else if (type == "expression" && child->FirstChild()) {
        std::cout <<"\nHeat source function: "
                  << child->FirstChild()->Value() << std::endl;
        prob.setSource(new EvalFunction(child->FirstChild()->Value()));
      }
      else
        std::cerr <<"  ** SIMPoisson2D::parse: Invalid source function "
                  << type << std::endl;
    }

    else if (!strcasecmp(child->Value(),"anasol")) {
      int code = 0;
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "square") {
        double L = 0.0;
        utl::getAttribute(child,"L",L);
        std::cout <<"\nAnalytical solution: Square L="<< L << std::endl;
        mySol = new AnaSol(NULL,new Square2D(L));
      }
      else if (type == "lshape") {
        std::cout <<"\nAnalytical solution: Lshape"<< std::endl;
        mySol = new AnaSol(NULL,new LshapePoisson());
      }
      else if (type == "sinussquare") {
        std::cout <<"\nAnalytical solution: SquareSinus"
                  <<"\nHeat source function: SquareSinusSource"<< std::endl;
        mySol = new AnaSol(NULL,new SquareSinus());
        prob.setSource(new SquareSinusSource());
      }
      else if (type == "interiorlayer") {
        std::cout <<"\nAnalytical solution: InteriorLayer"
                  <<"\nHeat source function: InteriorLayerSource"<< std::endl;
        mySol = new AnaSol(new PoissonInteriorLayerSol(),
                           new PoissonInteriorLayer());
        prob.setSource(new PoissonInteriorLayerSource());

        // Define the Dirichlet boundary condition from the analytical solution
        utl::getAttribute(child,"code",code);
        if (code < 1)
        {
          std::cerr <<" *** SIMPoisson2D::parse: Specify code > 0 for the"
                    <<" inhomogenous DIRICHLET boundary on InteriorLayer\n";
          return false;
        }
        this->setPropertyType(code,Property::DIRICHLET_INHOM);
        myScalars[code] = new PoissonInteriorLayerSol();
      }
      else if (type == "expression") {
        std::cout <<"\nAnalytical solution: Expression"<< std::endl;
        mySol = new AnaSol(child);
      }
      else
        std::cerr <<"  ** SIMPoisson2D::parse: Invalid analytical solution "
                  << type <<" (ignored)"<< std::endl;

      // Define the analytical boundary traction field
      if (code == 0 && utl::getAttribute(child,"code",code))
        if (code > 0 && mySol && mySol->getScalarSecSol())
        {
          setPropertyType(code,Property::NEUMANN);
          myVectors[code] = mySol->getScalarSecSol();
          myAFcode = code;
        }
    }

  return true;
}


bool SIMPoisson2D::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) return false;

  prob.setMaterial(mVec[propInd]);
  return true;
}


bool SIMPoisson2D::initNeumann (size_t propInd)
{
  VecFuncMap::const_iterator tit = myVectors.find(propInd);
  if (tit == myVectors.end()) return false;

  prob.setTraction(tit->second);
  return true;
}
