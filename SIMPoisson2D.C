// $Id$
//==============================================================================
//!
//! \file SIMPoisson2D.C
//!
//! \date Apr 16 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief 2D specific code for NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#include "SIMPoisson.h"
#include "PoissonSolutions.h"
#include "Functions.h"
#include "Utilities.h"
#include "AnaSol.h"
#include "tinyxml.h"

#ifdef HAS_LRSPLINE
#include "LR/ASMu2D.h"
#include "Profiler.h"
#endif


  template<>
bool SIMPoisson2D::parseDimSpecific(char* keyWord, std::istream& is)
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

  } else
    return false;

  return true;
}


  template<>
bool SIMPoisson2D::parseDimSpecific(const TiXmlElement* child)
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
  } else
    return false;

  return true;
}
