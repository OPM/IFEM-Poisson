// $Id$
//==============================================================================
//!
//! \file SIMPoisson3D.C
//!
//! \date May 25 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief 3D specific code for NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#include "SIMPoisson.h"
#include "PoissonSolutions.h"
#include "Functions.h"
#include "Utilities.h"
#include "AnaSol.h"
#include "tinyxml.h"


  template<>
bool SIMPoisson3D::parseDimSpecific(char* keyWord, std::istream& is)
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


  template<>
bool SIMPoisson3D::parseDimSpecific(const TiXmlElement* child)
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
