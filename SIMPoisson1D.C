// $Id$
//==============================================================================
//!
//! \file SIMPoisson1D.C
//!
//! \date Apr 16 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for 1D NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#include "SIMPoisson1D.h"
#include "PoissonSolutions.h"
#include "Functions.h"
#include "Utilities.h"
#include "AnaSol.h"
#include "tinyxml.h"


SIMPoisson1D::~SIMPoisson1D ()
{
  myProblem = 0;
  // To avoid that that SIMbase tries to delete already deleted functions
  if (myAFcode > 0) myVectors.erase(myAFcode);
}


bool SIMPoisson1D::parse (char* keyWord, std::istream& is)
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

  else if (!strncasecmp(keyWord,"SOURCE",6))
  {
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"LINE",4))
    {
      double L = atof(strtok(NULL," "));
      std::cout <<"\nHeat source function: Line L="<< L << std::endl;
      prob.setSource(new PoissonLineSource(L));
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      cline = strtok(NULL," ");
      std::cout <<"\nHeat source function: " << cline << std::endl;
      prob.setSource(new EvalFunction(cline));
    }
    else
      std::cerr <<"  ** SIMPoisson1D::parse: Unknown source function "
		<< cline << std::endl;
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    int code = -1;
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"LINE",4))
    {
      double L = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Line L="<< L << std::endl;
      mySol = new AnaSol(NULL,new PoissonLine(L));
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
      std::cerr <<"  ** SIMPoisson1D::parse: Unknown analytical solution "
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
    return this->SIM1D::parse(keyWord,is);

  return true;
}


bool SIMPoisson1D::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"poisson"))
    return this->SIM1D::parse(elem);

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
      if (type == "line") {
        double L = 0.0;
        utl::getAttribute(child,"L",L);
        std::cout <<"\nHeat source function: Line L="<< L << std::endl;
        prob.setSource(new PoissonLineSource(L));
      }
      else if (type == "expression" && child->FirstChild()) {
        std::cout <<"\nHeat source function: "
                  << child->FirstChild()->Value() << std::endl;
        prob.setSource(new EvalFunction(child->FirstChild()->Value()));
      }
      else
        std::cerr <<"  ** SIMPoisson1D::parse: Invalid source function "
                  << type << std::endl;
    }

    else if (!strcasecmp(child->Value(),"anasol")) {
      int code = 0;
      std::string type;
      utl::getAttribute(child,"type",type,true);
      if (type == "line") {
        double L = 0.0;
        utl::getAttribute(child,"L",L);
        std::cout <<"\nAnalytical solution: Line L="<< L << std::endl;
        mySol = new AnaSol(NULL,new PoissonLine(L));
      }
      else if (type == "expression") {
        std::cout <<"\nAnalytical solution: Expression"<< std::endl;
        mySol = new AnaSol(child);
      }
      else
        std::cerr <<"  ** SIMPoisson1D::parse: Invalid analytical solution "
                  << type <<" (ignored)"<< std::endl;

      // Define the analytical boundary traction field
      if (utl::getAttribute(child,"code",code))
	if (code > 0 && mySol && mySol->getScalarSecSol())
        {
          setPropertyType(code,Property::NEUMANN);
          myVectors[code] = mySol->getScalarSecSol();
          myAFcode = code;
        }
    }

  return true;
}


bool SIMPoisson1D::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) return false;

  prob.setMaterial(mVec[propInd]);
  return true;
}


bool SIMPoisson1D::initNeumann (size_t propInd)
{
  VecFuncMap::const_iterator tit = myVectors.find(propInd);
  if (tit == myVectors.end()) return false;

  prob.setTraction(tit->second);
  return true;
}
