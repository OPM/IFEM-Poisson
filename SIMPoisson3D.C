// $Id$
//==============================================================================
//!
//! \file SIMPoisson3D.C
//!
//! \date May 25 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for 3D NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#include "SIMPoisson3D.h"
#include "PoissonSolutions.h"
#include "Utilities.h"
#include "AnaSol.h"
#include <string.h>
#include "Functions.h"


SIMPoisson3D::~SIMPoisson3D ()
{
  myProblem = 0;
  // To avoid that that SIMbase tries to delete already deleted functions
  if (myAFcode > 0) myVectors.erase(myAFcode);
}


bool SIMPoisson3D::parse (char* keyWord, std::istream& is)
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
    if (!strncasecmp(cline,"CUBE",4))
    {
      std::cout <<"\nHeat source function: Cube"<< std::endl;
      prob.setSource(new PoissonCubeSource());
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      cline = strtok(NULL," ");
      std::cout <<"\nHeat source function: " << cline << std::endl;
      prob.setSource(new EvalFunction(cline));
    }
    else
      std::cerr <<"  ** SIMPoisson3D::parse: Unknown source function "
		<< cline << std::endl;
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    cline = strtok(keyWord+6," ");
    int code=-1;
    if (!strncasecmp(cline,"CUBE",4))
    {
      mySol = new AnaSol(NULL,new PoissonCube());
      std::cout <<"\nAnalytical solution: Cube"<< std::endl;
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      std::string primary, secondary, variables;
      int lines = atoi(strtok(NULL, " "));
      char* c = strtok(NULL, " ");
      if (c)
        code = atoi(c);
      else
        code = 0;
      RealFunc* s=NULL;
      VecFunc* v=NULL;
      for (int i = 0; i < lines; i++) {
        std::string function = utl::readLine(is);
        size_t pos;
        if ((pos = function.find("Variables=")) != std::string::npos) {
          variables += function.substr(pos+10);
          if (variables[variables.size()-1] != ';')
            variables += ";";
        }
        if ((pos = function.find("Primary=")) != std::string::npos) {
          primary = function.substr(pos+8);
          s = new EvalFunction((variables+primary).c_str());
        }
        if ((pos = function.find("Secondary=")) != std::string::npos) {
          secondary = function.substr(pos+10);
          v = new VecFuncExpr(secondary,variables);
        }
      }
      std::cout <<"\nAnalytical solution:" << std::endl;
      if (!variables.empty())
        std::cout << "\t Variables=" << variables << std::endl;
      if (s)
        std::cout << "\t Primary=" << primary << std::endl;
      if (v)
        std::cout << "\t Secondary=" << secondary << std::endl;
      mySol = new AnaSol(s, v);
    }
    else
    {
      std::cerr <<"  ** SIMPoisson3D::parse: Unknown analytical solution "
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
    return this->SIM3D::parse(keyWord,is);

  return true;
}


bool SIMPoisson3D::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) return false;

  prob.setMaterial(mVec[propInd]);
  return true;
}


bool SIMPoisson3D::initNeumann (size_t propInd)
{
  VecFuncMap::const_iterator tit = myVectors.find(propInd);
  if (tit == myVectors.end()) return false;

  prob.setTraction(tit->second);
  return true;
}
