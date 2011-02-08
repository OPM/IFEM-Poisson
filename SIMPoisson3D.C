// $Id: SIMPoisson3D.C,v 1.3 2010-09-05 12:44:18 kmo Exp $
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
#include "AnalyticSolutions.h"
#include "Utilities.h"
#include <string.h>


SIMPoisson3D::SIMPoisson3D (bool checkRHS) : SIM3D(checkRHS,1)
{
  myProblem = &poPrb;
  asol = 0;
}


SIMPoisson3D::~SIMPoisson3D ()
{
  myProblem = 0;
  if (asol) delete asol;
}


bool SIMPoisson3D::parse (char* keyWord, std::istream& is)
{
  char* cline = 0;

  if (!strncasecmp(keyWord,"ISOTROPHIC",10))
  {
    int nmat = atoi(keyWord+10);
    std::cout <<"\nNumber of isotrophic materials: "<< nmat << std::endl;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      int    code  = atoi(strtok(cline," "));
      double kappa = atof(strtok(NULL," "));
      std::cout <<"\tMaterial code "<< code <<":"<< kappa << std::endl;
      if (code == 0)
	poPrb.setMaterial(kappa);
      else for (unsigned int j = 0; j < myProps.size(); j++)
	if (myProps[j].pindx == code && myProps[j].pcode == Property::UNDEFINED)
	{
	  myProps[j].pindx = mVec.size();
	  myProps[j].pcode = Property::MATERIAL;
	}
      if (code > 0)
	mVec.push_back(kappa);
    }
  }

  else if (!strncasecmp(keyWord,"SOURCE",6))
  {
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"CUBE",4))
    {
      std::cout <<"\nHeat source function: Cube"<< std::endl;
      poPrb.setSource(new PoissonCubeSource());
    }
    else
      std::cerr <<"  ** SIMPoisson3D::parse: Unknown source function "
		<< cline << std::endl;
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"CUBE",4))
    {
      asol = new PoissonCube();
      std::cout <<"\nAnalytical solution: Cube"<< std::endl;
    }
    else
      std::cerr <<"  ** SIMPoisson3D::parse: Unknown analytical solution "
		<< cline << std::endl;

    // Define the analytical boundary traction field
    int code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
    if (code > 0 && asol)
    {
      for (unsigned int j = 0; j < myProps.size(); j++)
	if (myProps[j].pindx == code && myProps[j].pcode == Property::UNDEFINED)
	  myProps[j].pcode = Property::NEUMANN;

      myVectors[code] = asol;
    }
  }


  // The remaining keywords are retained for backward compatibility with the
  // prototype version. They enable direct specification of properties onto
  // the topological entities (blocks and faces) of the model.

  else if (!strncasecmp(keyWord,"MATERIAL",8))
  {
    int nmat = atoi(keyWord+8);
    std::cout <<"\nNumber of materials: "<< nmat << std::endl;
    for (int i = 0; i < nmat && (cline = utl::readLine(is)); i++)
    {
      double kappa = atof(strtok(cline," "));
      while (cline = strtok(NULL," "))
	if (!strncasecmp(cline,"ALL",3))
        {
	  std::cout <<"\tMaterial for all patches: "<< kappa << std::endl;
	  poPrb.setMaterial(kappa);
	}
	else
        {
	  int patch = atoi(cline);
	  if (patch < 1 || patch > myModel.size())
	  {
	    std::cerr <<" *** SIMPoisson3D::parse: Invalid patch index "
		      << patch << std::endl;
	    return false;
	  }
	  std::cout <<"\tMaterial for P"<< patch <<": "<< kappa  << std::endl;
	  myProps.push_back(Property(Property::MATERIAL,mVec.size(),patch,3));
	  mVec.push_back(kappa);
	}
    }
  }

  else if (!strncasecmp(keyWord,"EXNEUMANN",9))
  {
    Property neum;
    neum.pcode = Property::NEUMANN;
    neum.ldim = 2;

    int npres = atoi(keyWord+9);
    std::cout <<"\nNumber of Neumann integrals: "<< npres << std::endl;
    for (int i = 0; i < npres && (cline = utl::readLine(is)); i++)
    {
      neum.pindx = 1+i;
      neum.patch = atoi(strtok(cline," "));
      if (neum.patch < 1 || neum.patch > myModel.size())
      {
	std::cerr <<" *** SIMPoisson3D::parse: Invalid patch index "
		  << neum.patch << std::endl;
	return false;
      }

      neum.lindx = atoi(strtok(NULL," "));
      if (neum.lindx < 1 || neum.lindx > 6)
      {
	std::cerr <<" *** SIMPoisson3D::parse: Invalid face index "
		  << (int)neum.lindx << std::endl;
	return false;
      }

      if (asol)
      {
	std::cout <<"\tNeumann integral on P"<< neum.patch
		  <<" F"<< (int)neum.lindx << std::endl;
	myVectors[1+i] = asol;
      }

      myProps.push_back(neum);
    }
  }

  else
    return this->SIM3D::parse(keyWord,is);

  return true;
}


bool SIMPoisson3D::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) return false;

  poPrb.setMaterial(mVec[propInd]);
  return true;
}


bool SIMPoisson3D::initNeumann (size_t propInd)
{
  VecFuncMap::const_iterator tit = myVectors.find(propInd);
  if (tit == myVectors.end()) return false;

  poPrb.setTraction(tit->second);
  return true;
}
