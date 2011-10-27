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
#include "Utilities.h"
#include "AnaSol.h"
#ifdef HAS_LRSPLINE
#include "LR/ASMu2D.h"
#include "Profiler.h"
#endif
#include <string.h>


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
    else
      std::cerr <<"  ** SIMPoisson2D::parse: Unknown source function "
		<< cline << std::endl;
  }

  else if (!strncasecmp(keyWord,"ANASOL",6))
  {
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
      size_t code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
      if (code < 1)
      {
        std::cerr <<"  ** SIMPoisson2D::parse: Specify code > 0 for the"
		  <<" inhomogenous DIRICHLET boundary on InteriorLayer\n";
        return false;
      }
      for (size_t j = 0; j < myProps.size(); j++)
        if (myProps[j].pindx == code && myProps[j].pcode == Property::UNDEFINED)
          myProps[j].pcode = Property::DIRICHLET_INHOM;
      myScalars[code] = new PoissonInteriorLayerSol();
    }
    else
    {
      std::cerr <<"  ** SIMPoisson2D::parse: Unknown analytical solution "
		<< cline <<" (ignored)"<< std::endl;
      return true;
    }

    // Define the analytical boundary traction field
    int code = (cline = strtok(NULL," ")) ? atoi(cline) : 0;
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
