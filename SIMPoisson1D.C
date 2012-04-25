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
  myProblem = NULL; // Because it is not dynamically allocated

  // To prevent the SIMbase destructor try to delete already deleted functions
  if (aCode[0] > 0) myScalars.erase(aCode[0]);
  if (aCode[1] > 0) myVectors.erase(aCode[1]);
}


void SIMPoisson1D::clearProperties ()
{
  // To prevent SIMbase::clearProperties deleting the analytical solution
  if (aCode[0] > 0) myScalars.erase(aCode[0]);
  if (aCode[1] > 0) myVectors.erase(aCode[1]);

  mVec.clear();
  prob.setSource(NULL);
  prob.setTraction((RealFunc*)NULL);
  prob.setTraction((VecFunc*)NULL);
  this->SIMbase::clearProperties();
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
      if (code == 0)
        prob.setMaterial(kappa);
      else
        this->setPropertyType(code,Property::MATERIAL,mVec.size());
      mVec.push_back(kappa);
      std::cout <<"\tMaterial code "<< code <<": "<< kappa << std::endl;
    }
  }

  else if (!strncasecmp(keyWord,"SOURCE",6))
  {
    int code = -1; // Reserve negative code(s) for the source term function
    while (myScalars.find(code) != myScalars.end()) --code;

    cline = strtok(keyWord+6," ");
    if (!strncasecmp(cline,"LINE",4))
    {
      double L = atof(strtok(NULL," "));
      std::cout <<"\nHeat source function: Line L="<< L << std::endl;
      myScalars[code] = new PoissonLineSource(L);
    }
    else if (!strncasecmp(cline,"EXPRESSION",10))
    {
      cline = strtok(NULL," ");
      std::cout <<"\nHeat source function: " << cline << std::endl;
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
      double L = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Line L="<< L << std::endl;
      if (!mySol)
        mySol = new AnaSol(NULL,new PoissonLine(L));
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
      std::cerr <<"  ** SIMPoisson1D::parse: Invalid analytical solution "
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
      else
        this->setPropertyType(code,Property::MATERIAL,mVec.size());
      mVec.push_back(kappa);
    }

    else if (!strcasecmp(child->Value(),"source")) {
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
        continue;
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
          mySol = new AnaSol(NULL,new PoissonLine(L));
      }
      else if (type == "expression") {
        std::cout <<"\tAnalytical solution: Expression"<< std::endl;
        if (!mySol)
          mySol = new AnaSol(child);
      }
      else
        std::cerr <<"  ** SIMPoisson1D::parse: Invalid analytical solution "
                  << type <<" (ignored)"<< std::endl;

      // Define the analytical boundary traction field
      if (utl::getAttribute(child,"code",code))
        if (code > 0 && mySol && mySol->getScalarSecSol())
        {
          this->setPropertyType(code,Property::NEUMANN);
          myVectors[code] = mySol->getScalarSecSol();
          aCode[1] = code;
        }
    }

  return true;
}


bool SIMPoisson1D::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  bool ok = true;

  if (mySol) // Define analytical boundary condition fields
    for (PropertyVec::iterator p = myProps.begin(); p != myProps.end(); p++)
      if (p->pcode == Property::DIRICHLET_ANASOL)
      {
        if (mySol->getScalarSol() && (aCode[0]==0 || aCode[0]==abs(p->pindx)))
        {
          p->pcode = Property::DIRICHLET_INHOM;
          myScalars[abs(p->pindx)] = mySol->getScalarSol();
          aCode[0] = abs(p->pindx);
        }
        else
        {
          p->pcode = Property::UNDEFINED;
          std::cerr <<" *** SIMPoisson1D::preprocess: Analytic Dirichlet condit"
                    <<"ons\n     can only be assigned to one topology set.\n"
                    <<"     Ignoring specification for property code = "
                    << p->pindx << std::endl;
          ok = false;
        }
      }
      else if (p->pcode == Property::NEUMANN_ANASOL)
      {
        if (mySol->getScalarSecSol() && (aCode[1] == 0 || aCode[1] == p->pindx))
        {
          p->pcode = Property::NEUMANN;
          myVectors[p->pindx] = mySol->getScalarSecSol();
          aCode[1] = p->pindx;
        }
        else
        {
          p->pcode = Property::UNDEFINED;
          std::cerr <<" *** SIMPoisson1D::preprocess: Analytic Neumann conditio"
                    <<"ns\n     can only be assigned to one topology set.\n"
                    <<"     Ignoring specification for property code = "
                    << p->pindx << std::endl;
          ok = false;
        }
      }

  if (this->SIM1D::preprocess(ignored,fixDup) && ok) return true;

  std::cerr <<"\n *** SIMPoisson1D::preprocess: Aborting due to above error(s)."
            << std::endl;
  return false;
}


bool SIMPoisson1D::initMaterial (size_t propInd)
{
  if (propInd >= mVec.size()) return false;

  prob.setMaterial(mVec[propInd]);
  return true;
}


bool SIMPoisson1D::initNeumann (size_t propInd)
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
