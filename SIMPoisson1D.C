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
#include "Utilities.h"
#include "AnaSol.h"
#include <string.h>
#include "Functions.h"

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
      std::cout <<"\tMaterial code "<< code <<":"<< kappa << std::endl;
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
    cline = strtok(keyWord+6," ");
    int code = -1;
    if (!strncasecmp(cline,"LINE",4))
    {
      double L = atof(strtok(NULL," "));
      std::cout <<"\nAnalytical solution: Line L="<< L << std::endl;
      mySol = new AnaSol(NULL,new PoissonLine(L));
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
    return SIM1D::parse(elem);

  std::vector<const TiXmlElement*> parsed = handlePriorityTags(elem);
  const TiXmlElement* child = elem->FirstChildElement();
  while (child) {
    if (find(parsed.begin(),parsed.end(),child) != parsed.end()) {
      child = child->NextSiblingElement();
      continue;
    }
    if (!strcasecmp(child->Value(),"isotropic")) {
      int code = 0;
      if (child->Attribute("code"))
        code = atoi(child->Attribute("code"));
      double kappa=1000.f;
      if (child->Attribute("kappa"))
        kappa = atof(child->Attribute("kappa"));
      if (code == 0)
        prob.setMaterial(kappa);
      else if (setPropertyType(code,Property::MATERIAL,mVec.size()))
        mVec.push_back(kappa);
    }
    else if (!strcasecmp(child->Value(),"source")) {
      if (child->Attribute("type") &&
          !strcasecmp(child->Attribute("type"),"line")) {
        double L=0;
        if (child->Attribute("L"))
          L = atof(child->Attribute("L"));
        std::cout <<"\nHeat source function: Line L="<< L << std::endl;
        prob.setSource(new PoissonLineSource(L));
      } else if (child->Attribute("type") &&
          !strcasecmp(child->Attribute("type"),"expression")) {
        if (child->FirstChild() && child->FirstChild()->Value()) {
          std::cout << "\nHeat source function: " 
                    << child->FirstChild()->Value() << std::endl;
          prob.setSource(new EvalFunction(child->FirstChild()->Value()));
        }
      } else
        std::cerr <<"  ** SIMPoisson1D::parse: Unknown source function "
                  << (child->Attribute("type")?child->Attribute("type"):"") << std::endl;
    } else if (!strcasecmp(child->Value(),"anasol")) {
      if (child->Attribute("type") &&
          !strcasecmp(child->Attribute("type"),"line")) {
        double L=0;
        if (child->Attribute("L"))
          L = atof(child->Attribute("L"));
        std::cout <<"\nAnalytical solution: Line L="<< L << std::endl;
        mySol = new AnaSol(NULL,new PoissonLine(L));
      } else if (child->Attribute("type") &&
                 !strcasecmp(child->Attribute("type"),"expression")) {
        std::string variables, primary, secondary;
        const TiXmlElement* var = child->FirstChildElement("variables");
        if (var && var->FirstChild() && var->FirstChild()->Value()) {
          variables = var->FirstChild()->Value();
          if (variables[variables.size()-1] != ';')
            variables += ";";
        }
        const TiXmlElement* prim = child->FirstChildElement("primary");
        RealFunc* s = NULL;
        VecFunc* v = NULL;
        if (prim && prim->FirstChild() && prim->FirstChild()->Value()) {
          primary = prim->FirstChild()->Value();
          s = new EvalFunction((variables+primary).c_str());
        }
        const TiXmlElement* sec = child->FirstChildElement("secondary");
        if (sec && sec->FirstChild() && sec->FirstChild()->Value()) {
          secondary = sec->FirstChild()->Value();
          v = new VecFuncExpr(secondary,variables);
        }
        std::cout <<"\nAnalytical solution:" << std::endl;
        if (!variables.empty())
          std::cout << "\t Variables=" << variables << std::endl;
        if (s)
          std::cout << "\t Primary=" << primary << std::endl;
        if (v)
          std::cout << "\t Secondary=" << secondary << std::endl;
        mySol = new AnaSol(s, v);
      } else
        std::cerr <<"  ** SIMPoisson1D::parse: Unknown analytical solution "
          << (child->Attribute("type")?child->Attribute("type"):"") << std::endl;

        // Define the analytical boundary traction field
        size_t code=0;
        if (child->Attribute("code"))
          code = atoi(child->Attribute("code"));
        if (code > 0 && mySol && mySol->getScalarSecSol())
        {
          setPropertyType(code,Property::NEUMANN);
          myVectors[code] = mySol->getScalarSecSol();
          myAFcode = code;
        }
      }
      else
        SIM1D::parse(child);

      child = child->NextSiblingElement();
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
