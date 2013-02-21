// $Id$
//==============================================================================
//!
//! \file SIMPoisson.h
//!
//! \date May 25 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for NURBS-based FEM analysis of the Poisson equation.
//!
//==============================================================================

#ifndef _SIM_POISSON_H
#define _SIM_POISSON_H

#include "Poisson.h"
#include "AnaSol.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "Utilities.h"
#include "tinyxml.h"


/*!
  \brief Driver class for NURBS-based FEM analysis of Poisson problems.
*/

template<class Dim> class SIMPoisson : public Dim
{
public:
  //! \brief Default constructor.
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  SIMPoisson(bool checkRHS = false) : Dim(1,0,checkRHS), prob(Dim::dimension)
  {
    Dim::myProblem = &prob;
    aCode[0] = aCode[1] = 0;
  }

  //! \brief The destructor cleans up the dynamically allocated objects.
  virtual ~SIMPoisson()
  {
    Dim::myProblem = NULL; // Because it is not dynamically allocated

    // To prevent the SIMbase destructor try to delete already deleted functions
    if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
    if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
  }

  //! \brief Initializes the property containers of the model.
  virtual void clearProperties()
  {
    // To prevent SIMbase::clearProperties deleting the analytical solution
    if (aCode[0] > 0) Dim::myScalars.erase(aCode[0]);
    if (aCode[1] > 0) Dim::myVectors.erase(aCode[1]);
    aCode[0] = aCode[1] = 0;

    mVec.clear();
    prob.setSource(NULL);
    prob.setTraction((RealFunc*)NULL);
    prob.setTraction((VecFunc*)NULL);
    this->SIMbase::clearProperties();
  }

  //! \brief Performs some pre-processing tasks on the FE model.
  //! \details This method is reimplemented to resolve inhomogeneous boundary
  //! condition fields in case they are derived from the analytical solution.
  virtual bool preprocess(const std::vector<int>& ignored, bool fixDup)
  {
    if (Dim::mySol) // Define analytical boundary condition fields
      for (PropertyVec::iterator p = Dim::myProps.begin(); p != Dim::myProps.end(); p++)
        if (p->pcode == Property::DIRICHLET_ANASOL)
        {
          if (!Dim::mySol->getScalarSol())
            p->pcode = Property::UNDEFINED;
          else if (aCode[0] == abs(p->pindx))
            p->pcode = Property::DIRICHLET_INHOM;
          else if (aCode[0] == 0)
          {
            aCode[0] = abs(p->pindx);
            Dim::myScalars[aCode[0]] = Dim::mySol->getScalarSol();
            p->pcode = Property::DIRICHLET_INHOM;
          }
          else
            p->pcode = Property::UNDEFINED;
        }
        else if (p->pcode == Property::NEUMANN_ANASOL)
        {
          if (!Dim::mySol->getScalarSecSol())
            p->pcode = Property::UNDEFINED;
          else if (aCode[1] == p->pindx)
            p->pcode = Property::NEUMANN;
          else if (aCode[1] == 0)
          {
            aCode[1] = p->pindx;
            Dim::myVectors[aCode[1]] = Dim::mySol->getScalarSecSol();
            p->pcode = Property::NEUMANN;
          }
          else
            p->pcode = Property::UNDEFINED;
        }

    return Dim::preprocess(ignored,fixDup);
  }

protected:
  //! \brief Parses a data section from the input stream.
  //! \param[in] keyWord Keyword of current data section to read
  //! \param is The file stream to read from
  virtual bool parse(char* keyWord, std::istream& is)
  {
    char* cline = 0;

    if (parseDimSpecific(keyWord, is))
      return true;

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

    else
      return Dim::parse(keyWord,is);

    return true;
  }

  //! \brief Parses a data section from an XML element.
  //! \param[in] elem The XML element to parse
  virtual bool parse(const TiXmlElement* elem)
  {
    if (strcasecmp(elem->Value(),"poisson"))
      return Dim::parse(elem);

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement()) {
      if (parseDimSpecific(child))
        continue;

      if (!strcasecmp(child->Value(),"isotropic")) {
        int code = this->parseMaterialSet(child,mVec.size());
        double kappa = 1000.0;
        utl::getAttribute(child,"kappa",kappa);
        if (code == 0)
          prob.setMaterial(kappa);
        mVec.push_back(kappa);
        std::cout <<"\tMaterial code "<< code <<": "<< kappa << std::endl;
      } else
        Dim::parse(child);
    }

    return true;
  }

  //! \brief Initializes material properties for integration of interior terms.
  //! \param[in] propInd Physical property index
  virtual bool initMaterial(size_t propInd)
  {
    if (propInd >= mVec.size()) return false;

    prob.setMaterial(mVec[propInd]);
    return true;
  }

  //! \brief Initializes for integration of Neumann terms for a given property.
  //! \param[in] propInd Physical property index
  virtual bool initNeumann(size_t propInd)
  {
    typename Dim::SclFuncMap::const_iterator sit = Dim::myScalars.find(propInd);
    typename Dim::VecFuncMap::const_iterator vit = Dim::myVectors.find(propInd);

    if (sit != Dim::myScalars.end())
      prob.setTraction(sit->second);
    else if (vit != Dim::myVectors.end())
      prob.setTraction(vit->second);
    else
      return false;

    return true;
  }

private:
  //! \details This function allows for specialization of the template
  //!          while still reusing as much code as possible. Only put
  //!          dimension-specific code in here
  bool parseDimSpecific(char* keyWord, std::istream& is);
  bool parseDimSpecific(const TiXmlElement* child);

  Poisson   prob;     //!< Data and methods for the Poisson problem
  RealArray mVec;     //!< Material data
  int       aCode[2]; //!< Analytical BC code (used by destructor)
};


typedef SIMPoisson<SIM1D> SIMPoisson1D; //!< 1D specific driver
typedef SIMPoisson<SIM2D> SIMPoisson2D; //!< 2D specific driver
typedef SIMPoisson<SIM3D> SIMPoisson3D; //!< 3D specific driver

//! \brief Template specialization - 1D specific input parsing.
template<> bool SIMPoisson1D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 1D specific input parsing.
template<> bool SIMPoisson1D::parseDimSpecific(const TiXmlElement* elem);

//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMPoisson2D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 2D specific input parsing.
template<> bool SIMPoisson2D::parseDimSpecific(const TiXmlElement* elem);

//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMPoisson3D::parseDimSpecific(char* keyWord, std::istream& is);
//! \brief Template specialization - 3D specific input parsing.
template<> bool SIMPoisson3D::parseDimSpecific(const TiXmlElement* elem);

#endif
