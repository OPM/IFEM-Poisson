// $Id$
//==============================================================================
//!
//! \file PoissonArgs.C
//!
//! \date Jan 24 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the Poisson application.
//!
//==============================================================================

#include "PoissonArgs.h"
#include "Utilities.h"
#include "tinyxml.h"


bool PoissonArgs::parse(const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"poisson"))
    utl::getAttribute(elem,"adaptive",adap);

  return SIM::AppXMLInputBase::parse(elem);
}
