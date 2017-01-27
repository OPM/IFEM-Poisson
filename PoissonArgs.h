// $Id$
//==============================================================================
//!
//! \file PoissonArgs.h
//!
//! \date Jan 24 2017
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Preparsing of input files for the Poisson application.
//!
//==============================================================================

#include "AppCommon.h"

/*! \brief Struct holding application parameters.
 */

class PoissonArgs : public SIM::AppXMLInputBase
{
public:
  bool adap = false; //!< True to run an adaptive simulator

protected:
  //! \brief Parse an element from the input file
  bool parse(const TiXmlElement* elem) override;
};
