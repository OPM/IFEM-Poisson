// $Id$
//==============================================================================
//!
//! \file PoissonSource.C
//!
//! \date Jan 29 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Poisson source function.
//!
//==============================================================================

#include "Poisson.h"
#include "PoissonSource.h"

#include "AnaSol.h"


PoissonAnaSolSource::PoissonAnaSolSource (const AnaSol& aSol,
                                          const Poisson& prob)
  : anaSol(aSol), poisson(prob)
{
}


double PoissonAnaSolSource::evaluate (const Vec3& X) const
{
  return -poisson.getMaterial(X) * anaSol.getScalarSol()->hessian(X).trace();
}
