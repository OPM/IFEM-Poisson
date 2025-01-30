// $Id$
//==============================================================================
//!
//! \file PoissonSource.h
//!
//! \date Jan 29 2025
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Class for Poisson source function.
//!
//==============================================================================

#ifndef POISSON_SOURCE_H_
#define POISSON_SOURCE_H_

#include "Function.h"

class AnaSol;
class Poisson;


/*!
  \brief Class that derives the Poisson source function from the analytic solution.
 */

class PoissonAnaSolSource : public RealFunc
{
public:
  //! \brief Constructor for constant kappa
  //! \param aSol Analytic solution to use
  //! \param prob Reference to problem integrand (for material properties)
  PoissonAnaSolSource(const AnaSol& aSol, const Poisson& prob);

protected:
  //! \brief Evaluates the function.
  double evaluate(const Vec3& X) const override;

  const AnaSol& anaSol; //!< Reference to analytic solution
  const Poisson& poisson; //!< Reference to integrand
};

#endif
