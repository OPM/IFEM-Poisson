// $Id$
//==============================================================================
//!
//! \file PoissonSolutions.h
//!
//! \date Jul 1 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Analytic solutions for Poisson problems.
//!
//==============================================================================

#ifndef _POISSON_SOLUTIONS_H
#define _POISSON_SOLUTIONS_H

#include "Function.h"


/*!
  \brief Analytic solution for the Poisson equation on a square domain.
*/

class Square2D : public VecFunc
{
public:
  //! \brief Empty constructor.
  explicit Square2D(double = 1.0) {}
  //! \brief Empty destructor.
  virtual ~Square2D() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief Heat source for the Poisson equation on a square domain.
*/

class Square2DHeat : public RealFunc
{
public:
  //! \brief Empty constructor.
  explicit Square2DHeat(double = 1.0) {}
  //! \brief Empty destructor.
  virtual ~Square2DHeat() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;
};


/*!
  \brief Analytic solution for the Poisson equation on the L-shape domain.
*/

class LshapePoisson : public VecFunc
{
public:
  //! \brief Empty constructor.
  LshapePoisson() {}
  //! \brief Empty destructor.
  virtual ~LshapePoisson() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief Sinusoidal solution of the Poisson equation on a square domain.
*/

class SquareSinus : public VecFunc
{
public:
  //! \brief Empty constructor.
  SquareSinus() {}
  //! \brief Empty destructor.
  virtual ~SquareSinus() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief Heat source for the sinusoidal solution on a square domain.
*/

class SquareSinusSource : public RealFunc
{
public:
  //! \brief Empty constructor.
  SquareSinusSource() {}
  //! \brief Empty destructor.
  virtual ~SquareSinusSource() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;
};


/*!
  \brief Poisson problem with smooth solution and a steep interior layer.
*/

class PoissonInteriorLayer : public VecFunc
{
public:
  //! \brief Default constructor.
  explicit PoissonInteriorLayer(double s = 60.0) : SLOPE(s) {}
  //! \brief Empty destructor.
  virtual ~PoissonInteriorLayer() {}

protected:
  //! \brief Evaluates the heat flux at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;

private:
  double SLOPE; //!< layer SLOPE (large value gives problems in adaptive solver)
};


/*!
  \brief Analytical primary solution for PoissonInteriorLayer.
*/

class PoissonInteriorLayerSol : public RealFunc
{
public:
  //! \brief Default constructor.
  explicit PoissonInteriorLayerSol(double s = 60.0) : SLOPE(s) {}
  //! \brief Empty destructor.
  virtual ~PoissonInteriorLayerSol() {}

protected:
  //! \brief Evaluates the exact temperature distribution at the point \a X.
  virtual double evaluate(const Vec3& X) const;

private:
  double SLOPE; //!< layer SLOPE
};


/*!
  \brief Heat source for PoissonInteriorLayer.
*/

class PoissonInteriorLayerSource : public RealFunc
{
public:
  //! \brief Default constructor.
  explicit PoissonInteriorLayerSource(double s = 60.0) : SLOPE(s) {}
  //! \brief Empty destructor.
  virtual ~PoissonInteriorLayerSource() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;

private:
  double SLOPE; //!< layer SLOPE
};

/*!
  \brief Poisson 3D problem with smooth sharp solution and a steep interior layer.
*/

class PoissonWaterfall : public VecFunc
{
public:
  //! \brief Default constructor.
  explicit PoissonWaterfall(double eps = 0.002) : epsilon(eps) {}
  //! \brief Empty destructor.
  virtual ~PoissonWaterfall() {}

protected:
  //! \brief Evaluates the heat flux at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;

private:
  double epsilon; //!< layer epsilon (waterfall width approx sqrt(5eps) )
};


/*!
  \brief Analytical primary solution for PoissonWaterfall.
*/

class PoissonWaterfallSol : public RealFunc
{
public:
  //! \brief Default constructor.
  explicit PoissonWaterfallSol(double eps = 0.002) : epsilon(eps) {}
  //! \brief Empty destructor.
  virtual ~PoissonWaterfallSol() {}

protected:
  //! \brief Evaluates the exact temperature distribution at the point \a X.
  virtual double evaluate(const Vec3& X) const;

private:
  double epsilon; //!< layer epsilon (waterfall width approx sqrt(5eps) )
};


/*!
  \brief Heat source for PoissonWaterfall.
*/

class PoissonWaterfallSource : public RealFunc
{
public:
  //! \brief Default constructor.
  explicit PoissonWaterfallSource(double eps = 0.002) : epsilon(eps) {}
  //! \brief Empty destructor.
  virtual ~PoissonWaterfallSource() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;

private:
  double epsilon; //!< layer epsilon (waterfall width approx sqrt(5eps) )
};


/*!
  \brief Analytic solution for the Poisson equation on a cube domain.
*/

class PoissonCube : public VecFunc
{
public:
  //! \brief Empty Constructor.
  PoissonCube() {}
  //! \brief Empty destructor.
  virtual ~PoissonCube() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;
};


/*!
  \brief Heat source for the Poisson equation on a cube domain.
*/

class PoissonCubeSource : public RealFunc
{
public:
  //! \brief Empty constructor.
  PoissonCubeSource() {}
  //! \brief Empty destructor.
  virtual ~PoissonCubeSource() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;
};


/*!
  \brief Analytic solution for the Poisson equation on a line domain.
*/

class PoissonLine : public VecFunc
{
public:
  //! \brief Default constructor.
  explicit PoissonLine(double r = 1.0) : L(r) {}
  //! \brief Empty destructor.
  virtual ~PoissonLine() {}

protected:
  //! \brief Evaluates the analytic flux vector at the point \a X.
  virtual Vec3 evaluate(const Vec3& X) const;

private:
  double L; //!< Length parameter
};


/*!
  \brief Heat source for the Poisson equation on a line.
*/

class PoissonLineSource : public RealFunc
{
public:
  //! \brief Default constructor.
  explicit PoissonLineSource(double r = 1.0) : L(r) {}
  //! \brief Empty destructor.
  virtual ~PoissonLineSource() {}

protected:
  //! \brief Evaluates the heat field at the point \a X.
  virtual double evaluate(const Vec3& X) const;

private:
  double L; //!< Length parameter
};

#endif
