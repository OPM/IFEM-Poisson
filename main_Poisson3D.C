// $Id$
//==============================================================================
//!
//! \file main_Poisson3D.C
//!
//! \date 20 May 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Main program for the isogeometric solver for the Poisson equation.
//!
//==============================================================================

#include "SIMPoisson3D.h"
#include "SIMPoisson2D.h"
#include "SIMPoisson1D.h"
#include "AdaptiveSIM.h"
#include "LinAlgInit.h"
#include "HDF5Writer.h"
#include "XMLWriter.h"
#include "Utilities.h"
#include "Profiler.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>


/*!
  \brief Main program for the NURBS-based isogeometric Poisson equation solver.

  The input to the program is specified through the following
  command-line arguments. The arguments may be given in arbitrary order.

  \arg \a input-file : Input file with model definition
  \arg -dense :   Use the dense LAPACK matrix equation solver
  \arg -spr :     Use the SPR direct equation solver
  \arg -superlu : Use the sparse SuperLU equation solver
  \arg -samg :    Use the sparse algebraic multi-grid equation solver
  \arg -petsc :   Use equation solver from PETSc library
  \arg -lag : Use Lagrangian basis functions instead of splines/NURBS
  \arg -spec : Use Spectral basis functions instead of splines/NURBS
  \arg -LR : Use LR-spline basis functions instead of tensorial splines/NURBS
  \arg -nGauss \a n : Number of Gauss points over a knot-span in each direction
  \arg -vtf \a format : VTF-file format (-1=NONE, 0=ASCII, 1=BINARY)
  \arg -nviz \a nviz : Number of visualization points over each knot-span
  \arg -nu \a nu : Number of visualization points per knot-span in u-direction
  \arg -nv \a nv : Number of visualization points per knot-span in v-direction
  \arg -nw \a nw : Number of visualization points per knot-span in w-direction
  \arg -hdf5 : Write primary and projected secondary solution to HDF5 file
  \arg -ignore \a p1, \a p2, ... : Ignore these patches in the analysis
  \arg -eig \a iop : Eigenproblem solver to use (1...6)
  \arg -nev \a nev : Number of eigenvalues to compute
  \arg -ncv \a ncv : Number of Arnoldi vectors to use in the eigenvalue analysis
  \arg -shift \a shf : Shift value to use in the eigenproblem solver
  \arg -free : Ignore all boundary conditions (use in free vibration analysis)
  \arg -check : Data check only, read model and output to VTF (no solution)
  \arg -checkRHS : Check that the patches are modelled in a right-hand system
  \arg -vizRHS : Save the right-hand-side load vector on the VTF-file
  \arg -fixDup : Resolve co-located nodes by merging them into a single node
  \arg -1D : Use one-parametric simulation driver
  \arg -2D : Use two-parametric simulation driver
  \arg -adap : Use adaptive simulation driver
*/

int main (int argc, char** argv)
{
  Profiler prof(argv[0]);
  utl::profiler->start("Initialization");

  SystemMatrix::Type solver = SystemMatrix::SPARSE;
  int nGauss = 4;
  int format = -1;
  int n[3] = { 2, 2, 2 };
  std::vector<int> ignoredPatches;
  int iop = 0;
  int nev = 10;
  int ncv = 20;
  double shf = 0.0;
  bool checkRHS = false;
  bool vizRHS = false;
  bool fixDup = false;
  bool dumpHDF5 = false;
  char ndim = 3;
  char* infile = 0;

  LinAlgInit& linalg = LinAlgInit::Init(argc,argv);

  for (int i = 1; i < argc; i++)
    if (!strcmp(argv[i],"-dense"))
      solver = SystemMatrix::DENSE;
    else if (!strcmp(argv[i],"-spr"))
      solver = SystemMatrix::SPR;
    else if (!strcmp(argv[i],"-superlu"))
      solver = SystemMatrix::SPARSE;
    else if (!strcmp(argv[i],"-samg"))
      solver = SystemMatrix::SAMG;
    else if (!strcmp(argv[i],"-petsc"))
      solver = SystemMatrix::PETSC;
    else if (!strcmp(argv[i],"-nGauss") && i < argc-1)
      nGauss = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-vtf") && i < argc-1)
      format = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-hdf5"))
      dumpHDF5 = true;
    else if (!strcmp(argv[i],"-nviz") && i < argc-1)
      n[0] = n[1] = n[2] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nu") && i < argc-1)
      n[0] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nv") && i < argc-1)
      n[1] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nw") && i < argc-1)
      n[2] = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ignore"))
      while (i < argc-1 && isdigit(argv[i+1][0]))
	utl::parseIntegers(ignoredPatches,argv[++i]);
    else if (!strcmp(argv[i],"-eig") && i < argc-1)
      iop = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-nev") && i < argc-1)
      nev = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-ncv") && i < argc-1)
      ncv = atoi(argv[++i]);
    else if (!strcmp(argv[i],"-shift") && i < argc-1)
      shf = atof(argv[++i]);
    else if (!strcmp(argv[i],"-free"))
      SIMbase::ignoreDirichlet = true;
    else if (!strcmp(argv[i],"-check"))
      iop = 100;
    else if (!strcmp(argv[i],"-checkRHS"))
      checkRHS = true;
    else if (!strcmp(argv[i],"-vizRHS"))
      vizRHS = true;
    else if (!strcmp(argv[i],"-fixDup"))
      fixDup = true;
    else if (!strcmp(argv[i],"-1D"))
      ndim = 1;
    else if (!strcmp(argv[i],"-2D"))
      ndim = 2;
    else if (!strncmp(argv[i],"-adap",5))
      iop = 10;
    else if (!strncmp(argv[i],"-lag",4))
      SIMbase::discretization = ASM::Lagrange;
    else if (!strncmp(argv[i],"-spec",5))
      SIMbase::discretization = ASM::Spectral;
    else if (!strncmp(argv[i],"-LR",3))
      SIMbase::discretization = ASM::LRSpline;
    else if (!infile)
      infile = argv[i];
    else
      std::cerr <<"  ** Unknown option ignored: "<< argv[i] << std::endl;

  if (!infile)
  {
    std::cout <<"usage: "<< argv[0]
	      <<" <inputfile> [-dense|-spr|-superlu|-samg|-petsc]\n      "
	      <<" [-free] [-lag|-spec|-LR] [-1D|-2D] [-adap] [-nGauss <n>]\n"
	      <<"       [-vtf <format>] [-nviz <nviz>]"
	      <<" [-nu <nu>] [-nv <nv>] [-nw <nw>] [-hdf5]\n"
	      <<"       [-eig <iop>] [-nev <nev>] [-ncv <ncv] [-shift <shf>]\n"
	      <<"       [-ignore <p1> <p2> ...] [-fixDup]"
	      <<" [-checkRHS] [-check]\n";
    return 0;
  }

  // Load vector visualization is not available when using additional viz-points
  for (int d = 0; d < 3; d++)
    if (d >= ndim)
      n[d] = 1;
    else if (n[d] > 2)
      vizRHS = false;

  // Boundary conditions can be ignored only in generalized eigenvalue analysis
  if (iop != 4 && iop != 6) SIMbase::ignoreDirichlet = false;

  if (linalg.myPid == 0)
  {
    std::cout <<"\n >>> IFEM Poisson equation solver <<<"
	      <<"\n ==========================================\n"
	      <<"\nInput file: "<< infile
	      <<"\nEquation solver: "<< solver
	      <<"\nNumber of Gauss points: "<< nGauss;
    if (format >= 0)
    {
      std::cout <<"\nVTF file format: "<< (format ? "BINARY":"ASCII")
		<<"\nNumber of visualization points: "<< n[0];
      if (ndim > 1) std::cout <<" "<< n[1];
      if (ndim > 2) std::cout <<" "<< n[2];
    }

    if (iop > 0 && iop < 10)
      std::cout <<"\nEigenproblem solver: "<< iop
		<<"\nNumber of eigenvalues: "<< nev
		<<"\nNumber of Arnoldi vectors: "<< ncv
		<<"\nShift value: "<< shf;
    if (SIMbase::discretization == ASM::Lagrange)
      std::cout <<"\nLagrangian basis functions are used";
    else if (SIMbase::discretization == ASM::Spectral)
      std::cout <<"\nSpectral basis functions are used";
    else if (SIMbase::discretization == ASM::LRSpline)
      std::cout <<"\nLR-spline basis functions are used";
    if (SIMbase::ignoreDirichlet)
      std::cout <<"\nSpecified boundary conditions are ignored";
    if (fixDup)
      std::cout <<"\nCo-located nodes will be merged";
    if (checkRHS)
      std::cout <<"\nCheck that each patch has a right-hand coordinate system";
    if (!ignoredPatches.empty())
    {
      std::cout <<"\nIgnored patches:";
      for (size_t i = 0; i < ignoredPatches.size(); i++)
	std::cout <<" "<< ignoredPatches[i];
    }
    std::cout << std::endl;
  }
  utl::profiler->stop("Initialization");
  utl::profiler->start("Model input");

  // Read in model definitions and establish the FE data structures
  SIMbase* model;
  if (ndim == 1)
    model = new SIMPoisson1D();
  else if (ndim == 2)
    model = new SIMPoisson2D();
  else
    model = new SIMPoisson3D(checkRHS);

  SIMinput* theSim = model;
  AdaptiveSIM* aSim = 0;
  if (iop == 10)
    theSim = aSim = new AdaptiveSIM(model);

  if (!theSim->read(infile) || !model->preprocess(ignoredPatches,fixDup))
    return 1;

  utl::profiler->stop("Model input");

  model->setQuadratureRule(nGauss);

  Matrix eNorm;
  Vector gNorm, sol, load;
  std::vector<Mode> modes;
  int iStep = 1, nBlock = 0;

  double RMS_norm = 0;
  double avg_norm = 0;
  double min_err = 10000000;
  double max_err  = 0;

  DataExporter* exporter=NULL;
  if (dumpHDF5)
  {
    std::string foo = infile;
    size_t pos = foo.find_last_of('.');
    foo = foo.substr(0,pos);
    if (linalg.myPid == 0)
      std::cout <<"\nWriting HDF5 file "<< foo <<".hdf5"<< std::endl;
    exporter = new DataExporter(true);
    exporter->registerField("u","heat",DataExporter::SIM,
			   DataExporter::PRIMARY);
    exporter->setFieldValue("u",model,&sol);
    exporter->registerWriter(new HDF5Writer(foo));
    exporter->registerWriter(new XMLWriter(foo));
  }

  switch (iop) {
  case 0:
    model->setMode(SIM::STATIC);
    model->initSystem(solver,1,1);
    model->setAssociatedRHS(0,0);
    if (!model->assembleSystem())
      return 2;
    else if (vizRHS)
      model->extractLoadVec(load);

    // Solve the linear system of equations
    if (!model->solveSystem(sol,1))
      return 3;

    // Evaluate solution norms
    model->setMode(SIM::RECOVERY);
    if (!model->solutionNorms(Vectors(1,sol),eNorm,gNorm))
      return 4;

    // hard coding in the evaluation of the root-mean square of the error
    for(size_t i=1; i<=eNorm.cols(); i++) {
      avg_norm += eNorm(4,i);
      min_err = (min_err < eNorm(4,i)) ? min_err : eNorm(4,i);
      max_err = (max_err > eNorm(4,i)) ? max_err : eNorm(4,i);
    }
    avg_norm  /= eNorm.cols();
    for(size_t i=1; i<=eNorm.cols(); i++)
      RMS_norm += pow(eNorm(4,i)-avg_norm,2);
    RMS_norm = sqrt(RMS_norm/eNorm.cols())/avg_norm;
    gNorm.push_back(RMS_norm);
    gNorm.push_back(min_err);
    gNorm.push_back(max_err);
    gNorm.push_back(avg_norm);

    if (linalg.myPid == 0)
      AdaptiveSIM::printNorms(gNorm,std::cout);
    break;


  case 10:
    // Adaptive simulation
    while (true) {
      char iterationTag[256];
      sprintf(iterationTag, "Adaptive step #%03d", iStep);
      utl::profiler->start(iterationTag);
      if (!aSim->solveStep(infile,solver,iStep))
	return 5;
      else if (!aSim->writeGlv(infile,format,n,iStep,nBlock))
	return 6;
      else if (dumpHDF5)
        exporter->dumpTimeLevel(NULL,true);
      utl::profiler->stop(iterationTag);
      sprintf(iterationTag, "Refinement step #%03d", iStep);
      utl::profiler->start(iterationTag);
      if (!aSim->adaptMesh(++iStep))
	break;
      utl::profiler->stop(iterationTag);
    }

  case 100:
    break; // Model check

  default:
    // Free vibration: Assemble coefficient matrix [K]
    model->setMode(SIM::VIBRATION);
    model->initSystem(solver,1,0);
    if (!model->assembleSystem())
      return 5;

    if (!model->systemModes(modes,nev,ncv,iop,shf))
      return 6;
  }

  utl::profiler->start("Postprocessing");

  if (iop != 10 && format >= 0)
  {
    // Write VTF-file with model geometry
    if (!model->writeGlvG(n,nBlock,infile,format))
      return 7;

    // Write boundary tractions, if any
    if (!model->writeGlvT(iStep,nBlock))
      return 8;

    // Write Dirichlet boundary conditions
    if (!model->writeGlvBC(n,nBlock))
      return 8;

    // Write load vector to VTF-file
    if (!model->writeGlvV(load,"Load vector",n,iStep,nBlock))
      return 9;

    // Write solution fields to VTF-file
    if (!model->writeGlvS(sol,n,iStep,nBlock))
      return 10;

    // Write eigenmodes
    for (size_t j = 0; j < modes.size(); j++)
      if (!model->writeGlvM(modes[j], iop==3 || iop==4 || iop==6, n, nBlock))
	return 11;

    // Write element norms
    if (!model->writeGlvN(eNorm,iStep,nBlock))
      return 12;

    model->writeGlvStep(1,0.0,1);
  }
  model->closeGlv();
  if (dumpHDF5 && iop != 10)
    exporter->dumpTimeLevel();

  utl::profiler->stop("Postprocessing");
  delete theSim;
  delete exporter;
  return 0;
}
