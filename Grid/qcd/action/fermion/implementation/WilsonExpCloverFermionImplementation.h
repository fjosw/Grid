/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonExpCloverFermion.cc

    Copyright (C) 2017

    Author: paboyle <paboyle@ph.ed.ac.uk>
    Author: Guido Cossu <guido.cossu@ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
/*  END LEGAL */

#include <Grid/Grid.h>
#include <Grid/qcd/spin/Dirac.h>
#include <Grid/qcd/action/fermion/WilsonExpCloverFermion.h>

NAMESPACE_BEGIN(Grid);

// *NOT* EO
template <class Impl>
void WilsonExpCloverFermion<Impl>::M(const FermionField &in, FermionField &out)
{
  FermionField temp(out.Grid());

  // Wilson term
  out.Checkerboard() = in.Checkerboard();
  this->Dhop(in, out, DaggerNo);

  // Clover term
  Mooee(in, temp);

  out += temp;
}

template <class Impl>
void WilsonExpCloverFermion<Impl>::Mdag(const FermionField &in, FermionField &out)
{
  FermionField temp(out.Grid());

  // Wilson term
  out.Checkerboard() = in.Checkerboard();
  this->Dhop(in, out, DaggerYes);

  // Clover term
  MooeeDag(in, temp);

  out += temp;
}

template <class Impl>
void WilsonExpCloverFermion<Impl>::ImportGauge(const GaugeField &_Umu)
{
  WilsonFermion<Impl>::ImportGauge(_Umu);
  GridBase *grid = _Umu.Grid();
  typename Impl::GaugeLinkField Bx(grid), By(grid), Bz(grid), Ex(grid), Ey(grid), Ez(grid);

  // Compute the field strength terms mu>nu
  WilsonLoops<Impl>::FieldStrength(Bx, _Umu, Zdir, Ydir);
  WilsonLoops<Impl>::FieldStrength(By, _Umu, Zdir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Bz, _Umu, Ydir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Ex, _Umu, Tdir, Xdir);
  WilsonLoops<Impl>::FieldStrength(Ey, _Umu, Tdir, Ydir);
  WilsonLoops<Impl>::FieldStrength(Ez, _Umu, Tdir, Zdir);

  // Compute the Clover Operator acting on Colour and Spin
  // multiply here by the clover coefficients for the anisotropy
  CloverTerm  = fillCloverYZ(Bx) * csw_r;
  CloverTerm += fillCloverXZ(By) * csw_r;
  CloverTerm += fillCloverXY(Bz) * csw_r;
  CloverTerm += fillCloverXT(Ex) * csw_t;
  CloverTerm += fillCloverYT(Ey) * csw_t;
  CloverTerm += fillCloverZT(Ez) * csw_t;

  //
  // Make Exp Clover and inverse
  //

  int lvol = _Umu.Grid()->lSites();
  int DimRep = Impl::Dimension;
  {

   typedef iMatrix<ComplexD,6> mat;

   autoView(CTv,CloverTerm,CpuRead);
   autoView(CTExpv,ExpCloverTerm,CpuWrite);
   autoView(CTExpIv,ExpCloverTermInv,CpuWrite);

   thread_for(site, lvol, {
  	Coordinate lcoor;
  	grid->LocalIndexToLocalCoor(site, lcoor);

  	mat srcCloverOpUL(0.0); // upper left block
  	mat srcCloverOpLR(0.0); // lower right block
  	mat ExpCloverOp;

  	typename SiteCloverType::scalar_object Qx = Zero(), Qxexp = Zero();

  	peekLocalSite(Qx, CTv, lcoor);

  	// exp(A)

  	//
  	// upper left block
  	//

  	for (int j = 0; j < Ns/2; j++)
  	 for (int k = 0; k < Ns/2; k++)
  	  for (int a = 0; a < DimRep; a++)
  	   for (int b = 0; b < DimRep; b++){
  	    auto zz =  Qx()(j, k)(a, b);
  	    srcCloverOpUL(a + j * DimRep, b + k * DimRep) = std::complex<double>(zz);
  	   }

  	ExpCloverOp = Exponentiate(srcCloverOpUL,1.0/(diag_mass));

  	for (int j = 0; j < Ns/2; j++)
  	 for (int k = 0; k < Ns/2; k++)
  	  for (int a = 0; a < DimRep; a++)
  	   for (int b = 0; b < DimRep; b++)
  	    Qxexp()(j, k)(a, b) = ExpCloverOp(a + j * DimRep, b + k * DimRep);

  	//
    // lower right block
  	//

  	for (int j = 0; j < Ns/2; j++)
  	 for (int k = 0; k < Ns/2; k++)
  	  for (int a = 0; a < DimRep; a++)
  	   for (int b = 0; b < DimRep; b++){
  	  	auto zz =  Qx()(j+Ns/2, k+Ns/2)(a, b);
  	  	srcCloverOpLR(a + j * DimRep, b + k * DimRep) = std::complex<double>(zz);
  	   }


  	ExpCloverOp = Exponentiate(srcCloverOpLR,1.0/(diag_mass));

  	for (int j = 0; j < Ns/2; j++)
  	 for (int k = 0; k < Ns/2; k++)
  	  for (int a = 0; a < DimRep; a++)
  	   for (int b = 0; b < DimRep; b++)
  	    Qxexp()(j+Ns/2, k+Ns/2)(a, b) = ExpCloverOp(a + j * DimRep, b + k * DimRep);

  	// Now that the full 12x12 block is filled do poke!
  	pokeLocalSite(Qxexp, CTExpv, lcoor);

  	// exp(-A)

  	//
  	// upper left block
  	//

  	ExpCloverOp = Exponentiate(srcCloverOpUL,-1.0/(diag_mass));

  	for (int j = 0; j < Ns/2; j++)
  	 for (int k = 0; k < Ns/2; k++)
  	  for (int a = 0; a < DimRep; a++)
  	   for (int b = 0; b < DimRep; b++)
  	    Qxexp()(j, k)(a, b) = ExpCloverOp(a + j * DimRep, b + k * DimRep);

  	//
    // lower right block
  	//

  	ExpCloverOp = Exponentiate(srcCloverOpLR,-1.0/(diag_mass));

  	for (int j = 0; j < Ns/2; j++)
  	 for (int k = 0; k < Ns/2; k++)
  	  for (int a = 0; a < DimRep; a++)
  	   for (int b = 0; b < DimRep; b++)
  	    Qxexp()(j+Ns/2, k+Ns/2)(a, b) = ExpCloverOp(a + j * DimRep, b + k * DimRep);

  	// Now that the full 12x12 block is filled do poke!
  	pokeLocalSite(Qxexp, CTExpIv, lcoor);

   });
  }
  ExpCloverTerm *= diag_mass;
  ExpCloverTermInv *= (1 / diag_mass);

  //
  //
  //

  CloverTerm += diag_mass;

  //int lvol = _Umu.Grid()->lSites();
  //int DimRep = Impl::Dimension;

  {
    autoView(CTv,CloverTerm,CpuRead);
    autoView(CTIv,CloverTermInv,CpuWrite);
    thread_for(site, lvol, {
      Coordinate lcoor;
      grid->LocalIndexToLocalCoor(site, lcoor);
      Eigen::MatrixXcd EigenCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
      Eigen::MatrixXcd EigenInvCloverOp = Eigen::MatrixXcd::Zero(Ns * DimRep, Ns * DimRep);
      typename SiteCloverType::scalar_object Qx = Zero(), Qxinv = Zero();
      peekLocalSite(Qx, CTv, lcoor);
      //if (csw!=0){
      for (int j = 0; j < Ns; j++)
	for (int k = 0; k < Ns; k++)
	  for (int a = 0; a < DimRep; a++)
	    for (int b = 0; b < DimRep; b++){
	      auto zz =  Qx()(j, k)(a, b);
	      EigenCloverOp(a + j * DimRep, b + k * DimRep) = std::complex<double>(zz);
	    }
      //   if (site==0) std::cout << "site =" << site << "\n" << EigenCloverOp << std::endl;
      
      EigenInvCloverOp = EigenCloverOp.inverse();
      //std::cout << EigenInvCloverOp << std::endl;
      for (int j = 0; j < Ns; j++)
	for (int k = 0; k < Ns; k++)
	  for (int a = 0; a < DimRep; a++)
	    for (int b = 0; b < DimRep; b++)
	      Qxinv()(j, k)(a, b) = EigenInvCloverOp(a + j * DimRep, b + k * DimRep);
      //    if (site==0) std::cout << "site =" << site << "\n" << EigenInvCloverOp << std::endl;
      //  }
      pokeLocalSite(Qxinv, CTIv, lcoor);
    });
  }

  // Separate the even and odd parts
  pickCheckerboard(Even, CloverTermEven, CloverTerm);
  pickCheckerboard(Odd, CloverTermOdd, CloverTerm);

  pickCheckerboard(Even, CloverTermDagEven, adj(CloverTerm));
  pickCheckerboard(Odd, CloverTermDagOdd, adj(CloverTerm));

  pickCheckerboard(Even, CloverTermInvEven, CloverTermInv);
  pickCheckerboard(Odd, CloverTermInvOdd, CloverTermInv);

  pickCheckerboard(Even, CloverTermInvDagEven, adj(CloverTermInv));
  pickCheckerboard(Odd, CloverTermInvDagOdd, adj(CloverTermInv));

  pickCheckerboard(Even, ExpCloverTermEven, ExpCloverTerm);
  pickCheckerboard(Odd, ExpCloverTermOdd, ExpCloverTerm);

  pickCheckerboard(Even, ExpCloverTermDagEven, adj(ExpCloverTerm));
  pickCheckerboard(Odd, ExpCloverTermDagOdd, adj(ExpCloverTerm));

  pickCheckerboard(Even, ExpCloverTermInvEven, ExpCloverTermInv);
  pickCheckerboard(Odd, ExpCloverTermInvOdd, ExpCloverTermInv);

  pickCheckerboard(Even, ExpCloverTermInvDagEven, adj(ExpCloverTermInv));
  pickCheckerboard(Odd, ExpCloverTermInvDagOdd, adj(ExpCloverTermInv));
}

template <class Impl>
void WilsonExpCloverFermion<Impl>::Mooee(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseNo);
}

template <class Impl>
void WilsonExpCloverFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerYes, InverseNo);
}

template <class Impl>
void WilsonExpCloverFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseYes);
}

template <class Impl>
void WilsonExpCloverFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerYes, InverseYes);
}

template <class Impl>
void WilsonExpCloverFermion<Impl>::MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv)
{
  out.Checkerboard() = in.Checkerboard();
  CloverFieldType *Clover;
  assert(in.Checkerboard() == Odd || in.Checkerboard() == Even);

  if (dag)
  {
    if (in.Grid()->_isCheckerBoarded)
    {
      if (in.Checkerboard() == Odd)
      {
        Clover = (inv) ? &ExpCloverTermInvDagOdd : &ExpCloverTermDagOdd;
      }
      else
      {
        Clover = (inv) ? &ExpCloverTermInvDagEven : &ExpCloverTermDagEven;
      }
      out = *Clover * in;
    }
    else
    {
      Clover = (inv) ? &ExpCloverTermInv : &ExpCloverTerm;
      out = adj(*Clover) * in;
    }
  }
  else
  {
    if (in.Grid()->_isCheckerBoarded)
    {

      if (in.Checkerboard() == Odd)
      {
        //  std::cout << "Calling clover term Odd" << std::endl;
        Clover = (inv) ? &ExpCloverTermInvOdd : &ExpCloverTermOdd;
      }
      else
      {
        //  std::cout << "Calling clover term Even" << std::endl;
        Clover = (inv) ? &ExpCloverTermInvEven : &ExpCloverTermEven;
      }
      out = *Clover * in;
      //  std::cout << GridLogMessage << "*Clover.Checkerboard() "  << (*Clover).Checkerboard() << std::endl;
    }
    else
    {
      Clover = (inv) ? &ExpCloverTermInv : &ExpCloverTerm;
      out = *Clover * in;
    }
  }

} // MooeeInternal


// Derivative parts
template <class Impl>
void WilsonExpCloverFermion<Impl>::MooDeriv(GaugeField &mat, const FermionField &X, const FermionField &Y, int dag)
{
  assert(0);
}

// Derivative parts
template <class Impl>
void WilsonExpCloverFermion<Impl>::MeeDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  assert(0); // not implemented yet
}

NAMESPACE_END(Grid);
