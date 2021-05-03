/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/WilsonCloverTMFermion.cc

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
#include <Grid/qcd/action/fermion/WilsonCloverTMFermion.h>

NAMESPACE_BEGIN(Grid);

// *NOT* EO
template <class Impl>
void WilsonCloverTMFermion<Impl>::M(const FermionField &in, FermionField &out)
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
void WilsonCloverTMFermion<Impl>::Mdag(const FermionField &in, FermionField &out)
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
void WilsonCloverTMFermion<Impl>::ImportGauge(const GaugeField &_Umu)
{
  WilsonFermion<Impl>::ImportGauge(_Umu);
  GridBase *grid = _Umu.Grid();
  Gamma G5(Gamma::Algebra::Gamma5);
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
  CloverTerm += diag_mass;

  CloverFieldType T(CloverTerm.Grid());
  T = Zero();
  autoView(T_v,T,AcceleratorWrite);
  accelerator_for(i, CloverTerm.Grid()->oSites(),1,
  {
    T_v[i]()(0, 0) = +mu;
    T_v[i]()(1, 1) = +mu;
    T_v[i]()(2, 2) = -mu;
    T_v[i]()(3, 3) = -mu;
  });
  T = timesI(T);
  CloverTerm += T;

  int lvol = _Umu.Grid()->lSites();
  int DimRep = Impl::Dimension;

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
}

template <class Impl>
void WilsonCloverTMFermion<Impl>::Mooee(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseNo);
}

template <class Impl>
void WilsonCloverTMFermion<Impl>::MooeeDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerYes, InverseNo);
}

template <class Impl>
void WilsonCloverTMFermion<Impl>::MooeeInv(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerNo, InverseYes);
}

template <class Impl>
void WilsonCloverTMFermion<Impl>::MooeeInvDag(const FermionField &in, FermionField &out)
{
  this->MooeeInternal(in, out, DaggerYes, InverseYes);
}

template <class Impl>
void WilsonCloverTMFermion<Impl>::MooeeInternal(const FermionField &in, FermionField &out, int dag, int inv)
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
        Clover = (inv) ? &CloverTermInvDagOdd : &CloverTermDagOdd;
      }
      else
      {
        Clover = (inv) ? &CloverTermInvDagEven : &CloverTermDagEven;
      }
      out = *Clover * in;
    }
    else
    {
      Clover = (inv) ? &CloverTermInv : &CloverTerm;
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
        Clover = (inv) ? &CloverTermInvOdd : &CloverTermOdd;
      }
      else
      {
        //  std::cout << "Calling clover term Even" << std::endl;
        Clover = (inv) ? &CloverTermInvEven : &CloverTermEven;
      }
      out = *Clover * in;
      //  std::cout << GridLogMessage << "*Clover.Checkerboard() "  << (*Clover).Checkerboard() << std::endl;
    }
    else
    {
      Clover = (inv) ? &CloverTermInv : &CloverTerm;
      out = *Clover * in;
    }
  }

} // MooeeInternal


// Derivative parts
template <class Impl>
void WilsonCloverTMFermion<Impl>::MooDeriv(GaugeField &mat, const FermionField &X, const FermionField &Y, int dag)
{
  assert(0);
}

// Derivative parts
template <class Impl>
void WilsonCloverTMFermion<Impl>::MeeDeriv(GaugeField &mat, const FermionField &U, const FermionField &V, int dag)
{
  assert(0); // not implemented yet
}

NAMESPACE_END(Grid);
