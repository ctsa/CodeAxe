// -*- mode: c++; indent-tabs-mode: nil; -*-
//
//
// CodeAxe : phylogenetic analysis and simulation tools
//
//   http://www.phrap.org
//
//
// Copyright 2007 Christopher T Saunders (ctsa@u.washington.edu)
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
//
//

// $Id: nsaa_maps.cc 743 2007-08-14 15:47:12Z ctsa $

/// \file

#include "nsaa_maps.h"

const char NSAA_MAPS::LEHNINGER::syms[] = {'L','R','p','P','N'};
const char NSAA_MAPS::LANGE::syms[] = {'H','P','N','C'};
const char NSAA_MAPS::HP::syms[] = {'H','P'};
const char NSAA_MAPS::HPO::syms[] = {'H','P','O'};
const char NSAA_MAPS::BETA_BRANCH::syms[] = {'O','B'};
const char NSAA_MAPS::HELIX_FORM::syms[] = {'O','H'};
const char NSAA_MAPS::AROMATIC::syms[] = {'O','A'};
const char NSAA_MAPS::CYS::syms[] = {'O','C'};
const char NSAA_MAPS::GLY::syms[] = {'O','G'};
const char NSAA_MAPS::PRO::syms[] = {'O','P'};


const unsigned*
NSAA_MAPS::LEHNINGER::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    reduction_map[NSAA::G] = ALIPHATIC;
    reduction_map[NSAA::A] = ALIPHATIC;
    reduction_map[NSAA::V] = ALIPHATIC;
    reduction_map[NSAA::L] = ALIPHATIC;
    reduction_map[NSAA::I] = ALIPHATIC;
    reduction_map[NSAA::M] = ALIPHATIC;

    reduction_map[NSAA::F] = AROMATIC;
    reduction_map[NSAA::Y] = AROMATIC;
    reduction_map[NSAA::W] = AROMATIC;

    reduction_map[NSAA::S] = POLAR;
    reduction_map[NSAA::T] = POLAR;
    reduction_map[NSAA::C] = POLAR;
    reduction_map[NSAA::P] = POLAR;
    reduction_map[NSAA::N] = POLAR;
    reduction_map[NSAA::Q] = POLAR;

    reduction_map[NSAA::H] = POSITIVE;
    reduction_map[NSAA::K] = POSITIVE;
    reduction_map[NSAA::R] = POSITIVE;

    reduction_map[NSAA::E] = NEGATIVE;
    reduction_map[NSAA::D] = NEGATIVE;

    is_init = true;
  }
  return reduction_map;
}



const unsigned*
NSAA_MAPS::LANGE::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    reduction_map[NSAA::G] = NONPOLAR;
    reduction_map[NSAA::I] = NONPOLAR;
    reduction_map[NSAA::V] = NONPOLAR;
    reduction_map[NSAA::L] = NONPOLAR;
    reduction_map[NSAA::A] = NONPOLAR;
    reduction_map[NSAA::M] = NONPOLAR;
    reduction_map[NSAA::P] = NONPOLAR;
    reduction_map[NSAA::F] = NONPOLAR;
    reduction_map[NSAA::W] = NONPOLAR;

    reduction_map[NSAA::Q] = POSITIVE;
    reduction_map[NSAA::N] = POSITIVE;
    reduction_map[NSAA::Y] = POSITIVE;
    reduction_map[NSAA::H] = POSITIVE;
    reduction_map[NSAA::K] = POSITIVE;
    reduction_map[NSAA::R] = POSITIVE;

    reduction_map[NSAA::S] = NEGATIVE;
    reduction_map[NSAA::T] = NEGATIVE;
    reduction_map[NSAA::E] = NEGATIVE;
    reduction_map[NSAA::D] = NEGATIVE;

    reduction_map[NSAA::C] = CYS;

    is_init = true;
  }
  return reduction_map;
}


// this version somewhat places everything as either hydrophobic or
// polar (same categories as RASMOL)
//
const unsigned*
NSAA_MAPS::HP::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    reduction_map[NSAA::A] = H;
    reduction_map[NSAA::V] = H;
    reduction_map[NSAA::L] = H;
    reduction_map[NSAA::I] = H;
    reduction_map[NSAA::M] = H;
    reduction_map[NSAA::P] = H;

    reduction_map[NSAA::F] = H;
    reduction_map[NSAA::W] = H;

    reduction_map[NSAA::Y] = P;

    reduction_map[NSAA::G] = P;
    reduction_map[NSAA::S] = P;
    reduction_map[NSAA::T] = P;
    reduction_map[NSAA::C] = P;
    reduction_map[NSAA::N] = P;
    reduction_map[NSAA::Q] = P;

    reduction_map[NSAA::H] = P;
    reduction_map[NSAA::K] = P;
    reduction_map[NSAA::R] = P;

    reduction_map[NSAA::E] = P;
    reduction_map[NSAA::D] = P;

    is_init = true;
  }
  return reduction_map;
}


// Tyr arguably could go to the 'O'ther category here, I think moving
// it to 'P'olar, as in the the rasmol set, is going too far. Pro is
// left out because it has such a unique and extreme selection
// asymmetry in mammals
//
const unsigned*
NSAA_MAPS::HPO::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    reduction_map[NSAA::V] = H;
    reduction_map[NSAA::L] = H;
    reduction_map[NSAA::I] = H;
    reduction_map[NSAA::M] = H;
    reduction_map[NSAA::A] = H;

    reduction_map[NSAA::F] = H;
    reduction_map[NSAA::W] = H;
    reduction_map[NSAA::Y] = H;

    reduction_map[NSAA::N] = P;
    reduction_map[NSAA::S] = P;
    reduction_map[NSAA::T] = P;
    reduction_map[NSAA::Q] = P;

    reduction_map[NSAA::H] = P;
    reduction_map[NSAA::K] = P;
    reduction_map[NSAA::R] = P;

    reduction_map[NSAA::E] = P;
    reduction_map[NSAA::D] = P;

    reduction_map[NSAA::P] = O;
    reduction_map[NSAA::G] = O;
    reduction_map[NSAA::C] = O;

    is_init = true;
  }
  return reduction_map;
}



const unsigned*
NSAA_MAPS::BETA_BRANCH::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    for(unsigned i(0);i<NSAA::SIZE;++i) reduction_map[i] = OTHER;

    reduction_map[NSAA::V] = BETA_BRANCH;
    reduction_map[NSAA::I] = BETA_BRANCH;
    reduction_map[NSAA::T] = BETA_BRANCH;

    is_init = true;
  }
  return reduction_map;
}


const unsigned*
NSAA_MAPS::HELIX_FORM::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    for(unsigned i(0);i<NSAA::SIZE;++i) reduction_map[i] = OTHER;

    reduction_map[NSAA::A] = HELIX_FORM;
    reduction_map[NSAA::L] = HELIX_FORM;
    reduction_map[NSAA::K] = HELIX_FORM;
    reduction_map[NSAA::E] = HELIX_FORM;

    is_init = true;
  }
  return reduction_map;
}




const unsigned*
NSAA_MAPS::AROMATIC::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    for(unsigned i(0);i<NSAA::SIZE;++i) reduction_map[i] = OTHER;

    reduction_map[NSAA::F] = AROMATIC;
    reduction_map[NSAA::Y] = AROMATIC;
    reduction_map[NSAA::W] = AROMATIC;

    is_init = true;
  }
  return reduction_map;
}



const unsigned*
NSAA_MAPS::CYS::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    for(unsigned i(0);i<NSAA::SIZE;++i) reduction_map[i] = OTHER;

    reduction_map[NSAA::C] = CYS;

    is_init = true;
  }
  return reduction_map;
}



const unsigned*
NSAA_MAPS::GLY::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    for(unsigned i(0);i<NSAA::SIZE;++i) reduction_map[i] = OTHER;

    reduction_map[NSAA::G] = GLY;

    is_init = true;
  }
  return reduction_map;
}



const unsigned*
NSAA_MAPS::PRO::
get_reduction_map(){
  static bool is_init = false;
  static unsigned reduction_map[NSAA::SIZE];

  if(! is_init){
    for(unsigned i(0);i<NSAA::SIZE;++i) reduction_map[i] = OTHER;

    reduction_map[NSAA::P] = PRO;

    is_init = true;
  }
  return reduction_map;
}


const unsigned* NSAA_MAPS::converter::hpmap(HP::get_reduction_map());
