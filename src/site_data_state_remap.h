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

// $Id: site_data_state_remap.h 743 2007-08-14 15:47:12Z ctsa $

/// \file

#ifndef __SITE_DATA_STATE_REMAP_H
#define __SITE_DATA_STATE_REMAP_H

struct site_data;

void
site_data_nsc5_to_nsc4pre(site_data& sd);

void
site_data_nsc5_to_nsc4post(site_data& sd);

void
site_data_nsc5_to_nscodon(site_data& sd);

void
site_data_nsc4_to_nscodon(site_data& sd);


#endif
