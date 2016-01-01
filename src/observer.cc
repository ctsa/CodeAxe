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

// $Id: observer.cc 1153 2008-03-18 00:24:07Z ctsa $

/// \file


#include "observer.h"



observer::
~observer() {
  nots_t::iterator i(_nots.begin()),i_end(_nots.end());
  for(;i!=i_end;++i) (*i)->unregister_observer(this);
}



void
observer::
observe_notifier(const notifier& n) {
  n.register_observer(this);
  _nots.insert(&n);
}



void
notifier::
notify_observers(const EVENT_TYPE::index_t e) const {
  obss_t::iterator i(_obss.begin()),i_end(_obss.end());
  for(;i!=i_end;++i) (*i)->recieve_notifier_event(this,e);
}
