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

// $Id: observer.h 1153 2008-03-18 00:24:07Z ctsa $

/// \file

#ifndef __OBSERVER_H
#define __OBSERVER_H


#include <set>



namespace EVENT_TYPE {
  enum index_t {
    PDISTRO_CHANGE
  };
}



struct notifier;



struct observer {
  friend struct notifier;

  typedef observer self_t;

  observer() {}

  observer(const self_t&)  {} // do not copy notifier set

  virtual ~observer();

private:
  self_t& operator=(const self_t&);

  virtual void recieve_notifier_event(const notifier*,
                                      const EVENT_TYPE::index_t) = 0;

protected:
  void
  observe_notifier(const notifier& n);

private:
  typedef std::set<const notifier*> nots_t;

  mutable nots_t _nots;
};



struct notifier {
  friend struct observer;

  typedef notifier self_t;

  notifier() {}

  notifier(const self_t&)  {} // do not copy observer set

  virtual ~notifier() {}

private:
  self_t& operator=(const self_t&);

  void
  register_observer(observer* n) const { _obss.insert(n); }

  void
  unregister_observer(observer* n) const {
    const obss_t::iterator i(_obss.find(n));
    if(i != _obss.end()) _obss.erase(i);
  }

protected:
  void
  notify_observers(const EVENT_TYPE::index_t e) const;

private:
  typedef std::set<observer*> obss_t;

  mutable obss_t _obss;
};



#endif
