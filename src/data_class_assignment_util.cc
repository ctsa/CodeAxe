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

// $Id: data_class_assignment_util.cc 1222 2008-05-22 23:10:06Z ctsa $

/// \file

#include "cat_manager.h"
#include "data_class_assignment_util.h"
#include "site_data.h"
#include "util/general/log.h"

#include <cstdlib>

#include <string>
#include <ostream>



data_class_id_assignment_map::
data_class_id_assignment_map(const cat_manager& cm,
                             const site_data& sd)
   : _assigned_data_set_size(cm.assigned_data_set_size()) {

  const unsigned n_data_classes(sd.data_class_size());

  // default state is one assigned data set:
  //
  _dci_ads_map.resize(n_data_classes,0);
  _is_dci_ads_map.resize(n_data_classes,true);

  const data_class_label_assignment_map_type& dc_ads_map(cm.data_class_label_assignment_map());
  if(! dc_ads_map.empty()){
    std::fill(_is_dci_ads_map.begin(),_is_dci_ads_map.end(),false);

    data_class_label_assignment_map_type::const_iterator i=dc_ads_map.begin(),i_end=dc_ads_map.end();
    for(;i!=i_end;++i){
      const std::string& data_class_label(i->first);
      const unsigned assigned_data_set_id(i->second);
      if( sd.data_class_labels.testid(data_class_label) ){
        const unsigned data_class_id(sd.data_class_labels.getid(data_class_label));
        _dci_ads_map[data_class_id]=assigned_data_set_id;
        _is_dci_ads_map[data_class_id]=true;
      } else {
        log_os << "ERROR:: input site data file does not contain data class: " << data_class_label << "\n";
        exit(EXIT_FAILURE);
      }
    }
  }
}
