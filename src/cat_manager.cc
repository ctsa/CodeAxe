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

// $Id: cat_manager.cc 1222 2008-05-22 23:10:06Z ctsa $

/// \file

#include "bi_tree.h"
#include "cat_expression_parser.h"
#include "cat_info.h"
#include "cat_manager.h"
#include "color_graph.h"
#include "name_id_lup.h"
#include "prob_gtor_param_simple.h"
#include "util/general/die.h"
#include "util/general/io_util.h"
#include "util/general/log.h"
#include "util/math/random_util.h"

#include <functional>
#include <iostream>
#include <iomanip>
#include <sstream>



static
bool
is_same_group(const unsigned* const * a,
              const unsigned* const * b,
              const unsigned n_branches){

  const CAT_MIX_TYPE::index_t mt(CAT_MIX_TYPE::GROUP);
  for(unsigned pti(0);pti<CAT_PARAM_TYPE::SIZE;++pti){
    const CAT_PARAM_TYPE::index_t pt(static_cast<CAT_PARAM_TYPE::index_t>(pti));
    const unsigned cat_type(typed_cat_index(pt,mt));
    for(unsigned i(0);i<n_branches;++i){
      if(a[i][cat_type] != b[i][cat_type]) return false;
    }
  }
  return true;
}


static
CAT_PARAM_TYPE::index_t
parse_cattype_param(const std::string& s){
  using namespace CAT_PARAM_TYPE;

  const std::string s2(s.substr(1,2));

  for(unsigned i(0);i<SIZE;++i){
    if(s2 == short_syms[i]) return static_cast<index_t>(i);
  }
  if(s2 == "mo") return SIZE;
  else {
    log_os << "ERROR:: invalid param_type in cat-model: " << s2 << "\n";
    abort();
  }
}



static
CAT_MIX_TYPE::index_t
parse_cattype_mix(const std::string& s){
  if     (s[0] == 's') return CAT_MIX_TYPE::SITE;
  else if(s[0] == 'g') return CAT_MIX_TYPE::GROUP;
  else pass_away("invalid mix type in cat model");
}



namespace MDS_BOGUS { // here b/c the enum can't go into the following function
  enum cat_g_t { DEP, INDEP };
}

/// - find "sets", which are connected dep<->indep cat graphs
/// - associate each dep cat no with a set
///
static
void
make_dep_set(const cat_manager& cm,
             const CAT_PARAM_TYPE::index_t dep_pt,
             const CAT_MIX_TYPE::index_t dep_mt,
             const CAT_PARAM_TYPE::index_t indep_pt,
             const CAT_MIX_TYPE::index_t indep_mt,
             std::vector<unsigned>& dep_set_no){

  using namespace MDS_BOGUS;

  color_graph<std::pair<cat_g_t,unsigned> > g;

  const unsigned n_cats(cm.cat_size());
  const unsigned n_branches(cm.branch_size());

  unsigned dntc(0);
  if(dep_pt==CAT_PARAM_TYPE::TIME) dntc=cm.typed_cat_size(CAT_PARAM_TYPE::TIME);
  unsigned intc(0);
  if(indep_pt==CAT_PARAM_TYPE::TIME) intc=cm.typed_cat_size(CAT_PARAM_TYPE::TIME);

  for(unsigned i(0);i<n_cats;++i){
    for(unsigned b(0);b<n_branches;++b){
      const unsigned di((dntc*b)+cm.typed_cat_no_from_cat_no_y_branch_id(i,b,dep_pt,dep_mt));
      const unsigned ii((intc*b)+cm.typed_cat_no_from_cat_no_y_branch_id(i,b,indep_pt,indep_mt));

      g.add_edge(std::make_pair(DEP,di),std::make_pair(INDEP,ii));
    }
  }

  const unsigned n_dep_cats(cm.typed_cat_size(dep_pt,dep_mt));
  dep_set_no.resize(n_dep_cats);

  // color connected graphs by set_no
  unsigned set_no(0);
  for(unsigned i(0);i<n_dep_cats;++i){
    const bool is_colored(g.vertex_is_colored(std::make_pair(DEP,i)));
    if(! is_colored){
      g.color_connected(std::make_pair(DEP,i),set_no);
      ++set_no;
    }
    dep_set_no[i] = g.vertex_color(std::make_pair(DEP,i));
  }
}



/// sets up global cat no sorting such that group type_cats change
/// more slowly in global cat order than site type_cats
///
struct cat_sorter : public std::binary_function<bool, const unsigned* const *, const unsigned* const *> {

  cat_sorter(const unsigned b) : _n_branches(b) {}

  bool operator()(const unsigned* const * x, const unsigned* const * y) {
    static const unsigned assumed_mix_size(2);

    assert(CAT_MIX_TYPE::SIZE == assumed_mix_size); // a new mix type will break the sort assumptions

    static const CAT_MIX_TYPE::index_t mix_order[] = {CAT_MIX_TYPE::GROUP,CAT_MIX_TYPE::SITE};

    for(unsigned i(0);i<assumed_mix_size;++i){
      const CAT_MIX_TYPE::index_t mt(mix_order[i]);
      for(unsigned pti(0);pti<CAT_PARAM_TYPE::SIZE;++pti){
        const CAT_PARAM_TYPE::index_t pt(static_cast<CAT_PARAM_TYPE::index_t>(pti));
        const unsigned cat_type(typed_cat_index(pt,mt));
        for(unsigned b(0);b<_n_branches;++b){
          if(x[b][cat_type] != y[b][cat_type])
            return (x[b][cat_type] < y[b][cat_type]);
        }
      }
    }

    // cats are equal!!
    pass_away("At least two categories are equivalent");

    return false;
  }

private:
  unsigned _n_branches;
};



static
bool
is_same_branch_paramcat(const unsigned* bc1,
                        const unsigned* bc2){
  for(unsigned i(0);i<CAT_TYPE_SIZE;++i) if(bc1[i]!=bc2[i]) return false;
  return true;
}



void
cat_manager::
init_cat_info(){

  const unsigned n_branches(branch_size());

  std::sort(_data.typed_cat_no.ptr(),
            _data.typed_cat_no.ptr()+_data.typed_cat_no.dim1(),
            cat_sorter(n_branches));

  const unsigned n_cats(cat_size());

  // setup sizes
  //
  for(unsigned i(0);i<CAT_TYPE_SIZE;++i){
    for(unsigned c(0);c<n_cats;++c){
      for(unsigned b(0);b<n_branches;++b){
        _data.typed_cat_size[i]=std::max(_data.typed_cat_size[i],_data.typed_cat_no[c][b][i]);
      }
    }
    _data.typed_cat_size[i] += 1;
  }

  // setup global_cat to group_cat map
  //
  _data.cat_to_group.resize(n_cats);
  unsigned group_cat(0);
  for(unsigned i(0);i<n_cats;++i){
    if(i>0 && ! is_same_group(_data.typed_cat_no[i],_data.typed_cat_no[i-1],n_branches)){
      group_cat++;
    }
    _data.cat_to_group[i] = group_cat;
  }

  // check for invalid category combinations:
  //
  for(unsigned pti(0);pti<CAT_PARAM_TYPE::SIZE;++pti){
    if(CAT_PARAM_TYPE::is_mixable[pti]) continue;
    const CAT_PARAM_TYPE::index_t pt(static_cast<CAT_PARAM_TYPE::index_t>(pti));
    if(typed_cat_size(pt,CAT_MIX_TYPE::GROUP) > 1 && typed_cat_size(pt,CAT_MIX_TYPE::SITE) > 1){
      log_os << "ERROR:: incompatible cats: can't have site and group versions of: " << CAT_PARAM_TYPE::syms[pti] << "\n";
      exit(EXIT_FAILURE);
    }
  }

  // create branch_cat_sets: branches with equal parameterizations
  // (excepting time) at each global cat number
  //
  _data.branch_cat_set_first_branch_id.resize(n_cats);
  _data.branch_cat_set_no.init(n_cats,n_branches);
  for(unsigned c(0);c<n_cats;++c){
    std::vector<unsigned>& set_1branch(_data.branch_cat_set_first_branch_id[c]);
    unsigned* set_no(_data.branch_cat_set_no[c]);
    const unsigned* const * tcnb(_data.typed_cat_no[c]);

    for(unsigned b(0);b<n_branches;++b){
      bool is_in_set(false);
      for(unsigned s(0);s<set_1branch.size();++s){
        if(is_same_branch_paramcat(tcnb[b],tcnb[set_1branch[s]])){
          set_no[b]=s;
          is_in_set=true;
          break;
        }
      }
      if(!is_in_set){
        set_no[b]=set_1branch.size();
        set_1branch.push_back(b);
      }
    }

  }

  // find the number of mutation-model<->tree normalization sets:
  make_dep_set(*this,CAT_PARAM_TYPE::MUT_MODEL,CAT_MIX_TYPE::EITHER,
               CAT_PARAM_TYPE::TIME,CAT_MIX_TYPE::EITHER,
               _data.mut_model_set_by_mut_model_cat);

  for(unsigned mt(0);mt<CAT_MIX_TYPE::SIZE;++mt){
    const CAT_MIX_TYPE::index_t mti(static_cast<CAT_MIX_TYPE::index_t>(mt));
    // find the number of mut_rate<->mut_model normalization sets:
    make_dep_set(*this,CAT_PARAM_TYPE::MUT_RATE,mti,
                 CAT_PARAM_TYPE::MUT_MODEL,CAT_MIX_TYPE::EITHER,
                 _data.mut_rate_set_by_mut_rate_cat[mt]);

    // find the number of selection_strength<->selection_model normalization sets:
    make_dep_set(*this,CAT_PARAM_TYPE::SEL_STRENGTH,mti,
                 CAT_PARAM_TYPE::SEL_MATRIX,CAT_MIX_TYPE::EITHER,
                 _data.sel_strength_set_by_sel_strength_cat[mt]);
  }
}



static
void
cascade_typed_cat_no(const bi_tree_node* r,
                     const unsigned tcnp,
                     unsigned** tcn,
                     const bool* const * is_tcn,
                     const unsigned cid){

  const unsigned bi(r->branch_id());
  if(! is_tcn[bi][cid]) tcn[bi][cid]=tcnp;

  if(r->is_leaf()) return;
  cascade_typed_cat_no(r->child1(),tcn[bi][cid],tcn,is_tcn,cid);
  cascade_typed_cat_no(r->child2(),tcn[bi][cid],tcn,is_tcn,cid);
}



void
cat_manager::
process_seq_cat_definitions(const cat_expression_parser::seq_cat_definitions& scd){

  typedef cat_expression_parser cep_t;

  // get the global cat size as the number of category labels &
  // use this num to setup dependent data structures:
  //
  const unsigned n_scdefs(scd.size());
  const unsigned n_cats_init(std::max(unsigned(1),n_scdefs));
  const unsigned n_branches(branch_size());

  _data.typed_cat_no.init(n_cats_init,n_branches,CAT_TYPE_SIZE,0);
  simple_matrix3d<bool> is_typed_cat_no(n_scdefs,n_branches,CAT_TYPE_SIZE,false);

  simple_matrix<unsigned> typed_cat_no_root(n_scdefs,CAT_TYPE_SIZE,0);
  simple_matrix<bool> is_typed_cat_no_root(n_scdefs,CAT_TYPE_SIZE,false);
  simple_matrix<bool> is_typed_cat_no_any_branch(n_scdefs,CAT_TYPE_SIZE,false);

  const bi_tree& tree(_smls.tree());

  for(unsigned i(0);i<n_scdefs;++i){

    const std::string& seq_cat_label(scd[i].first);
    if(_data.cat_label_lup.testid(seq_cat_label)){
      pass_away("cat_manager.process_seq_cat_definitions(): cannot accept replicated seq cat labels in seq cat definitions");
    }
    const unsigned cat_no(_data.cat_label_lup.assignid(seq_cat_label));

    const cep_t::param_set_definitions& psd(scd[i].second);
    const unsigned n_psdefs(psd.size());

    for(unsigned j(0);j<n_psdefs;++j){

      const std::string& param_set_code(psd[j].first);
      const CAT_PARAM_TYPE::index_t pt(parse_cattype_param(param_set_code));
      const CAT_MIX_TYPE::index_t mt(parse_cattype_mix(param_set_code));

      const cep_t::param_set_mappings& psms(psd[j].second);
      for(unsigned k(0);k<psms.size();++k){
        if( ! is_branch_cat_param_type(pt) ){
          if (k>0 || ! psms[k].second.empty()){
            pass_away("cat_manager.process_seq_cat_definitions(): non-branch param type used with branch cats");
          }
        }

        const std::string& param_set_label(psms[k].first);

        if(pt == CAT_PARAM_TYPE::SIZE){ // full model categories
          for(unsigned pi(0);pi<CAT_PARAM_TYPE::SIZE;++pi){
            if(CAT_PARAM_TYPE::is_mixable[pi]) continue;

            const CAT_PARAM_TYPE::index_t pti(static_cast<CAT_PARAM_TYPE::index_t>(pi));
            add_param_set_mapping_rule(pti,mt,param_set_label,
                                       typed_cat_no_root[cat_no],is_typed_cat_no_root[cat_no]);
          }

        } else if(is_branch_cat_param_type(pt)){

          const std::string& node_label(psms[k].second);

          if(! node_label.empty()){
            if(! tree.is_node_label(node_label)){
              pass_away((std::string("unknown tree node label in cat-model: ")+node_label).c_str());
            }
          }

          // identify the tree node associated with branch_label (or
          // root for no label)
          //
          if(node_label.empty() || tree.node(node_label)->is_root()){
            add_param_set_mapping_rule(pt,mt,param_set_label,
                                       typed_cat_no_root[cat_no],
                                       is_typed_cat_no_root[cat_no]);

          } else {
            is_typed_cat_no_any_branch[cat_no][typed_cat_index(pt,mt)]=true;

            const unsigned param_set_branch_id(tree.node(node_label)->branch_id());
            add_param_set_mapping_rule(pt,mt,param_set_label,
                                       _data.typed_cat_no[cat_no][param_set_branch_id],
                                       is_typed_cat_no[cat_no][param_set_branch_id]);
          }

        } else {
          add_param_set_mapping_rule(pt,mt,param_set_label,
                                     typed_cat_no_root[cat_no],is_typed_cat_no_root[cat_no]);
        }
      }
    }
  }

  // make categories assignments cascade down the tree from their
  // specified nodes:
  {
    // check for valid branch specification:
    //
    const bi_tree_node* root(tree.root());
    if(root->is_leaf()) pass_away("No single node trees.");
    const unsigned rc1_id(root->child1()->branch_id());
    const unsigned rc2_id(root->child2()->branch_id());

    for(unsigned i(0);i<n_scdefs;++i){
      for(unsigned j(0);j<CAT_TYPE_SIZE;++j){
        if(! is_typed_cat_no_any_branch[i][j]) continue;
        if(! ( is_typed_cat_no_root[i][j] ||
               (is_typed_cat_no[i][rc1_id][j] && is_typed_cat_no[i][rc2_id][j]))){
          pass_away("ambigous branch category pattern");
        }
      }
    }


    // do the category cascade:
    //
    for(unsigned i(0);i<n_scdefs;++i){
      for(unsigned j(0);j<CAT_TYPE_SIZE;++j){
        cascade_typed_cat_no(root->child1(),typed_cat_no_root[i][j],
                             _data.typed_cat_no[i],is_typed_cat_no[i],j);
        cascade_typed_cat_no(root->child2(),typed_cat_no_root[i][j],
                             _data.typed_cat_no[i],is_typed_cat_no[i],j);
      }
    }
  }

  // add default labels:
  //
  name_id_lup& cl(_data.cat_label_lup);
  if(cl.size() == 0) cl.assignid("default_seq_cat");

  for(unsigned i(0);i<CAT_PARAM_TYPE::SIZE;++i){
    for(unsigned j(0);j<CAT_MIX_TYPE::SIZE;++j){
      const unsigned ti(typed_cat_index(static_cast<CAT_PARAM_TYPE::index_t>(i),
                                        static_cast<CAT_MIX_TYPE::index_t>(j)));
      name_id_lup& tcl(_data.typed_cat_label_lup[ti]);
      if(tcl.size()==0) tcl.assignid(std::string("default_")+CAT_MIX_TYPE::short_syms[j]+CAT_PARAM_TYPE::short_syms[i]+"_param_set");
    }
  }
}



/// assumed to be called from ctor only
///
void
cat_manager::
parse_cat_expression(const std::string& ce_str){

  typedef cat_expression_parser cep_t;

  const cep_t cep(ce_str.c_str());
  const cep_t::cat_expression& ce(cep.ce());

  process_seq_cat_definitions(ce.first);

  init_cat_info(); /// \todo this method needs a better name

#if 0
  // for now, any data-cat info is an error:
  //
  const std::auto_ptr<cep_t::assigned_cat_expression> ace_new_ptr(new cep_t::assigned_cat_expression);
  const cep_t::assigned_cat_expression* ace_ptr(0);
  if       (ce.size()==0) {
    ace_ptr=ace_new_ptr.get();
  } else if(ce.size()>1 || (ce[0].first.first != UNASSIGNED_CAT_LABEL)){
    pass_away("cat_manager: assigned data cats are not currently allowed in the cat-expression");
  } else {
    ace_ptr=&(ce[0].second);
  }
  const cep_t::assigned_cat_expression& ace(*ace_ptr);
#endif
}



cat_manager::
cat_manager(const cat_manager_options& cmo,
            const cat_manager_sml_share& smls) : base_t(), _smls(smls) {

  parse_cat_expression(cmo.cat_model_str);

  const unsigned n_cats(cat_size());

  // one param per state -- simplest cat param setup:
  //
  _ptor.reset(new prob_gtor_param_simple(n_cats));
  register_param(ptor());

  // setup param locks
  //
  set_is_train_param_state( ! (cmo.is_lock_cat_prob || n_cats==1));
}



const char* const section_id = "cat_manager";
const char* const data_class_label_assignment_map_label ="data_class_label_assignment_map";
const char* const end_label = "END";

const iliner il(section_id);



void
cat_manager::
store_state(std::ostream& os) const {

  {
    data_class_label_assignment_map_type::const_iterator i=data_class_label_assignment_map().begin(),i_end=data_class_label_assignment_map().end();
    for(;i!=i_end;++i){
      os << section_id << " " << data_class_label_assignment_map_label << " " << i->first << " " << i->second << "\n";
    }
  }

  const unsigned n_cats(cat_size());
  const unsigned n_branches(branch_size());
  os << section_id << " n_cats " << n_cats << "\n";

  for(unsigned i(0);i<n_cats;++i){
    for(unsigned b(0);b<n_branches;++b){
      for(unsigned mt(0);mt<CAT_MIX_TYPE::SIZE;++mt){
        os << section_id << " cat_" << i << "_branch_" << b << "_"
           << CAT_MIX_TYPE::syms[mt];
        for(unsigned pt(0);pt<CAT_PARAM_TYPE::SIZE;++pt){
          os << "_" << CAT_PARAM_TYPE::short_syms[pt];
        }
        for(unsigned pt(0);pt<CAT_PARAM_TYPE::SIZE;++pt){
          os << " " << typed_cat_no_from_cat_no_y_branch_id(i,b,
                                                            static_cast<CAT_PARAM_TYPE::index_t>(pt),
                                                            static_cast<CAT_MIX_TYPE::index_t>(mt));
        }
        os << "\n";
      }
    }
  }

  for(unsigned i(0);i<n_cats;++i){
    os << section_id << " cat_label_" << i << " " << _data.cat_label_lup.getstr(i) << "\n";
  }

  for(unsigned i(0);i<CAT_PARAM_TYPE::SIZE;++i){
    for(unsigned j(0);j<CAT_MIX_TYPE::SIZE;++j){
      const unsigned ti(typed_cat_index(static_cast<CAT_PARAM_TYPE::index_t>(i),
                                        static_cast<CAT_MIX_TYPE::index_t>(j)));
      const name_id_lup& tcl(_data.typed_cat_label_lup[ti]);
      for(unsigned k(0);k< tcl.size();++k){
        os << section_id << " typed_cat_label_"
           << CAT_MIX_TYPE::short_syms[j] << CAT_PARAM_TYPE::short_syms[i]
           << "_" << k << " " << tcl.getstr(k) << "\n";
      }
    }
  }

  std::ostringstream oss;
  oss << section_id << " cat_prob_gtor";
  ptor().store_state(oss.str().c_str(),os);

  os << section_id << " " << end_label << "\n";
}



cat_manager::
cat_manager(std::istream& is,
            const cat_manager_sml_share& smls) : base_t(), _smls(smls) {

  const unsigned n_branches(branch_size());

  unsigned n_cats;

  // clear out the end of a previous line:
  clear_whitespace_line(is);
  std::string buf;
  while(getline(is,buf)){
    std::istringstream sbuf(buf);

    std::string tmpstr;
    sbuf >> tmpstr;
    if(tmpstr != section_id){
      log_os << "ERROR:: -- Invalid model file format. --\n"
             << "ERROR:: expected_section: " << section_id << "\n"
             << "ERROR:: section: " << tmpstr << "\n"
             << "ERROR:: line buffer: " << buf << "\n";
      abort();
    }

    sbuf >> tmpstr;

    unsigned tmpu;
    if       (tmpstr == data_class_label_assignment_map_label){
      sbuf >> tmpstr >> tmpu;
      _data.data_class_label_assignment_map[tmpstr] = tmpu;
    } else if(tmpstr == "n_cats"){
      sbuf >> n_cats;
      _data.typed_cat_no.init(n_cats,n_branches,CAT_TYPE_SIZE);
      break;
    } else {
      log_os << "ERROR:: invalid model file format. section: " << section_id << " subsection: " << tmpstr << "\n";
      abort();
    }
  }

  for(unsigned i(0);i<n_cats;++i){
    for(unsigned b(0);b<n_branches;++b){
      for(unsigned mt(0);mt<CAT_MIX_TYPE::SIZE;++mt){
        il.advance(is);
        for(unsigned pt(0);pt<CAT_PARAM_TYPE::SIZE;++pt){
          is >> typed_cat_no_from_cat_no_y_branch_id(i,b,
                                                     static_cast<CAT_PARAM_TYPE::index_t>(pt),
                                                     static_cast<CAT_MIX_TYPE::index_t>(mt));
        }
      }
    }
  }

  init_cat_info();


  { // read labels for global cats and typed cats:
    std::string label;

    _data.cat_label_lup.clear();
    for(unsigned i(0);i<n_cats;++i){
      il.advance(is);
      is >> label;

      if(_data.cat_label_lup.testid(label)) die("duplicate cat label in modelfile");
      _data.cat_label_lup.assignid(label);
    }

    for(unsigned i(0);i<CAT_PARAM_TYPE::SIZE;++i){
      for(unsigned j(0);j<CAT_MIX_TYPE::SIZE;++j){
        const unsigned ti(typed_cat_index(static_cast<CAT_PARAM_TYPE::index_t>(i),
                                          static_cast<CAT_MIX_TYPE::index_t>(j)));

        name_id_lup& tcl(_data.typed_cat_label_lup[ti]);
        const unsigned tis(_data.typed_cat_size[ti]);

        for(unsigned k(0);k<tis;++k){
          il.advance(is);
          is >> label;
          if(tcl.testid(label)) die("duplicate typed cat label in modelfile\n");
          tcl.assignid(label);
        }
      }
    }
  }

  // read model parameters
  _ptor.reset(new prob_gtor_param_simple(n_cats));
  register_param(ptor());

  prob_gtor_param& pt(ptor());
  const unsigned ps(pt.param_size());
  simple_array<smlfloat> param(ps);
  simple_array<bool> itp(ps);

  for(unsigned i(0);i<ps;++i){
    il.advance(is);
    bool btmp;
    is >> btmp >> param[i];
    itp[i] = btmp;
  }

  pt.set_param_state(param.begin());
  pt.set_is_train_param_state(itp.begin());

  il.advance(is,end_label);
}



void
cat_manager::
branch_typed_cat_pdistro(prob_t* tcp,
                         const unsigned branch_id,
                         const CAT_PARAM_TYPE::index_t pt,
                         const CAT_MIX_TYPE::index_t mt) const {

  const unsigned n_typed_cats(typed_cat_size(pt,mt));
  for(unsigned tc(0);tc<n_typed_cats;++tc) tcp[tc] = 0.;

  const unsigned n_cats(cat_size());
  const prob_t* cp(ptor().pdistro());

  for(unsigned c(0);c<n_cats;++c){
    const unsigned typed_cat(typed_cat_no_from_cat_no_y_branch_id(c,branch_id,pt,mt));
    tcp[typed_cat] += cp[c];
  }

  pdistro_norm(tcp,tcp+n_typed_cats);
}



const std::string&
cat_manager::
typed_cat_label(const unsigned tcn,
                const CAT_PARAM_TYPE::index_t pt,
                const CAT_MIX_TYPE::index_t mt) const{

  CAT_MIX_TYPE::index_t mt2(mt);
  if(mt == CAT_MIX_TYPE::EITHER){
    if(CAT_PARAM_TYPE::is_mixable[pt])  die("invalid call to typed_cat_label");
    if(typed_cat_size(pt,CAT_MIX_TYPE::SITE) < typed_cat_size(pt,CAT_MIX_TYPE::GROUP)){
      mt2=CAT_MIX_TYPE::GROUP;
    } else {
      mt2=CAT_MIX_TYPE::SITE;
    }
  }
  const unsigned n_typed_cats(typed_cat_size(pt,mt2));
  if(tcn>=n_typed_cats) die("invalid typed_cat_no");
  return _data.typed_cat_label_lup[typed_cat_index(pt,mt2)].getstr(tcn);
}



const std::string&
cat_manager::
cat_label(const unsigned cn) const {
  return _data.cat_label_lup.getstr(cn);
}


void
cat_manager::
reset_param(const PARAM_INIT_TYPE::index_t pinit) {

  prob_gtor_param& pt(ptor());

  if(pinit == PARAM_INIT_TYPE::RANDOM){
    static const smlfloat rcprob_unorm_min(1e-3);

    const unsigned ps(param_size());
    simple_array<smlfloat> param(ps);
    for(unsigned i(0);i<(ps);++i){
      param[i] = random_uniform()+rcprob_unorm_min;
    }
    pt.set_param_state(param.ptr());

  } else if(pinit == PARAM_INIT_TYPE::START){
    const unsigned ss(pt.state_size());
    simple_array<smlfloat> distro(ss);
    pdistro_unif(distro.ptr(),distro.ptr()+ss);
    pt.set_pdistro(distro.ptr());

  } else { die("unknown param_init_type"); }
}


#if 0
bool
cat_manager::
is_missing_obs_partition() const {
  /// \todo fix this to work with the assigned cat case too
  ///
  if(assigned_cat_size()>0)
    die("obs partitions not yet setup for assigned cats");

  return (typed_cat_size(CAT_PARAM_TYPE::OBS)>1);
}
#endif



unsigned
cat_manager::
branch_size() const {
  return _smls.tree().branch_size();
}



const std::string&
cat_manager::
branch_label(const unsigned branch_id) const {
  return _smls.tree().branch_node(branch_id)->label();
}



void
cat_manager::
branch_cat_set_label(const unsigned cat_no,
                     const unsigned branch_cat_set,
                     std::string& label) const {

  static const char sep_char='+';

  const unsigned n_branches(branch_size());
  bool is_first(true);

  std::ostringstream oss;
  for(unsigned b(0);b<n_branches;++b){
    if(get_branch_cat_set(cat_no,b)==branch_cat_set){
      if(! is_first) oss << sep_char;
      oss << branch_label(b);
      is_first=false;
    }
  }
  label=oss.str();
}



void
cat_manager::
add_param_set_mapping_rule(const CAT_PARAM_TYPE::index_t& pt,
                           const CAT_MIX_TYPE::index_t& mt,
                           const std::string& label,
                           unsigned* typed_cat_no,
                           bool* is_typed_cat_no){

  const unsigned cindex(typed_cat_index(pt,mt));
  const unsigned tcn(_data.typed_cat_label_lup[cindex].assignid(label.c_str()));

  assert(cindex<CAT_TYPE_SIZE);
  if(is_typed_cat_no[cindex]) {
    std::ostringstream oss;
    oss << "conflicting cat type assignment: param_type,mix_type,label,typed_cat_id: "
        << CAT_PARAM_TYPE::syms[pt] << " " << CAT_MIX_TYPE::syms[mt] << " "
        << label << " " << tcn << "\n";
    throw substk_exception(oss.str().c_str());
  }
  typed_cat_no[cindex]=tcn;
  is_typed_cat_no[cindex]=true;
}



bool
cat_manager::
is_missing_obs_partition() const {
  /// \todo fix this to work with the assigned cat case too
  ///
  if(assigned_data_set_size()>1)
    die("obs partitions not yet setup for multiple assigned data sets");

  return (typed_cat_size(CAT_PARAM_TYPE::OBS)>1);
}



static
void
cm_reporter(const cat_manager& cm,
            std::ostream& os) {

  os << "category_report:\n\n";


  const unsigned n_cats(cm.cat_size());
  os << "sequential_cat_count: " << n_cats << "\n\n";

  std::vector<prob_t> cat_pdistro(n_cats);
  cm.cat_pdistro(cat_pdistro.begin());
  for(unsigned c(0);c<n_cats;++c){
    os << "seq_cat_prob " << c << " " << cat_pdistro[c] << "\n";
  }
  os << "\n";

  for(unsigned c(0);c<n_cats;++c){
    os << "seq_cat_branch_set_count " << c << " " << cm.branch_cat_set_size(c) << "\n";
  }
  os << "\n";

  const unsigned n_group_cats(cm.group_cat_size());
  os << "group_sequential_cat_count: " << n_group_cats << "\n\n";

  std::vector<prob_t> gcat_pdistro(n_group_cats);
  cm.group_cat_pdistro(gcat_pdistro.begin());
  for(unsigned gc(0);gc<n_group_cats;++gc){
    os << "group_seq_cat_prob " << gc << " " << gcat_pdistro[gc] << "\n";
  }
  os << "\n";

  for(unsigned i(0);i<CAT_PARAM_TYPE::SIZE;++i){
    const CAT_PARAM_TYPE::index_t pt(static_cast<CAT_PARAM_TYPE::index_t>(i));

    if(CAT_PARAM_TYPE::is_mixable[pt]){
      for(unsigned j(0);j<CAT_MIX_TYPE::SIZE;++j){
        const CAT_MIX_TYPE::index_t mt(static_cast<CAT_MIX_TYPE::index_t>(j));
        os << "typed_cat_count " << CAT_PARAM_TYPE::syms[pt]
           << " " << CAT_MIX_TYPE::syms[mt] << " " << cm.typed_cat_size(pt,mt) << "\n";
      }
    } else {
      os << "typed_cat_count " <<  CAT_PARAM_TYPE::syms[pt]
         << " nomix " << cm.typed_cat_size(pt) << "\n";
    }
  }

  os << "\n\n";
}



void
cat_manager::
report(std::ostream& os) const {

  cm_reporter(*this,os);
}
