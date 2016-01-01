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

// $Id: CodeAxe.cc 1222 2008-05-22 23:10:06Z ctsa $

/// \file

#include "bg_gtor_options.h"
#include "cat_info.h"
#include "cat_post_prob.h"
#include "lhood_model.h"
#include "rate_gtor_nscodon.h"
#include "rate_gtor_nuc.h"
#include "root_gtor.h"
#include "process_fsa.h"
#include "site_data_fastlup_util.h"
#include "site_data_report.h"
#include "site_data_state_remap.h"
#include "site_data_util.h"
#include "sim.h"
#include "subs_ml_confidence.h"
#include "substk_exception.h"
#include "subs_ml_print_util.h"
#include "subs_ml_model_ci_info.h"
#include "subs_ml_model_init_options.h"
#include "subs_ml_model_min_options.h"
#include "subs_ml_model_train.h"
#include "tang_ei.h"
#include "util/bio/bioseq_util.h"
#include "util/general/die.h"
#include "util/general/log.h"
#include "util/general/io_util.h"
#include "util/general/metatags.h"
#include "util/math/math_util_exception.h"
#include "util/math/random_util.h"

#include <boost/lexical_cast.hpp>

#include <cstdlib>
#include <cstring>

#include <fstream>
#include <ostream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

using boost::lexical_cast;


const unsigned DEFAULT_MIN_MATCH(5);
const char stdstream_tag[] = "-";

const string progname("CodeAxe");
string binname;
//const string prog_version("1.0.0");
const string prog_version("D92");

/// \todo -- find a way to get the parent dir revision number to go here automatically:
///
const string svn_revision("$Rev: 1222 $");


namespace MODE {
  enum index_t {
    NONE,
    PROCESS,
    TRAIN,
    LHOOD,
    SIM,
    CREATE_MODEL,
    NORM,
    REPORT_MODEL,
    REPORT_SEQ,
    REPORT_DATA,
    CONVERT_DATA,
    CONFIDENCE,
    SITE_CAT_POSTP,
    GROUP_CAT_POSTP,
    EI
  };
}



static
std::string
version_string(){
  static const unsigned offset(6);
  return progname+" version "+prog_version+"-r"+svn_revision.substr(offset,svn_revision.size()-offset-2);
  //+" /  ("+binname+")";
}



void
usage_error(const char* xmessage = 0)  NORETURN_TAG;

void
usage_error(const char* xmessage) {

  static const char emphasis_tag[] = "-----------";

  std::ostream& os(log_os);

  os << "\n"
     << binname << " :  phylogenetic analysis and simulation tools\n"
     << "\n" << version_string() << "\n"
     << "\n\n"
     << emphasis_tag << " MODES: " << emphasis_tag << "\n"
     << "\n"
     << "-process-seq - process sequence files into site data file\n"
     << "-report-seq  - report sequence file info\n"
     << "  -site-model             - (see below)\n"
     << "  -in-seq {dir|file}      - file or dir of files for alignment in fasta format (.fsa)\n"
     << "  [-out {file|-}]         - write site data to file (default stdout)\n"
     << "  [-max-align-count N]    - read only the first N alignments in directory\n"
     << "  [-4fold]                - if site-model is nuc, filter for 4fold synon nuc sites\n"
     << "  [-gc-max x]             - filter any sequence with %gc > x\n"
     << "  [-gc-min x]             - filter any sequence with %gc <= x\n"
     << "  [-no-window-filter]     - turn off window quality filter\n"
     << "  [-no-ambig-gap-filter]  - allow ambiguous bases to be treated as a single type ambiguous site\n"
     << "  [-no-filter]            - no sequence window, ambig-gap, or coding sequence filtration\n"
     << "  [-min-window-match N]   - N matches of 20 in quality window (default << " << DEFAULT_MIN_MATCH << " <<)\n"
     << "  [-no-codon-border]      - use all possible codons (default syncs site count for codon/c4/c5)\n"
     << "  [-in-class {dir|file}]  - data class assignments for positions from -in-seq (s/.fsa/.class/) \n"
     << "  [-class-labels file]    - data class labels\n"
     << "  [-no-adjacent-nuc-diff] - filter sites to have no adjacent differences\n"
     << "  [-only-single-nuc-diff] - filter sites to have only single nuc differences\n"
     << "\n"
     << "-create-model - generate a new model\n"
     << "  [-in-data {file|-}]   - read site data from file\n"
     << "  [-out {file|-}]       - write model to file (default stdout)\n"
     << "  [-random-param]       - randomize starting parameters in a biologically plausible range\n"
     << "  [-seed N]             - set random number seed\n"
     << "  (any model flags)\n"
     << "\n"
     << "-ml - search for max lhood model parameters\n"
     << "  -in-model {file|-}    - read model from file\n"
     << "  -in-data {file|-}     - read site data from file\n"
     << "  [-out {file|-}]       - write model to file (default: stdout)\n"
     << "  (any model scoring flags)\n"
     << "\n"
     << "-lhood - report log data prob/model likelihood, ln(P(D|M))\n"
     << "  -in-model {file|-}    - read model from file\n"
     << "  -in-data {file|-}     - read site data from file\n"
     << "  (any model scoring flags)\n"
     << "\n"
     << "-confidence - report (1.-alpha) confidence intervals on model parameters (alpha fixed@0.05)\n"
     << "  -in-model {file|-}  - read model from file (assumed mle of in-data)\n"
     << "  -in-data {file|-}   - read site data from file\n"
     << "  [-ci-search]        - solve ci's by an inversion of the likelihood ratio statistic,\n"
     << "                         assuming the profile likelihood function is normal wrt each parameter (default)\n"
     << "  [-ci-fisher]        - solve ci's by the inverse fisher information matrix\n"
     << "  [-ci-alpha x]       - 100*(1-x)% confidence range (default: x=" << DEFAULT_ALPHA << ")\n"
     << "  [-ci-mutation-only] - evaluate confidence intervals for mutation parameters only (search mode only)\n"
     << "  [-ci-selection-only]- evaluate confidence intervals for aa selection parameters only (search mode only)\n"
     << "  [-ci-exclude-root]  - exclude root parameters from ci calc (search mode only)\n"
     << "  [-ci-single]        - evaluate the ci for a single parameter only (search mode only)\n"
     << "  -out file           - write ci info to file, multiple processes given same filename will partition job\n"
     << "\n"
     << "-norm - normalize model parameters inc. set bg to stationary distro\n"
     << "  -in-model {file|-}  - read model from file\n"
     << "  [-out {file|-}]     - write normalized model to file (default stdout)\n"
     << "\n"
     << "-sim - simulate model data\n"
     << "  -in-model {file|-}         - read model from file\n"
     << "  -sim-size N                - size of simulated data\n"
     << "  [-sim-model arg]           - set simulator type, arg is one of:\n"
     << "    continuous = continuous time simulator (default)\n"
     << "    discrete   = discrete time simulator\n"
     << "    iss        = independent site sampler: alternate continuous time method; indy mutation models only\n"
     << "  [-sim-report-time]         - use non-iss simulator to estimate branch time in synon subs per nuc\n"
     << "  [-sim-group-size-min N]    - min of uniform distribution for group size (default: 100)\n"
     << "  [-sim-group-size-max N]    - max of uniform distribution for group size (default: 200)\n"
     << "  [-sim-assigned-cat-prob p] - prob of prior category knowledge at each site (default: 0.)\n"
     << "  [-sim-assigned-cat]        - category at each site is known (shortcut for -sim-assigned-cat-prob 1)\n"
     << "  [-site-model x]            - site model for simulated data (default: in-model site type)\n"
     << "  [-sim-sequence-output]     - output nucleotide sequence data instead of site data\n"
     << "  [-out {file|-}]            - write simulated site data to file (default stdout)\n"
     << "  [-seed N]                  - set random number seed\n"
     << "\n"
     << "-report-model - summarize model parameters\n"
     << "  -in-model {file|-}        - read model from file\n"
     << "  [-in-confidence {file|-}] - read free parameter confidence intervals from file\n"
     << "  [-out {file|-}]           - write report to file (default stdout)\n"
     << "\n"
     << "-report-data  - summarize model data\n"
     << "  -in-data {file|-}    - read site data from file\n"
     << "  [-pretty-print-data] - decode data and pretty print\n"
     << "  [-out {file|-}]      - write report to file (default stdout)\n"
     << "\n"
     << "-convert-data - convert data to new basic model type\n"
     << "  -site-model        - (see below)\n"
     << "  -in-data {file|-}  - read site data from file\n"
     << "  [-out {file|-}]    - write site data to file (default stdout)\n"
     << "\n"
     << "-site-cat-post-prob - get site category posterior prob for each site\n"
     << "  -in-model {file|-}           - read model from file\n"
     << "  -in-data {file|-}            - read site data from file\n"
     << "\n"
     << "-group-cat-post-prob - get group category posterior prob for each group\n"
     << "  -in-model {file|-}           - read model from file\n"
     << "  -in-data {file|-}            - read site data from file\n"
     << "\n"
     << "-ei - calculate (an approximation of) the evolutionary index matrix of Tang et al.\n"
     << "  -in-data {file|-}  - read site data from file, assumed to be only two species\n"
     << "  -org1\n"
     << "  -org2\n"
     << "\n\n";
  os << emphasis_tag << " MODEL FLAGS (used in create-model mode): " << emphasis_tag << "\n"
     << "\n"
     << "{-tree arg | -tree-file file}\n"
     << "  arg    = newick tree string\n"
     << "  file   = file containing newick tree string\n"
     << "    note:\n"
     << "      - tree must be bifurcating and rooted\n"
     << "      - leaf node names must match datafile seq names\n"
     << "      - interior node names are optional\n"
     << "      - branch times will be used for models with a single time category\n"
     << "    examples:\n"
     << "      -tree \"((zebra,donkey),unicorn);\"\n"
     << "      -tree \"((zebra:0.2,donkey:0.3):1.1,unicorn:1.4);\"\n"
     << "      -tree-file mammals.tree\n"
     << "\n"
     << "-reversible-tree - lock one time parameter so that reversible models\n"
     << "                   are not overdetermined\n"
     << "\n"
     << "-site-model arg  - set model site type, arg is one of:\n"
     << "  nuc     = nucleotide\n"
     << "  dinuc   = overlapping dinuc conditioned to single nuc\n"
     << "  trinuc  = overlapping trinuc conditioned to single nuc\n"
     << "  codon   = codon\n"
     << "  c4-pre  = 'extended codon' 4-mer conditioned on lst nuc to single codon\n"
     << "  c4-post = 'extended codon' 4-mer conditioned on 4th nuc to single codon\n"
     << "  c5      = 'extended codon' 5-mer conditioned on 1st & 5th nuc\n"
     << "  c5-pre  = 'extended codon' 5-mer conditioned on 1st & 2nd nuc\n"
     << "  c5-post = 'extneded codon' 5-mer conditioned on 4th & 5th nuc\n"
     << "  binary  = 2-state abstract model (for debugging)\n"
     << "\n"
     << "[-rate-model arg] - set mutation rate type, arg is one of:\n"
     << "  jc69    = jukes-cantor-1969\n"
     << "  k80     = kimura-1980\n"
     << "  f81     = felsenstein-1981\n"
     << "  hky85   = hasegawa-kishino-yano-1985\n"
     << "            - rev for nuc site only\n"
     << "  gtr     = general time reversible\n"
     << "            - rev for nuc site only\n"
     << "  gy94    = goldman-yang-1994\n"
     << "            - for codon site only\n"
     << "            - rev when a symmetric selection model is chosen\n"
     << "  nonrev  = unrestricted (default)\n"
     << "\n"
     << "[-context-model arg] - set mutation context, arg is one of:\n"
     << "  indy             = independent nucleotide parameters (default)\n"
     << "  pre-doublet      = f(N-1) context modifier\n"
     << "  post-doublet     = g(N+1) context modifier\n"
     << "  factored-triplet = f(N-1)*g(N+1) context modifier\n"
     << "  triplet          = h(N-1,N+1) context modifier\n"
     << "  cpg-only         = add two rates for CpG transitions and tranversions\n"
     << "  tpa-only         = add two rates for TpA transitions and tranversions\n"
     << "  cpg-1ti          = add one rate for CpG transition\n"
     << "  cpg-2ti          = add two rates; one for each CpG transition\n"
     << "  cpg-nonrev       = add 8 rates, one for every CpG context\n"
     << "\n"
     << "[-no-edge-correction]   - turn off context edge correction\n"
     << "[-free-edge-correction] - free edge-correction strength parameters (default: locked at 1)\n"
     << "\n"
     << "[-select-model arg] - set aa selection for coding models, arg is one of:\n"
     << "  none              = no selection\n"
     << "  single            = one aa selection parameter\n"
     << "  hp                = hydrophobic-polar selection\n"
     << "  from              = one selection parameter per source aa (20 fp)\n"
     << "  to                = one selection parameter per sink aa (20 fp)\n"
     << "  from-to           = product of source and sink aa parameters (40 fp)\n"
     << "  symmetric         = general aa symmetric selection (75 fp)\n"
     << "  from-symmetric    = product of source aa and symmetric exchange parameters (95 fp)\n"
     << "  symmetric-to      = product of sink aa and symmetric exchange parameters (95 fp)\n"
     << "  from-symmetric-to = product of source aa, sink aa, and symmetric exchange parameters (115 fp)\n"
     << "  asymmetric        = general aa selection (150 fp) (default)\n"
     << "\n"
     << "[-codon-bias-model arg] - set selection on codon usage, arg is one of:\n"
     << "  none             = (default)\n"
     << "  synon-ratio      = ratio for synonymous changes only (41 fp)\n"
     << "  synon-norm-ratio = ratio for all changes, but parameters are normed for each aa group (41 fp)\n"
     << "\n"
     << "[-c5-approx-modelar arg] - set approximation option for c5 model, arg is one of:\n"
     << "  none     = (default)\n"
     << "  c4-joint = approx c5 lhood by conditioning 2 c4 blocks\n"
     << "  c4-gmean = use geometric mean of 2 c4 blocks\n"
     << "\n"
     << "[-root-model arg] - set root state distribution, arg is one of:\n"
     << "  full         = one parameter per state, minus restrictions imposed by overlap (default)\n"
     << "  codon        = derive from an nscodon distro (c4/c5 site-models)\n"
     << "  codon-dinuc  = derive from an nscodon & codon edge dinuc distro (c4/c5 site-models)\n"
     << "  codon-trinuc =\n"
     << "  nuc          = derive from nuc distro (dinuc/trinuc/codon/c4 site-models)\n"
     << "  nuc-pos      = derive from codon position nuc distro (codon/c4 site-models)\n"
     << "  obs-avg      = average observed leaf state distribution\n"
     << "  obs-time-avg = branchtime averaged observed leaf state distribution\n"
     << "  least-sq     = per-branch least-squares error (BLS) root\n"
     << "  gc-shift     = full distro in one cat w/ gc moded derivatives in other root cats\n"
     << "\n"
     << "[-subs-rate-bg-model arg] - set rate matrix background distro type, arg is one of:\n"
     << "  obs-avg      = average observed leaf state distribution (default)\n"
     << "  full         = one parameter per state, minus restrictions imposed by overlap\n"
     << "  codon-dinuc  = XXX\n"
     << "\n"
     << "[{-cat-model arg | -cat-model-file file}] - specify model categories\n"
     << "  arg    = category expression\n"
     << "  file   = file containing category expression\n"
     << "\n"
     << "  category expression syntax:\n"
     << "    cat_expression:\n"
     << "      seq_cat_definitions | seq_cat_definitions data_set_cat_mapping\n"
     << "\n"
     << "    seq_cat_definitions:\n"
     << "      seq_cat_definition | seq_cat_definition seq_cat_definitions\n"
     << "    seq_cat_definition:\n"
     << "      seq_cat_label \"{\" param_set_definitions \"}\"\n"
     << "    param_set_definitions:\n"
     << "      param_set_definition | param_set_definition param_set_definitions\n"
     << "    param_set_definition:\n"
     << "      param_set_code \"{\" param_set_mappings \"}\"\n"
     << "    param_set_mappings:\n"
     << "      param_set_mapping_rule | param_set_mapping_rule \",\" param_set_mappings\n"
     << "    param_set_mapping_rule:\n"
     << "      param_set_label | param_set_label \":\" tree_node_label\n"
     << "\n"
     << "    data_set_cat_mapping:\n"
     << "      data_set_cat_map_rule | data_set_cat_map_rule data_set_cat_mapping\n"
     << "    data_set_cat_map_rule:\n"
     << "     '[' data_set_definition data_set_seq_cats ']'\n"
     << "    data_set_definition:\n"
     << "      assigned_data_set_label '{' data_class_labels '}'\n"
     << "    data_class_labels:\n"
     << "      data_class_label | data_class_label ',' data_class_labels\n"
     << "    data_set_seq_cats:\n"
     << "      seq_cat_label | seq_cat_label data_set_seq_cats\n"
     << "\n"
     << "  where:\n"
     << "    param_set_code is a 3 character code: \"XYY\", where:\n"
     << "                  X is the mixture_type_code; one of...\n"
     << "                    s = site mixture\n"
     << "                    g = group mixture\n"
     << "                  YY is the param_type_code; one of...\n"
     << "                    mr = mutation rate\n"
     << "                    mm = mutation model\n"
     << "                    ss = selection strength\n"
     << "                    sm = selection matrix\n"
     << "                    ro = root params (meaningless when using obs-x/least-sq)\n"
     << "                    ti = branch time\n"
     << "                    ob = observed state distro\n"
     << "                    bg = substitution rate background distro\n"
     << "                    mo = shortcut for complete model: mm+sm+ro+ti+ob\n"
     << "    param_set_label is a label used to re-reference the same parameters in any param_set_definition\n"
     << "                    with matching param_set_code, allowing a very general form of parameter tying\n"
     << "    tree_node_label is a label from the model's phylogenetic tree. The corresponding param_set_label \n"
     << "                    (and the parameters tied to it) is applied to the branch above this node and all\n"
     << "                    branches below it, until a node with its own param_set_label is found. The root\n"
     << "                    node is used when tree_node_label is unspecified.\n"
     << "\n"
     << " defaults: none for seq_cat_definitions, if -cat-model[-file] is not used,\n"
     << "             then a single seq_cat is applied to all data\n"
     << "           default data_set_cat_mapping merges all data classes into one set\n"
     << "             to which all seq_cats are assigned\n"
     << "\n"
     << "  examples: (w/o data_cat_mapping)\n"
     << "    \"helix {ssm{0}sro{0}} non_helix {ssm{1}sro{1}}\"\n"
     << "      (separate selection matrix and root for two model cats)\n\n"
     << "    \"\"\"\n"
     << "    0 { gmr{0} gss{0} }\n"
     << "    1 { gmr{1} gss{0} }\n"
     << "    2 { gmr{2} gss{0} }\n"
     << "    3 { gmr{0} gss{1} }\n"
     << "    4 { gmr{1} gss{1} }\n"
     << "    5 { gmr{2} gss{1} }\n"
     << "    \"\"\"\n"
     << "      (full cross of 3 group mutation rate cats with 2 group selection strength cats)\n\n"
     << "    \"0{gmo{0}} 1{gmo{1}}\"\n"
     << "      (run a mixture of two full models on the the data\n\n"
     << "    \"0 { gss{ 0:mm, 1:rn, 2:mm_rn, 3:hs } }\"\n"
     << "      (run a separate selection strength param on each branch of the tree: \"((mm,rn)mm_rn,hs);\" )\n\n"
     << "    \"0 { gss{ nonrat_select, rat_select:rn } }\"\n"
     << "      (run a separate selection strength param on the branch terminating at rn)\n\n"
     << "    \"sel_each { gss{ 0:mm, 1:rn, 2:mm_rn, 3:hs } } sel_just_rat { gss{ nonrat_select, rat_select:rn } }\"\n"
     << "      (a group mixture model of the two trees described above)\n\n"
     << "\n"
     << "[-unlock-cat-prob] - train the category prior probs (default is fixed uniform distribution)\n"
#if 0
     << "    helix:H+G+I,non_helix:B+E+T+S+X\n"
     << "      (this makes a helix and non-helix model category,if data categories are DSSP codes)\n"
#endif
     << "\n\n";
  os << emphasis_tag << " PERSISTENT MODEL SCORING FLAGS (allowed in create/ml/lhood modes, saved in model file, used in ml/lhood modes): " << emphasis_tag << "\n"
     << "\n"
     << "[-tol]                - absolute convergence tolerance in lnP units (default: " << DEFAULT_CONVERGE_TOLERANCE << ")\n"
     << "[-prefer-diag-expm]   - prefer matrix diagonalization to series for exp(M)\n"
     << "[-use-category-em]    - use EM to find the partition of the analyzed data in categories,\n"
     << "                        in addition to regular parameters. cannot change model file.\n"
     << "[-equil-obs-update]   - \n"
     << "\n"
     << "[-min arg] - minimizer type, arg is one of:\n"
     << "  conj-dir      = conjugate direction set (powell 1964)\n"
     << "  conj-grad     = conjugate gradient (polak-ribiere)\n"
     << "  conj-hybrid   = " << DEFAULT_CONJ_GRAD_START_ITER << " iterations of conj-grad, followed by conj-dir (default)\n"
#ifdef SHOW_EXPERIMENTAL_OPTIONS
     << "  praxis        = praxis conjugate direction set method (brent 1973)\n"
     << "  praxis-hybrid = " << DEFAULT_CONJ_GRAD_START_ITER << " iterations of conj-grad, followed by praxis\n"
    //     << "  em            = use EM (using conj-dir for m-step min)\n"
#endif
     << "\n\n";
  os << emphasis_tag << " NONPERSISTENT MODEL SCORING FLAGS (used in ml/lhood modes, not saved in model file): " << emphasis_tag << "\n"
     << "\n"
     << "[-refine]             - start minimization at full tol\n"
//       << "[-start-tol]          - scale tol from start-tol to tol during minimization (default: 10)\n"
     << "[-rootcycle-mode]     -\n"
     << "[-max-steps n]        - stop minimization and report result after n iterations\n"
     << "[-pingpong file]      - periodically backup model state to file.{ping,pong}\n"
     << "\n";

  if(xmessage != 0) {
    os << "\n"
       << "******** COMMAND-LINE ERROR:: " << xmessage << " ********\n"
       << "\n";
  }

  exit(EXIT_FAILURE);
}




/// not for large files...
///
static
void
get_file_buffer(const std::string& filename,
                std::string& buffer){
  if(filename.empty()) {
    pass_away("get_file_buffer(): empty filename");
  }
  ifstream fis(filename.c_str());
  check_nonempty_istream(fis,filename.c_str());

  ostringstream oss;
  oss << fis.rdbuf();
  buffer=oss.str();
}



static
void
codeaxe_sdf_init(site_data_fastlup& sdf,
                 subs_ml_model& mdl,
                 std::istream& is){

  const cat_manager& cm(mdl.get_cat_manager());
  const bool is_read_group_info(cm.group_cat_size()>1);
  site_data sd;

  sd.load_state(is,is_read_group_info);
  const SITE_MODEL::index_t msm(mdl.get_rate_gtor().site_model());
  if(msm != sd.sm){
    log_os << "mdlsm/datasm: " << msm << " " << sd.sm << "\n";
    pass_away("Data file and model file site types conflict!");
  }

  // set untrained parameters from data, such as the background distro:
  mdl.init_untrained_data_dependencies(sd);

  mdl_site_data_fastlup_init(sd,mdl,sdf);

  // if we're reading in a model with stored group posterior probabilities, use
  // these to update the obs distribution:
  //
  if(mdl.opt().is_cat_em){
    mdl.check_cat_post_prob_data_match(sdf);
    mdl.update_model_obs_distro(sdf);
  }
}



static
bool
is_eq(const char* a,const char* b){ return strcmp(a,b)==0; }



void unknown_arg_error(const char*,const char*) NORETURN_TAG;

void
unknown_arg_error(const char* label, const char* arg){
  usage_error((string("unknown ")+label+" arg: "+arg).c_str());
}



static
MODE::index_t
mode_arg(const char* arg){
  using namespace MODE;

  if     (is_eq(arg,"-process-seq"))         return PROCESS;
  else if(is_eq(arg,"-ml"))                  return TRAIN;
  else if(is_eq(arg,"-create-model"))        return CREATE_MODEL;
  else if(is_eq(arg,"-norm"))                return NORM;
  else if(is_eq(arg,"-lhood"))               return LHOOD;
  else if(is_eq(arg,"-sim"))                 return SIM;
  else if(is_eq(arg,"-report-model"))        return REPORT_MODEL;
  else if(is_eq(arg,"-report-data"))         return REPORT_DATA;
  else if(is_eq(arg,"-report-seq"))          return REPORT_SEQ;
  else if(is_eq(arg,"-convert-data"))        return CONVERT_DATA;
  else if(is_eq(arg,"-confidence"))          return CONFIDENCE;
  else if(is_eq(arg,"-site-cat-post-prob"))  return SITE_CAT_POSTP;
  else if(is_eq(arg,"-group-cat-post-prob")) return GROUP_CAT_POSTP;
  else if(is_eq(arg,"-ei"))                  return EI;
  else                                       return NONE;
}



static
RATE_GTOR_MODEL::index_t
site_model_arg(const char* arg){
  using namespace RATE_GTOR_MODEL;

  if     (is_eq(arg,"nuc"))     return RATE_GTOR_MODEL::NUC;
  else if(is_eq(arg,"dinuc"))   return RATE_GTOR_MODEL::DINUC;
  else if(is_eq(arg,"trinuc"))  return RATE_GTOR_MODEL::TRINUC;
  else if(is_eq(arg,"codon"))   return RATE_GTOR_MODEL::CODON;
  else if(is_eq(arg,"c4-pre"))  return C4PRE;
  else if(is_eq(arg,"c4-post")) return C4POST;
  else if(is_eq(arg,"c5"))      return C5;
  else if(is_eq(arg,"c5-pre"))  return C5PRE;
  else if(is_eq(arg,"c5-post")) return C5POST;
  else if(is_eq(arg,"binary"))  return BINARY;
  else
    unknown_arg_error("-site_model",arg);
}



static
MIN::index_t
min_arg(const char* arg){
  using namespace MIN;

  if     (is_eq(arg,"conj-dir"))      return CONJ_DIR;
  else if(is_eq(arg,"conj-grad"))     return CONJ_GRAD;
  else if(is_eq(arg,"conj-hybrid"))   return CG_CD_HYBRID;
  else if(is_eq(arg,"praxis"))        return PRAXIS;
  else if(is_eq(arg,"praxis-hybrid")) return CG_PRAXIS_HYBRID;
  else if(is_eq(arg,"em"))            return EM;
  else
    unknown_arg_error("-min",arg);
}



static
RATE_MODEL_NUC::index_t
rate_model_arg(const char* arg){
  using namespace RATE_MODEL_NUC;

  if     (is_eq(arg,"jc69"))   return JC69;
  else if(is_eq(arg,"k80"))    return K80;
  else if(is_eq(arg,"f81"))    return F81;
  else if(is_eq(arg,"hky85"))  return HKY85;
  else if(is_eq(arg,"gtr"))    return REV;
  else if(is_eq(arg,"gy94"))   return GY94;
  else if(is_eq(arg,"nonrev")) return NONREV;
  else
    unknown_arg_error("-rate-model",arg);
}



static
CONTEXT_MODEL_NUC::index_t
context_model_arg(const char* arg){
  using namespace CONTEXT_MODEL_NUC;

  if     (is_eq(arg,"indy"))             return INDY;
  else if(is_eq(arg,"pre-doublet"))      return PRE_DOUBLET;
  else if(is_eq(arg,"post-doublet"))     return POST_DOUBLET;
  else if(is_eq(arg,"factored-triplet")) return FACTORED_TRIPLET;
  else if(is_eq(arg,"triplet"))          return TRIPLET;
  else if(is_eq(arg,"cpg-only"))         return CPG_ONLY;
  else if(is_eq(arg,"tpa-only"))         return TPA_ONLY;
  else if(is_eq(arg,"cpg-1ti"))          return CPG_1TI;
  else if(is_eq(arg,"cpg-2ti"))          return CPG_2TI;
  else if(is_eq(arg,"cpg-nonrev"))       return CPG_NONREV;
  else
    unknown_arg_error("-context-model",arg);
}



static
SELECT_MODEL::index_t
select_model_arg(const char* arg){
  using namespace SELECT_MODEL;

  if     (is_eq(arg,"none"))      return NONE;
  else if(is_eq(arg,"single"))    return SINGLE;
  else if(is_eq(arg,"hp"))        return HP;
  else if(is_eq(arg,"from"))      return FROM;
  else if(is_eq(arg,"to"))        return TO;
  else if(is_eq(arg,"from-to"))   return FROM_TO;
  else if(is_eq(arg,"symmetric")) return SYMM;
  else if(is_eq(arg,"from-symmetric")) {
    /// \todo -- check from-symmetric for a normalization trap!!
    usage_error("from-symmetric selection suspected of minimization trap due to incorrect normalization.");
    return FROM_SYMM;
  }
  else if(is_eq(arg,"symmetric-to"))      return SYMM_TO;
  else if(is_eq(arg,"from-symmetric-to")) return FROM_SYMM_TO;
  else if(is_eq(arg,"asymmetric"))        return ASYMM;
  else
    unknown_arg_error("-select-model",arg);
}



static
ROOT_GTOR_MODEL::index_t
root_model_arg(const char* arg){

  if     (is_eq(arg,"full"))         return ROOT_GTOR_MODEL::FULL;
  else if(is_eq(arg,"codon"))        return ROOT_GTOR_MODEL::CODON;
  else if(is_eq(arg,"codon-dinuc"))  return ROOT_GTOR_MODEL::CODON_DINUC;
  else if(is_eq(arg,"codon-trinuc")) return ROOT_GTOR_MODEL::CODON_TRINUC;
  else if(is_eq(arg,"nuc"))          return ROOT_GTOR_MODEL::NUC;
  else if(is_eq(arg,"nuc-pos"))      return ROOT_GTOR_MODEL::NUC_POS;
  else if(is_eq(arg,"obs-avg"))      return ROOT_GTOR_MODEL::OBS_AVG;
  else if(is_eq(arg,"obs-time-avg")) return ROOT_GTOR_MODEL::OBS_TIME_AVG;
  else if(is_eq(arg,"least-sq"))     return ROOT_GTOR_MODEL::LSPROB;
  else if(is_eq(arg,"gc-shift"))     return ROOT_GTOR_MODEL::FULL_GC_SHIFT;
  else
    unknown_arg_error("-root-model",arg);
}



static
SUBS_RATE_BG_MODEL::index_t
subs_rate_bg_model_arg(const char* arg){

  if     (is_eq(arg,"full"))         return SUBS_RATE_BG_MODEL::FULL;
  else if(is_eq(arg,"codon-dinuc"))  return SUBS_RATE_BG_MODEL::CODON_DINUC;
  else if(is_eq(arg,"obs-avg"))      return SUBS_RATE_BG_MODEL::OBS_AVG;
  else
    unknown_arg_error("-subs-rate-bg-model",arg);
}



struct mode_data {

  mode_data(int init_argc,char* init_argv[])
    : argc(init_argc), argv(init_argv), argmark(argc,false), argstr(argc),
      model_isp(0), data_isp(0),osp(&std::cout) {

    for(int i(0);i<argc;++i) argstr[i]=argv[i];
    if(argc>0) argmark[0]=true;

    ai.bininfo=version_string();
    for(int i(0);i<argc;++i){
      if(i) ai.cmdline += ' ';
      ai.cmdline += argv[i];
    }
  }

  void
  finalize_args(){
    for(int i(0);i<argc;++i){
      if(! argmark[i]) usage_error((string("Invalid argument: ")+argv[i]).c_str());
    }
  }

  void
  need_model_data(const bool is_m,const bool is_d){
    if(is_m) { if(model_isp==0) usage_error("mode requires model file"); }
    else     { if(model_isp!=0) usage_error("mode does not use model file"); }

    if(is_d) { if(data_isp==0) usage_error("mode requires data file"); }
    else     { if(data_isp!=0) usage_error("mode does not use data file"); }
  }

  std::istream& model_is() { return *model_isp; }
  std::istream& data_is()  { return *data_isp; }
  std::ostream& os()  { return *osp; }

  int argc;
  char** argv;
  std::vector<bool> argmark;
  std::vector<std::string> argstr;
  audit_info ai;

  std::istream* model_isp;
  std::istream* data_isp;
  std::ostream* osp;
};



static
RATE_GTOR_MODEL::index_t
get_site_model_arg(mode_data& md){

  RATE_GTOR_MODEL::index_t rgm(RATE_GTOR_MODEL::NONE);

  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;

    if       (md.argstr[i]=="-site-model"){
      if(rgm != RATE_GTOR_MODEL::NONE) usage_error("Multiple site-model args");
      if(++i>=md.argc) { usage_error("No site-model arg"); }
      rgm=site_model_arg(md.argv[i]);
    } else { continue; }

    md.argmark[i-1] = true;
    md.argmark[i] = true;
  }

  return rgm;
}



static
void
codeaxe_norm(mode_data& md){
  md.finalize_args();
  md.need_model_data(true,false);

  subs_ml_model mdl(md.model_is(),md.ai);

  mdl.norm();

  mdl.store_state(md.os());
}



static
void
codeaxe_postp(mode_data& md,
              const MODE::index_t mode){
  md.finalize_args();
  md.need_model_data(true,true);

  subs_ml_model mdl(md.model_is(),md.ai);

  site_data_fastlup sdf;
  codeaxe_sdf_init(sdf,mdl,md.data_is());

  //  train_subs_ml_model_obs_partition(sdf,mdl);

  const CPPM::index_t pm(mode == MODE::SITE_CAT_POSTP ? CPPM::SITE : CPPM::GROUP);
  smlfloat lnp,lnp_norm;
  simple_init_matrix<prob_t> ppi;
  get_lnprob_from_param(mdl,sdf,lnp,lnp_norm,pm,&ppi);

  make_cat_post_prob(mdl.get_cat_manager().cat_size(),sdf,pm,ppi.ptr());

  cat_post_prob_report(mdl.rate_gtor_model(),mdl.tree(),
                       mdl.get_cat_manager(),sdf,ppi.ptr(),pm,md.os());
}



static
void
codeaxe_convert(mode_data& md){
  RATE_GTOR_MODEL::index_t rgm(get_site_model_arg(md));
  if(rgm == RATE_GTOR_MODEL::NONE) {
    usage_error("no destination site-model type selected in convert-data mode");
  }

  md.finalize_args();
  md.need_model_data(false,true);

  site_data sd(md.ai);
  sd.load_state(md.data_is());

  const SITE_MODEL::index_t r_in(sd.sm);
  const SITE_MODEL::index_t r_out(RATE_GTOR_MODEL::convert_to_site_model(rgm));

  if(r_in == r_out){
    // trivial_case -- do nothing
  } else if(r_in == SITE_MODEL::NSC5) {
    if(r_out == SITE_MODEL::NSC4PRE) {
      site_data_nsc5_to_nsc4pre(sd);
    } else if (r_out == SITE_MODEL::NSC4POST) {
      site_data_nsc5_to_nsc4post(sd);
    } else if (r_out == SITE_MODEL::NSCODON) {
      site_data_nsc5_to_nscodon(sd);
    } else {
      pass_away("Invalid data conversion");
    }
  } else if((r_in == SITE_MODEL::NSC4PRE || r_in == SITE_MODEL::NSC4POST)
            && r_out == SITE_MODEL::NSCODON){
    site_data_nsc4_to_nscodon(sd);
  } else {
    pass_away("Invalid data conversion");
  }

  sd.sm = r_out;
  sd.store_state(md.os());
}



static
void
codeaxe_report_model(mode_data& md){

  std::string ci_file;

  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;
    if(md.argstr[i]=="-in-confidence"){
      if(++i>=md.argc) usage_error();
      ci_file = md.argstr[i];
      md.argmark[i-1] = true;
      md.argmark[i] = true;
      break;
    }
  }

  md.finalize_args();
  md.need_model_data(true,false);

  subs_ml_model mdl(md.model_is(),md.ai);

  if(!ci_file.empty()){
    subs_ml_model_ci_info mdl_ci;

    std::ifstream infp(ci_file.c_str());
    mdl_ci.load_state(infp);
    const unsigned vsize(mdl_ci.val.size());
    simple_array<smlfloat> var(vsize);
    for(unsigned i(0);i<vsize;++i) var[i] = mdl_ci.val[i].variance;
    mdl.attach_ci(var.begin());
  }

  mdl.report(md.os());
}



static
void
codeaxe_report_data(mode_data& md){

  bool is_pretty_print_data(false);

  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;
    if(md.argstr[i]=="-pretty-print-data"){
      is_pretty_print_data = true;
      md.argmark[i] = true;
      break;
    }
  }

  md.finalize_args();
  md.need_model_data(false,true);

  site_data sd;
  sd.load_state(md.data_is());
  if(sd.sm != SITE_MODEL::NONE){
    site_data_report(sd,is_pretty_print_data,md.os());
  }
}



static
void
codeaxe_ei(mode_data& md){

  string org1,org2;

  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;
    const int starti(i);
    if       (md.argstr[i]=="-org1"){
      if(++i>=md.argc) usage_error("No org1 argument");
      org1 = md.argstr[i];
    } else if(md.argstr[i]=="-org2"){
      if(++i>=md.argc) usage_error("No org2 argument");
      org2 = md.argstr[i];
    } else { continue; }

    md.argmark[starti] = true;
    md.argmark[i] = true;
  }

  if(org1.empty() || org2.empty()) usage_error("org1 and org2 required in ei mode");

  md.finalize_args();
  md.need_model_data(false,true);

  site_data sd;
  sd.load_state(md.data_is());

  get_ei(sd,org1,org2);
}



static
void
codeaxe_sim(mode_data& md){

  const RATE_GTOR_MODEL::index_t rate_gtor_model(get_site_model_arg(md));
  sim_options sim_opt;

  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;
    int starti(i);

    if       (md.argstr[i]=="-sim-size"){
      if(++i>=md.argc) usage_error();
      const int tmpi(lexical_cast<int>(md.argv[i]));
      if(tmpi <= 0 ) usage_error("Invalid sim size");
      sim_opt.size = tmpi;

    } else if(md.argstr[i]=="-sim-model"){
      if(++i>=md.argc) { usage_error("No site-model arg"); }
      if       (md.argstr[i]=="continuous"){
        sim_opt.method = SIM_MODEL::CONTINUOUS;
      } else if(md.argstr[i]=="discrete"){
        sim_opt.method = SIM_MODEL::DISCRETE;
      } else if(md.argstr[i]=="iss"){
        sim_opt.method = SIM_MODEL::ISS;
      } else { usage_error((string("Unknown sim-model arg: ")+md.argstr[i]).c_str()); }

    } else if(md.argstr[i]=="-sim-group-size-min"){
      if(++i>=md.argc) usage_error();
      const int tmpi(lexical_cast<int>(md.argv[i]));
      if(tmpi < 0 ) usage_error("Invalid group-sim-size-min");
      sim_opt.group_size_min = tmpi;

    } else if(md.argstr[i]=="-sim-group-size-max"){
      if(++i>=md.argc) usage_error();
      const int tmpi(lexical_cast<int>(md.argv[i]));
      if(tmpi < 0 ) usage_error("Invalid group-sim-size-max");
      sim_opt.group_size_max = tmpi;

    } else if(md.argstr[i]=="-sim-assigned-cat-prob"){
      if(++i>=md.argc) usage_error();
      sim_opt.assigned_cat_prob = lexical_cast<prob_t>(md.argv[i]);
      if(sim_opt.assigned_cat_prob < 0.) usage_error("Invalid sim-assigned-cat-prob");

    } else if(md.argstr[i]=="-sim-assigned-cat"){
      sim_opt.assigned_cat_prob = 2.;

    } else if(md.argstr[i]=="-sim-report-time"){
      sim_opt.is_report_time = true;

    } else if(md.argstr[i]=="-sim-sequence-output"){
      sim_opt.is_nuc_seq_output  = true;

    } else { continue; }

    md.argmark[starti] = true;
    md.argmark[i] = true;
  }

  md.finalize_args();
  md.need_model_data(true,false);

  if(sim_opt.size <= 0) {
    usage_error("Invalid/Missing sim size");
  }

  if(sim_opt.group_size_min > sim_opt.group_size_max){
    usage_error("Invalid sim-group-size values");
  }

  subs_ml_model mdl(md.model_is(),md.ai);

  const RATE_GTOR_MODEL::index_t rgm(mdl.rate_gtor_model());

  if(rate_gtor_model != RATE_GTOR_MODEL::NONE) {
    sim_opt.output_site_model = rate_gtor_model;
  } else {
    sim_opt.output_site_model = rgm;
  }

  nuc_seq_data nsd(md.ai);
  simulate_data(sim_opt,mdl,nsd);

  if(sim_opt.is_nuc_seq_output){
    nsd.store_state(md.os());
  } else {
    site_data sd;
    nuc_seq_to_site_data(sd,nsd,sim_opt.output_site_model,false,false,1,-1);
    sd.store_state(md.os());
  }
}



static
void
codeaxe_handle_seq(mode_data& md,
                   const MODE::index_t mode){

  md.need_model_data(false,false);

  string seq_file_node,class_file_node,class_label_file;

  double gc_max(1.0);
  double gc_min(0.0);

  bool is_nuc_4fold_filter(false);
  bool is_codon_border(true);

  bool is_no_adjacent_nuc_diff(false);
  bool is_single_nuc_diff(false);

  seq_filter_options fopt;

  fopt.min_match = DEFAULT_MIN_MATCH;

  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;
    const int starti(i);

    if(md.argstr[i]=="-max-align-count"){
      if(++i>=md.argc) usage_error();
      fopt.max_align_count = lexical_cast<unsigned>(md.argv[i]);

    } else if(md.argstr[i]=="-4fold"){
      is_nuc_4fold_filter = true;

    } else if(md.argstr[i]=="-gc-max"){
      if(++i>=md.argc) usage_error();
      gc_max = lexical_cast<double>(md.argv[i]);

    } else if(md.argstr[i]=="-gc-min"){
      if(++i>=md.argc) usage_error();
      gc_min = lexical_cast<double>(md.argv[i]);

    } else if(md.argstr[i]=="-no-window-filter"){
      fopt.is_use_window_filter = false;

    } else if(md.argstr[i]=="-no-ambig-gap-filter"){
      fopt.is_use_ambi_filter = false;
      fopt.is_use_gap_filter = false;

    } else if(md.argstr[i]=="-no-filter"){
      fopt.is_skip_filters = true;

    } else if(md.argstr[i]=="-min-window-match"){
      if(++i>=md.argc) usage_error();
      fopt.min_match = lexical_cast<unsigned>(md.argv[i]);

    } else if(md.argstr[i]=="-no-codon-border"){
      is_codon_border = false;

    } else if(md.argstr[i]=="-in-seq"){
      if(++i>=md.argc) usage_error();
      seq_file_node = md.argv[i];

    } else if(md.argstr[i]=="-in-class"){
      if(++i>=md.argc) usage_error();
      class_file_node = md.argv[i];

    } else if(md.argstr[i]=="-no-adjacent-nuc-diff"){
      is_no_adjacent_nuc_diff = true;

    } else if(md.argstr[i]=="-only-single-nuc-diff"){
      is_single_nuc_diff = true;

    } else if(md.argstr[i]=="-class-labels"){
      if(++i>=md.argc) usage_error();
      class_label_file = md.argv[i];

    } else { continue; }

    md.argmark[starti] = true;
    md.argmark[i] = true;
  }

  RATE_GTOR_MODEL::index_t rate_gtor_model(RATE_GTOR_MODEL::NONE);
  if( mode == MODE::PROCESS) rate_gtor_model=get_site_model_arg(md);

  md.finalize_args();

  bool is_readseq_ready(! seq_file_node.empty());
  if( ! is_readseq_ready ) usage_error("no -in-seq argument");

  if( mode == MODE::PROCESS && rate_gtor_model == RATE_GTOR_MODEL::NONE ){
    usage_error("no site-model selected in process-seq mode");
  }

  nuc_seq_data in_seq;
  process_seq_data(in_seq,fopt,seq_file_node.c_str(),
                   class_file_node.c_str(),class_label_file.c_str());

  if(mode == MODE::REPORT_SEQ){
    nuc_seq_data_report(in_seq,md.os());
  } else {
    site_data sd(md.ai);
    // convert filtered nucleotide data into site-model specific data for further analysis
    nuc_seq_to_site_data(sd,in_seq,rate_gtor_model,is_codon_border,
                         is_nuc_4fold_filter,gc_max,gc_min,
                         is_no_adjacent_nuc_diff,is_single_nuc_diff);
    sd.store_state(md.os());
  }
}



static
void
codeaxe_confidence(mode_data& md,
                   const std::string& outfile){

  md.need_model_data(true,true);

  conf_options co;
  bool is_method_set(false);
  bool is_selection_only(false);
  bool is_mutation_only(false);
  bool is_exclude_root(false);

  // set any new model training options:
  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;
    const int starti(i);
    if       (md.argstr[i]=="-ci-alpha"){
      if(++i>=md.argc) usage_error();
      co.alpha = lexical_cast<smlfloat>(md.argv[i]);

    } else if(md.argstr[i]=="-ci-fisher"){
      if(is_method_set) usage_error("multiple ci methods defined");
      co.is_fisher_method=true;
      is_method_set=true;

    } else if(md.argstr[i]=="-ci-search"){
      if(is_method_set) usage_error("multiple ci methods defined");
      co.is_fisher_method=false;
      is_method_set=true;

    } else if(md.argstr[i]=="-ci-mutation-only"){
      is_mutation_only=true;

    } else if(md.argstr[i]=="-ci-selection-only"){
      is_selection_only=true;

    } else if(md.argstr[i]=="-ci-exclude-root"){
      is_exclude_root=true;

    } else if(md.argstr[i]=="-ci-single"){
      co.is_single_param=true;

    } else { continue; }

    md.argmark[starti] = true;
    md.argmark[i] = true;
  }

  md.finalize_args();

  if(outfile.empty()) usage_error("mode requires explicit filename for output");

  if((is_mutation_only+is_selection_only+is_exclude_root) > 1){
    usage_error("only one ci subset flag allowed");
  }

  if(co.is_fisher_method){
    if(co.is_single_param){
      usage_error("single param option is meaningless for fisher conf interval calculation");
    }
    if(is_mutation_only){
      usage_error("mutation only param option is meaningless for fisher conf interval calculation");
    }
    if(is_selection_only){
      usage_error("selection only param option is meaningless for fisher conf interval calculation");
    }
    if(is_exclude_root){
      usage_error("no ci-exclude-root");
    }
  }

  // read model:
  subs_ml_model mdl(md.model_is(),md.ai);

  // read site data, store in fastlup format required by lhood-calling functions:
  site_data_fastlup sdf;
  codeaxe_sdf_init(sdf,mdl,md.data_is());

  // setup mask on ci calculations:
  const unsigned ps(mdl.param_size());
  simple_array<bool> ci_mask(ps,false);

  /// \todo improve brittle hacks used to set ci subsets
  ///
  if(is_mutation_only){
    std::fill(ci_mask.begin(),ci_mask.end(),true);

    const rate_gtor_nuc_base* rgns(dynamic_cast<const rate_gtor_nuc_base*>(&mdl.get_rate_gtor()));
    if(rgns){
      const unsigned st(mdl.get_time_gtor().param_size()+rgns->param_start_mut());
      const unsigned s(rgns->param_size_mut());
      for(unsigned i(st);i<(st+s);++i){
        ci_mask[i] = false;
      }
    }
  } else if(is_selection_only){
    std::fill(ci_mask.begin(),ci_mask.end(),true);

    const rate_gtor_nscodon_base* rgns(dynamic_cast<const rate_gtor_nscodon_base*>(&mdl.get_rate_gtor()));
    if(rgns){
      const unsigned st(mdl.get_time_gtor().param_size()+rgns->param_start_aa());
      const unsigned s(rgns->param_size_aa());
      for(unsigned i(st);i<(st+s);++i){
        ci_mask[i] = false;
      }
    }
  } else if(is_exclude_root){
    subs_ml_model mdl_copy(mdl);
    mdl_copy.set_is_train_param_state(false);
    mdl_copy.get_root_gtor_nonconst().set_is_train_param_state(true);

    mdl_copy.is_train_param_state(ci_mask.begin());
  }

  /// \todo unkludgify...
  ///
  const rate_gtor_nuc_base* rgns(dynamic_cast<const rate_gtor_nuc_base*>(&mdl.get_rate_gtor()));
  if(rgns){
    const unsigned nmmc(mdl.get_cat_manager().typed_cat_size(CAT_PARAM_TYPE::MUT_MODEL));
    for(unsigned i(0);i<nmmc;++i){
      if(rgns->context_model_nuc(i) == CONTEXT_MODEL_NUC::FACTORED_TRIPLET){
        warning("cannot calculate valid confidence intervals for factored-triplet mutation parameters");
        break;
      }
    }
  }

  get_subs_ml_model_conf_interval(sdf,mdl,co,ci_mask,outfile);
}



static
void
codeaxe_lhood(mode_data& md,
              const MODE::index_t mode){

  md.need_model_data(true,true);

  // read model:
  bool is_refine(false);
  subs_ml_model mdl(md.model_is(),md.ai);

  // set any new model training options:
  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;
    const int starti(i);
    if(md.argstr[i]=="-tol"){
      if(++i>=md.argc) usage_error();
      mdl.opt().converge_value = lexical_cast<smlfloat>(md.argv[i]);

    } else if(md.argstr[i]=="-max-steps"){
      if(++i>=md.argc) usage_error();
      mdl.opt().max_steps = lexical_cast<unsigned>(md.argv[i]);

    } else if(md.argstr[i]=="-pingpong"){
      if(++i>=md.argc) usage_error();
      mdl.opt().is_pingpong = true;
      mdl.opt().pingpongfile = md.argv[i];

    } else if(md.argstr[i]=="-min"){
      if(++i>=md.argc) usage_error();
      mdl.opt().min = min_arg(md.argv[i]);

    } else if(md.argstr[i]=="-prefer-diag-expm"){
      mdl.opt().is_prefer_diag_expm = true;

    } else if(md.argstr[i]=="-use-category-em"){
      mdl.opt().is_cat_em = true;

    } else if(md.argstr[i]=="-equil-obs-update"){
      mdl.opt().obs_update_mode = OBS_UPDATE_MODE::EQUIL;

    } else if(md.argstr[i]=="-rootcycle-mode"){
      mdl.opt().is_rootcycle = true;

    } else if(md.argstr[i]=="-refine"){
      is_refine = true;

    } else { continue; }

    md.argmark[starti] = true;
    md.argmark[i] = true;
  }

  md.finalize_args();

  // read site data, store in fastlup format required by lhood-calling functions:
  site_data_fastlup sdf;
  codeaxe_sdf_init(sdf,mdl,md.data_is());

  if(mode == MODE::TRAIN){
    try {
      train_subs_ml_model(sdf,mdl,is_refine);
    } catch (substk_exception& e) {
      log_os << "FATAL:: SubsTK EXCEPTION: " << e.what() << "\n"
             << "...caught in " << __FILE__ << ":" << __LINE__ << "\n"
             << "...dumping model state:\n";
      mdl.store_state(log_os);
      exit(EXIT_FAILURE);
    } catch(math_util_exception& e) {
      log_os << "FATAL:: MATH UTIL EXCEPTION: " << e.what() << "\n"
             << "...caught in " << __FILE__ << ":" << __LINE__ << "\n"
             << "...dumping model state:\n";
      mdl.store_state(log_os);
      exit(EXIT_FAILURE);
    }

    mdl.store_state(md.os());
  } else if(mode == MODE::LHOOD){

    smlfloat lnp,lnp_norm;
    get_lnprob_from_param(mdl,sdf,lnp,lnp_norm);

    const unsigned site_count(sdf.total_count);
    lnp_line_report(lnp,lnp_norm,site_count,md.os());
  }
}



static
void
codeaxe_create(mode_data& md){

  RATE_GTOR_MODEL::index_t rate_gtor_model(get_site_model_arg(md));

  if(md.model_isp) usage_error("mode does not use model file");

  subs_ml_model_init_options sml_opt(rate_gtor_model,md.ai);

#ifdef OBS_SWITCH
  OBS_TYPE::index_t obs_model(OBS_TYPE::NORMAL);
#endif

  enum tree_input_type {
    TPT_NONE,
    TPT_STR,
    TPT_FILE
  };
  tree_input_type tpt(TPT_NONE);
  string treefile;

  enum cme_input_type {
    CME_NONE,
    CME_STR,
    CME_FILE
  };
  cme_input_type cme(CME_NONE);
  string cat_model_file;

  for(int i(0);i<md.argc;++i){
    if(md.argmark[i]) continue;
    const int starti(i);

    if       (md.argstr[i]=="-tree"){
      if(++i>=md.argc) { usage_error("no -tree arg"); }
      if(tpt != TPT_NONE) usage_error("multiple tree arguments");
      tpt = TPT_STR;
      sml_opt.tree_buffer = md.argv[i];

    } else if(md.argstr[i]=="-tree-file"){
      if(++i>=md.argc) { usage_error("no -tree-file arg"); }
      if(tpt != TPT_NONE) usage_error("multiple phylogenetic tree arguments");
      tpt = TPT_FILE;
      treefile = md.argv[i];

    } else if(md.argstr[i]=="-reversible-tree"){
      sml_opt.tgo.is_reversible=true;

    } else if(md.argstr[i]=="-rate-model"){
      if(++i>=md.argc) usage_error("no -rate-model arg");
      sml_opt.nopt.set_all_rate_model(rate_model_arg(md.argv[i]));

    } else if(md.argstr[i]=="-context-model"){
      if(++i>=md.argc) usage_error("no -context-model arg");
      sml_opt.nopt.set_all_context_model(context_model_arg(md.argv[i]));

    } else if(md.argstr[i]=="-no-edge-correction"){
      sml_opt.nopt.is_use_edge_correction=false;

    } else if(md.argstr[i]=="-free-edge-correction"){
      sml_opt.nopt.is_free_edge_correction=true;

    } else if(md.argstr[i]=="-select-model"){
      if(++i>=md.argc) usage_error("no -select-model arg");
      sml_opt.copt.set_all_select_model(select_model_arg(md.argv[i]));

    } else if(md.argstr[i]=="-codon-bias-model"){
      if(++i>=md.argc) usage_error();
      if       (md.argstr[i]=="none"){
        sml_opt.copt.codon_bias_model = CODON_BIAS_MODEL::NONE;
      } else if(md.argstr[i]=="synon-ratio"){
        sml_opt.copt.codon_bias_model = CODON_BIAS_MODEL::SYNON_RATIO;
      } else if(md.argstr[i]=="synon-norm-ratio"){
        sml_opt.copt.codon_bias_model = CODON_BIAS_MODEL::SYNON_NORM_RATIO;
      } else { usage_error("unknown -codon-bias-model arg"); }

    } else if(md.argstr[i]=="-c5-approx-model"){
      if(++i>=md.argc) usage_error();
      if       (md.argstr[i]=="none"){
        sml_opt.c5t = C5_APPROX::NONE;
      } else if(md.argstr[i]=="c4-joint"){
        sml_opt.c5t = C5_APPROX::C4_JOINT;
      } else if(md.argstr[i]=="c4-gmean"){
        sml_opt.c5t = C5_APPROX::C4_GMEAN;
      } else { usage_error("unknown -c5-approx arg"); }

    } else if(md.argstr[i]=="-root-model"){
      if(++i>=md.argc) usage_error();
      sml_opt.rogm = root_model_arg(md.argv[i]);

    } else if(md.argstr[i]=="-subs-rate-bg-model"){
      if(++i>=md.argc) usage_error();
      sml_opt.bopt.bg = subs_rate_bg_model_arg(md.argv[i]);

    } else if(md.argstr[i]=="-cat-model"){
      if(++i>=md.argc || (md.argstr[i].empty() || md.argstr[i][0]=='-')) usage_error("no -cat-model arg");
      if(cme != CME_NONE) usage_error("multiple cat-model args");
      cme = CME_STR;
      sml_opt.catman_opt.cat_model_str=md.argv[i];

    } else if(md.argstr[i]=="-cat-model-file"){
      if(++i>=md.argc || (md.argstr[i].empty() || md.argstr[i][0]=='-')) usage_error("no -cat-model-file arg");
      if(cme != CME_NONE) usage_error("multiple cat-model args");
      cme = CME_FILE;
      cat_model_file = md.argv[i];

    } else if(md.argstr[i]=="-unlock-cat-prob"){
      sml_opt.catman_opt.is_lock_cat_prob = false;

    } else if(md.argstr[i]=="-random-param"){
      sml_opt.pinit = PARAM_INIT_TYPE::RANDOM;

    } else { continue;}

    md.argmark[starti] = true;
    md.argmark[i] = true;
  }

  md.finalize_args();

  // check input flags:
  if( rate_gtor_model == RATE_GTOR_MODEL::NONE ){
    usage_error("no site-model selected in create-model mode");
  }

  if(sml_opt.nopt.is_free_edge_correction && (! sml_opt.nopt.is_use_edge_correction)){
    usage_error("incompatible edge correction options");
  }

  // take care of the tree:
  if      (tpt == TPT_NONE) {
    usage_error("no tree specified");
  } else if(tpt == TPT_FILE) {
    get_file_buffer(treefile,sml_opt.tree_buffer);
  }

  if(cme == CME_FILE){
    get_file_buffer(cat_model_file,sml_opt.catman_opt.cat_model_str);
  }

  // initialize model
  subs_ml_model mdl(sml_opt);

  /////////////////////////////////////////////////
  /// after ctor mdl mods:
  ///


  // read (optional) site data:
  if(md.data_isp){

    site_data sd;
    sd.load_state(md.data_is());
    if(mdl.get_rate_gtor().site_model() != sd.sm){
      log_os << "rgm: " << mdl.get_rate_gtor().site_model() << " sd: " << sd.sm << "\n";
      pass_away("Data file and model file site types conflict!");
    }


    // check data_catmap labels, and set priors to model_cats
    const data_class_label_assignment_map_type& dc(mdl.get_cat_manager().data_class_label_assignment_map());
    data_class_label_assignment_map_type::const_iterator i=dc.begin(),i_end=dc.end();
    for(;i!=i_end;++i){
      if( i->first != UNASSIGNED_CAT_LABEL &&
          (! sd.data_class_labels.testid(i->first)) ){
        pass_away("data class assignment map references unknown data class label");
      }
    }

    // set trained model parameters from data, such as the starting
    // root distribution, and data category priors:
    //
    mdl.init_trained_data_dependencies(sd);

    // set untrained parameters from data (data parameters), such as
    // the background distribution. These will be reset for a new data
    // file associated with any model
    //
    mdl.init_untrained_data_dependencies(sd);
  }

  mdl.store_state(md.os());
}



static
void
try_main(int argc,char* argv[]){

  //////////////////////////////////////
  // mode independent setup:
  //
  binname=argv[0];

  std::ios_base::sync_with_stdio(false);

  mode_data md(argc,argv);

  /////////////////////////////////////
  // setup other modeless flags:
  //
  {
    bool seed_set(false);
    for(int i(0);i<argc;++i){
      if(md.argmark[i]) continue;

      if       (md.argstr[i]=="-seed"){
        if(++i>=argc) { usage_error("no -seed arg"); }
        if(seed_set) { usage_error("Seed already set"); }
        md.ai.seed=lexical_cast<long>(argv[i]);
        seed_set=true;
      } else { continue; }

      md.argmark[i-1] = true;
      md.argmark[i] = true;
    }

    md.ai.seed=random_init(md.ai.seed);
  }


  ///////////////////////////////////////
  // parse primary i/o file info:
  //
  string outfile;
  string indatafile,inmodelfile;

  for(int i(0);i<argc;++i){
    if(md.argmark[i]) continue;
    string* str_ptr(0);
    if     (md.argstr[i]=="-out")     { str_ptr=&outfile; }
    else if(md.argstr[i]=="-in-data") { str_ptr=&indatafile; }
    else if(md.argstr[i]=="-in-model"){ str_ptr=&inmodelfile; }
    else { continue; }

    if(++i>=argc || (md.argstr[i].empty() || (md.argstr[i].size()>1 && md.argstr[i][0]=='-'))){
      usage_error(("no filename following "+md.argstr[i-1]).c_str());
    }
    if( ! str_ptr->empty() ) usage_error((md.argstr[i-1]+" assigned multiple times").c_str());
    str_ptr->operator=(md.argstr[i]);
    md.argmark[i-1] = true;
    md.argmark[i] = true;
  }

  if(! inmodelfile.empty() && inmodelfile == indatafile){
    usage_error("-in-data and -in-model args match");
  }


  //////////////////////////////////////////
  // determine program mode:
  //
  using namespace MODE;

  index_t mode(NONE);

  for(int i(0);i<argc;++i){
    if(md.argmark[i]) continue;
    const MODE::index_t new_mode(mode_arg(argv[i]));
    if(new_mode==NONE) continue;

    if(mode != NONE) usage_error("multiple modes selected");
    mode = new_mode;
    md.argmark[i]=true;
  }


  ///////////////////////////////////////////////
  // setup primary file streams:
  //

  // setup input-model stream:
  ifstream model_isf;
  if(! inmodelfile.empty() ){
    if( inmodelfile == stdstream_tag ){
      md.model_isp = &cin;
    } else {
      model_isf.open(inmodelfile.c_str());
      md.model_isp = &model_isf;
    }
    check_nonempty_istream(md.model_is(),inmodelfile.c_str());
  }

  // setup input-data stream
  ifstream data_isf;
  if(! indatafile.empty() ){
    if( indatafile == stdstream_tag ){
      md.data_isp = &cin;
    } else {
      data_isf.open(indatafile.c_str());
      md.data_isp = &data_isf;
    }
    check_nonempty_istream(md.data_is(),indatafile.c_str());
  }

  const bool is_no_stream_out_mode(mode == CONFIDENCE);

  // setup output stream
  ofstream main_osf;
  if( ! (outfile.empty() || outfile == stdstream_tag || is_no_stream_out_mode) ){
    main_osf.open(outfile.c_str());
    md.osp = &main_osf;
  }


  /////////////////////////////////////////
  // take care o business!!
  //
  switch(mode){
  case NORM:            return codeaxe_norm(md);
  case SITE_CAT_POSTP:
  case GROUP_CAT_POSTP: return codeaxe_postp(md,mode);
  case CONVERT_DATA:    return codeaxe_convert(md);
  case REPORT_MODEL:    return codeaxe_report_model(md);
  case REPORT_DATA:     return codeaxe_report_data(md);
  case EI:              return codeaxe_ei(md);
  case SIM:             return codeaxe_sim(md);
  case PROCESS:
  case REPORT_SEQ:      return codeaxe_handle_seq(md,mode);
  case CONFIDENCE:      return codeaxe_confidence(md,outfile);
  case TRAIN:
  case LHOOD:           return codeaxe_lhood(md,mode);
  case CREATE_MODEL:    return codeaxe_create(md);
  default:
    usage_error("invalid/no mode selected");
  }
}



int
main(int argc,char* argv[]){

  // last chance to catch exceptions...
  //
  try{
    try_main(argc,argv);

  } catch (const substk_exception& e) {
    log_os << "FATAL:: SubsTK EXCEPTION: " << e.what() << "\n";
    log_os << "...caught in main()\n";
    exit(EXIT_FAILURE);

  } catch(const std::exception& e) {
    log_os << "FATAL:: EXCEPTION: " << e.what() << "\n";
    log_os << "...caught in main()\n";
    exit(EXIT_FAILURE);

  } catch(...) {
    log_os << "FATAL:: UNKNOWN EXCEPTION\n";
    log_os << "...caught in main()\n";
    exit(EXIT_FAILURE);
  }
  return EXIT_SUCCESS;
}
