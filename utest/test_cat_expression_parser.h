// -*- mode: c++; indent-tabs-mode: nil; -*-
// 


#include <cat_expression_parser.h>

#include <cxxtest/TestSuite.h>

#include <iostream>
#include <sstream>
#include <string>


class test_cat_expression_parser : public CxxTest::TestSuite {
public:

  // test that some correct cat-expressions are accepted:
  void testParse1() {
    TS_ASSERT_THROWS_NOTHING(cat_expression_parser("cat1{gmm{foo1}gti{bar}}"));
  }
  void testParse2() {
    TS_ASSERT_THROWS_NOTHING(cat_expression_parser("cat1 {gmm{foo1 } gti{bar}} cat2 {gmm{foo}gti{bar2}}"));
  }
  void testParse3() {
    TS_ASSERT_THROWS_NOTHING(cat_expression_parser("0{gmm{0}gss{0:mm,1:rn,2:mm_rn,3:hs}} 1{gmm{1}}"));
  }
  void testParse4() {
    TS_ASSERT_THROWS_NOTHING(cat_expression_parser("0 {gmm{0}gss{0:mm,1:rn,2:mm_rn,3:hs}} 2 {gss{4}} 1{gmm{1}} [non_helix{X,H} 0 2] [helix{E,L} 1]"));
  }

  // test that the parser rejects several bad expressions:
  void testParseReject1() {
    TS_ASSERT_THROWS_ANYTHING(cat_expression_parser("{gmm{foo1}gti{bar}}"));
  }

  void testParseReject2() {
    TS_ASSERT_THROWS_ANYTHING(cat_expression_parser("cat1{}"));
  }

  void testParseReject3() {
    TS_ASSERT_THROWS_ANYTHING(cat_expression_parser("[non helix{X,H}]0{gmm{0}}"));
  }
};

