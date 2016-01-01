// -*- mode: c++; indent-tabs-mode: nil; -*-
// 


#include <bi_tree.h>
#include <newick_tree_parser.h>

#include <cxxtest/TestSuite.h>

#include <iostream>
#include <sstream>
#include <string>


class test_bi_tree : public CxxTest::TestSuite {
public:

  // test that newick_tree_parser accepts any well-formed tree:
  void testParse1() {
    TS_ASSERT_THROWS_NOTHING(newick_tree_parser("(A,(B,C,D,F));"));
  }
  void testParse2() {
    TS_ASSERT_THROWS_NOTHING(newick_tree_parser("(A:0.1,(B:0.1,C:1.3): 3.0);"));
  }
  void testParse3() {
    bi_tree("(,(,));");
    TS_ASSERT_THROWS_NOTHING(newick_tree_parser("(,,(,));"));
  }
  void testParse4() {
    TS_ASSERT_THROWS_NOTHING(newick_tree_parser("(,  (/t,/n /t)  : 21.5 )  rooooot   /n ;"));
  }
  void testParse5() {
    TS_ASSERT_THROWS_NOTHING(newick_tree_parser("(A big dog, ( in the woods : 0.0004 , eats a ) lot of) food;"));
  }


  // test that newick_tree_parser rejects several bad trees:
  void testParseReject1() {
    TS_ASSERT_THROWS_ANYTHING(newick_tree_parser("(A,(B,C))"));
  }
  void testParseReject2() {
    TS_ASSERT_THROWS_ANYTHING(newick_tree_parser("(A,B,C));"));
  }
  void testParseReject3() {
    TS_ASSERT_THROWS_ANYTHING(newick_tree_parser("(A,(B,C:LENGTH));"));
  }

  // test that bi_tee rejects additional (valid) trees, as expected:
  void testParseReject4() {
    TS_ASSERT_THROWS_ANYTHING(bi_tree("(A,(B,C,D));"));
  }
  void testParseReject5() {
    TS_ASSERT_THROWS_ANYTHING(bi_tree("(A C,(B,C));"));
  }


  void testWrite() {
    static const char tree[]="(A,(B,C)BC)ABC;";
    bi_tree t(tree);
    std::ostringstream oss;
    oss << t;
    TS_ASSERT_EQUALS(std::string(tree),oss.str());
  }


  void testOp1(){
    const bi_tree t("((A,B),C)R;");
    static const bi_tree_node* z(0);
    TS_ASSERT_EQUALS(t.node("R")->parent(),z);
  }
  void testOp2(){
    const bi_tree t("((A,B),C)R;");
    TS_ASSERT_EQUALS(t.root()->is_leaf(),false);
  }
  void testOp3(){
    const bi_tree t("((A,B),C)R;");
    TS_ASSERT_EQUALS(t.node("A")->is_leaf(),true);
  }
  void testOp4(){
    const bi_tree t("((A,B),C)R;");
    TS_ASSERT_EQUALS(t.node("A")->sister()->label(),std::string("B"));
  }

};

