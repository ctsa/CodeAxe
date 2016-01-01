// -*- mode: c++; indent-tabs-mode: nil; -*-
// 

#include <newick_tree_parser.h>

#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>


int main() {
  std::string s(std::istream_iterator<char>(std::cin), std::istream_iterator<char>());
  std::cout << std::fixed << std::setprecision(3) << newick_tree_parser(s.c_str()) << "\n";
};

