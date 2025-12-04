cmake_minimum_required(VERSION 3.13.4)
project(learncpp LANGUAGES CXX)

find_package(deal.II 9.6.0 REQUIRED
  HINTS /usr/share/cmake/deal.II $ENV{DEAL_II_DIR}
)

deal_ii_initialize_cached_variables()

add_executable(main main.cc)

deal_ii_setup_target(main)

# include(CTest)
# enable_testing()