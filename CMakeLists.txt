project(network_tools)

add_executable(network_tools main.cpp)

include_directories(/home/mit/csb/fischer/libs/include)
link_directories(/home/mit/csb/fischer/libs/lib)
   
target_link_libraries(network_tools)
#target_link_libraries(network_tools iml cblas atlas gmp)