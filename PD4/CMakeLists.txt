cmake_minimum_required(VERSION 3.6)
project("PD4")

set( CMAKE_INSTALL_PREFIX /Volumes/Text/Documents-working/3-Current-Projects/3\ PD/3.pr.projects/papers/paper1/results  )
message("CIP = ${CMAKE_INSTALL_PREFIX}" )

set(CMAKE_CXX_STANDARD 11)
add_library(pdDeps STATIC PD.cpp phylogeny.cpp TaxonLists.cpp taxon.cpp Utils.cpp )

add_executable(QTree QTree.cpp  )
target_link_libraries(QTree pdDeps)

add_executable(BalTree BalTree.cpp )
target_link_libraries(BalTree pdDeps)

add_executable(YuleTree YuleTree.cpp )
target_link_libraries(YuleTree pdDeps)

add_executable(BDTree BDTree.cpp )
target_link_libraries(BDTree pdDeps)
 
add_executable( YuleSim YuleSim.cpp )
target_link_libraries(YuleSim pdDeps)
 

add_executable( BDSim BDSim.cpp )
target_link_libraries(BDSim pdDeps)


add_executable( BBDTreeSim BBDTreeSim.cpp )
target_link_libraries(BBDTreeSim pdDeps)


add_executable( QtreeSim QtreeSim.cpp )
target_link_libraries(QtreeSim pdDeps)

 
add_executable(PD PDtree.cpp  )
target_link_libraries(PD pdDeps)

add_executable(pdlist pdlist.cpp )
target_link_libraries(pdlist pdDeps)
  

add_executable(PDsampled PDsampled.cpp )
target_link_libraries(PDsampled pdDeps)
 

add_executable(PDsampleExpanding PDsampleExpanding.cpp )
target_link_libraries(PDsampleExpanding pdDeps)

add_executable(PDSampleSingleExpanding PDSampleSingleExpanding.cpp )
target_link_libraries(PDSampleSingleExpanding pdDeps)
 
 add_executable( Morphology Morphology.cpp )
target_link_libraries(Morphology pdDeps)

 
 install(TARGETS QTree BalTree    YuleTree BDTree YuleSim BDSim BBDTreeSim QtreeSim PD pdlist PDsampled PDsampleExpanding PDSampleSingleExpanding 	Morphology  DESTINATION ${CMAKE_INSTALL_PREFIX})
 