#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(G4SBS)

#----------------------------------------------------------------------------
# Force use of c++11 standard
#
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
# Update header file
add_custom_target( _gitinfo ALL COMMAND ${PROJECT_SOURCE_DIR}/pullgitinfo.py ${PROJECT_SOURCE_DIR} COMMENT "Creating gitinfo.hh" )


#ROOT path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/src/dss2007 ${CMAKE_CURRENT_BINARY_DIR}/include)

ROOT_GENERATE_DICTIONARY(
    G4SBSDict                         # path to dictionary to generate
    "include/G4SBSRunData.hh"         # list of classes to process
    "include/G4SBSTextFile.hh"
    LINKDEF "include/G4SBSLinkDef.h"  # ROOT linkDef file
)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/G4SBS*.cc ${PROJECT_SOURCE_DIR}/src/dss2007/*.cc ${PROJECT_SOURCE_DIR}/include/Pythia6_tree.C)
file(GLOB headers ${CMAKE_CURRENT_BINARY_DIR}/include/gitinfo.hh ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/include/wiser_pion.h ${PROJECT_SOURCE_DIR}/src/cteq/*.h ${PROJECT_SOURCE_DIR}/src/dss2007/*.hh ${PROJECT_SOURCE_DIR}/include/Pythia6_tree.h)

file(GLOB cteqsources ${PROJECT_SOURCE_DIR}/src/cteq/*.c  )
file(GLOB cteqheaders ${PROJECT_SOURCE_DIR}/src/cteq/*.h )

file(GLOB libsources ${PROJECT_SOURCE_DIR}/src/G4SBSRunData.cc ${PROJECT_SOURCE_DIR}/src/G4SBSTextFile.cc)
file(GLOB libheaders ${CMAKE_CURRENT_BINARY_DIR}/include/gitinfo.hh ${PROJECT_SOURCE_DIR}/include/*.hh )

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

include(${Geant4_USE_FILE})
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_library(sbscteq STATIC ${cteqsources} ${cteqheaders})
SET_TARGET_PROPERTIES(sbscteq PROPERTIES COMPILE_FLAGS "-DCTEQ_TBL_PATH=\\\"${PROJECT_SOURCE_DIR}/cteq-tbls\\\"")

add_executable(g4sbs g4sbs.cc ${sources} ${headers} G4SBSDict.cxx)
add_library(g4sbsroot SHARED ${libsources} ${libheaders} G4SBSDict.cxx)

target_link_libraries(g4sbs ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} sbscteq )
target_link_libraries(g4sbsroot ${ROOT_LIBRARIES} )

option(WITH_G4SBS_GDML "Build g4sbs with GDML output" OFF)
if( WITH_G4SBS_GDML)
    set_target_properties(g4sbs 
    PROPERTIES COMPILE_FLAGS "-DG4SBS_USE_GDML")
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build g4sbs. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G4SBS_SCRIPTS
  vis.mac
  run_example.mac
  scripts/gen_1.5GeV2.mac
  scripts/gen_4GeV2.mac
  scripts/gen_7GeV2.mac
  scripts/gen_10GeV2.mac
  scripts/gmn_3.5GeV2.mac
  scripts/gmn_4.5GeV2.mac
  scripts/gmn_5.7GeV2.mac
  scripts/gmn_8.1GeV2.mac
  scripts/gmn_10.2GeV2.mac
  scripts/gmn_12GeV2.mac
  scripts/gmn_13.5GeV2.mac
  scripts/genrp_4.5GeV2.mac
  scripts/gep_5GeV2.mac
  scripts/gep_8GeV2.mac
  scripts/gep_12GeV2.mac
  #scripts/sidis_template.mac
  #scripts/make_SIDIS_scripts.C
  scripts/vis_sidis.mac
  scripts/sidis_background.mac
  scripts/vis_gep.mac
  scripts/vis_gmn.mac
  scripts/vis_a1n.mac
  scripts/vis_genrp.mac
  scripts/preinit_ckov_scint.mac
  scripts/preinit_ckov_noscint.mac
  scripts/preinit_nockov_scint.mac
  scripts/preinit_nockov_noscint.mac
  scripts/preinit_ckov_noscint_nocalorimeters.mac
  scripts/sidis_example_script.mac
  scripts/default_thresholds.mac
  database/ecal_map.txt
  scripts/transversity.mac
  scripts/a1n_bigbite30deg_disgen.mac
  scripts/default_thresholds.mac
  scripts/sidis_template.mac
  )

foreach(_script ${G4SBS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

set(G4SBS_ROOT_MACROS
root_macros/C16_analysis.C
root_macros/C16_tree.C
root_macros/C16_tree_with_flux.C
root_macros/DIS_rates_DIS.C
root_macros/DIS_rates_PYTHIA.C
root_macros/ECAL_background_dose_rate.C
root_macros/Fit_SBS_optics_and_spin.C
root_macros/G4SBS_optics_fit.C
root_macros/G4SBS_spin_transport.C
root_macros/Optics_A1n_SIDIS.C
root_macros/Pythia6_minbias_gen.C
root_macros/Pythia6_tree.C
root_macros/RICH_analysis.C
root_macros/SIDIS_DST.C
root_macros/SIDIS_reduced_data_file.C
root_macros/TriggerLogicECalGEP.C
root_macros/a1n_tree.C
root_macros/background_rates_GMN.C
root_macros/flux.C
root_macros/g4sbs_a1n_tree.C
root_macros/g4sbs_tree.C
root_macros/gep_optics_tree.C
root_macros/gep_tree.C
root_macros/gep_tree_July2015.C
root_macros/gep_tree_elastic.C
root_macros/gep_tree_with_spin.C
root_macros/gep_trigger_analysis.C
root_macros/gep_trigger_analysis_elastic.C
root_macros/gep_trigger_analysis_elastic_L2.C
root_macros/gep_trigger_analysis_pythia.C
root_macros/gep_trigger_analysis_pythia_L2.C
root_macros/gmn_tree.C
root_macros/rate_table.C
root_macros/rate_table_DST.C
root_macros/C16_tree.h
root_macros/C16_tree_with_flux.h
root_macros/Pythia6_tree.h
root_macros/SIDIS_DST.h
root_macros/a1n_tree.h
root_macros/g4sbs_a1n_tree.h
root_macros/g4sbs_tree.h
root_macros/gep_optics_tree.h
root_macros/gep_tree.h
root_macros/gep_tree_July2015.h
root_macros/gep_tree_elastic.h
root_macros/gep_tree_with_spin.h
root_macros/gmn_tree.h
)

foreach(_macro ${G4SBS_ROOT_MACROS} )
	       configure_file( ${PROJECT_SOURCE_DIR}/${_macro}
	       ${PROJECT_BINARY_DIR}/${_macro} COPYONLY )
endforeach()

configure_file( ${PROJECT_SOURCE_DIR}/rootlogon.C ${PROJECT_BINARY_DIR}/rootlogon.C )

set(DSS2007_GRIDS HLO.GRID HNLO.GRID KALO.GRID KNLO.GRID PILO.GRID PROLO.GRID PLO.GRID PNLO.GRID)

foreach(_grid ${DSS2007_GRIDS})
  configure_file( ${PROJECT_SOURCE_DIR}/src/dss2007/${_grid} ${PROJECT_BINARY_DIR}/${_grid} COPYONLY )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#add_custom_target(G4SBS DEPENDS g4sbs)

add_dependencies(g4sbs _gitinfo)
add_dependencies(g4sbsroot _gitinfo)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS g4sbs DESTINATION bin)


