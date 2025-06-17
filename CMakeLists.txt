cmake_minimum_required ( VERSION 3.12 )
project ( mountain LANGUAGES CXX Fortran )

set ( CMAKE_BUILD_TYPE Release )

# — Compiler flags (match your Makefile) —
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}   -O -Wall -fPIC -m64" )
set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O -Wall -fPIC -m64" )

# — Find dependencies —
find_package ( CGAL REQUIRED )
find_package ( ROOT REQUIRED )
include ( ${ROOT_USE_FILE} )

# — CERNLIB directory (mathlib, kernlib, packlib) —
link_directories ( /home/cupsoft/Products/CERNLIB/V2006/lib )

add_executable ( flat flat.cc music-sr.f )
target_link_libraries ( flat PRIVATE
  ${ROOT_LIBRARIES}
  mathlib
  kernlib
  packlib
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES} )

add_executable ( mountain
  mountain.cc
  music-sr.f
)

# — Link libraries: CGAL, ROOT, CERNLIB, and the Fortran runtime —
target_link_libraries ( mountain PRIVATE
  CGAL::CGAL
  ${ROOT_LIBRARIES}
  mathlib
  kernlib
  packlib
  ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}
)
