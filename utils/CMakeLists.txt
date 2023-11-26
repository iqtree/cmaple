add_library(cmaple_utils
tools.cpp tools.h
timeutil.h
operatingsystem.cpp operatingsystem.h
gzstream.h gzstream.cpp
matrix.h
logstream.h logstream.cpp
)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    if(ZLIB_FOUND)
  		target_link_libraries(cmaple_utils PUBLIC OpenMP::OpenMP_CXX ${ZLIB_LIBRARIES})
	else(ZLIB_FOUND)
  		target_link_libraries(cmaple_utils PUBLIC OpenMP::OpenMP_CXX zlibstatic)
	endif(ZLIB_FOUND)
else(OpenMP_CXX_FOUND)
	if(ZLIB_FOUND)
  		target_link_libraries(cmaple_utils ${ZLIB_LIBRARIES})
	else(ZLIB_FOUND)
  		target_link_libraries(cmaple_utils zlibstatic)
	endif(ZLIB_FOUND)
endif(OpenMP_CXX_FOUND)
