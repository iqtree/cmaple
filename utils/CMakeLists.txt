add_library(cmaple_utils
tools.cpp tools.h
timeutil.h
operatingsystem.cpp operatingsystem.h
gzstream.h gzstream.cpp
matrix.h
logstream.h logstream.cpp
)

if(ZLIB_FOUND)
  target_link_libraries(cmaple_utils ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  target_link_libraries(uticmaple_utilsls zlibstatic)
endif(ZLIB_FOUND)
