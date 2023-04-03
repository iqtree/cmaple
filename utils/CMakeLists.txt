add_library(utils
tools.cpp tools.h
timeutil.h
operatingsystem.cpp operatingsystem.h
gzstream.h gzstream.cpp
matrix.h
logstream.h logstream.cpp
)

target_link_libraries(utils)

if(ZLIB_FOUND)
  target_link_libraries(utils ${ZLIB_LIBRARIES})
else(ZLIB_FOUND)
  target_link_libraries(utils zlibstatic)
endif(ZLIB_FOUND)
