project("cmaple-simde")

## SIMDe was obtained from https://github.com/simd-everywhere/simde/tags
##  -- unzip the archive and copy the 'simde' subfolder here

# export include directory
set(SIMDE_INCLUDE_DIRECTORY
${PROJECT_SOURCE_DIR}
  CACHE
  INTERNAL "SIMDe include directories" FORCE
)