
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (ossim-demosaic)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

ADD_EXECUTABLE(ossim-demosaic
ossim-demosaic.cpp
)
INCLUDE_DIRECTORIES(/usr/local/include)

TARGET_LINK_LIBRARIES(ossim-demosaic ossim)


