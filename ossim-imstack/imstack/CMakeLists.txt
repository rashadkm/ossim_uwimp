
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (ossim-gsoc)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
SET(APPNAME ossim-colorcorrection)
ADD_EXECUTABLE(${APPNAME}
imagestack_nhv/src/Calculus.cpp
imagestack_nhv/src/Color.cpp
imagestack_nhv/src/Complex.cpp
imagestack_nhv/src/Control.cpp
imagestack_nhv/src/Convolve.cpp
imagestack_nhv/src/DFT.cpp
#imagestack_nhv/src/Display.cpp
#imagestack_nhv/src/DisplayWindow.cpp
imagestack_nhv/src/Exception.cpp
imagestack_nhv/src/Filter.cpp
imagestack_nhv/src/GaussTransform.cpp
imagestack_nhv/src/Geometry.cpp
imagestack_nhv/src/HDR.cpp
imagestack_nhv/src/Image.cpp
imagestack_nhv/src/LAHBPCG.cpp
imagestack_nhv/src/LightField.cpp
imagestack_nhv/src/LocalLaplacian.cpp
imagestack_nhv/src/Arithmetic.cpp
imagestack_nhv/src/Alignment.cpp
imagestack_nhv/src/NetworkOps.cpp
imagestack_nhv/src/Network.cpp
imagestack_nhv/src/Operation.cpp
imagestack_nhv/src/Paint.cpp
#imagestack_nhv/src/Panorama.cpp
#imagestack_nhv/src/PatchMatch.cpp
imagestack_nhv/src/Parser.cpp
imagestack_nhv/src/Prediction.cpp
imagestack_nhv/src/Projection.cpp
imagestack_nhv/src/Stack.cpp
imagestack_nhv/src/Statistics.cpp
imagestack_nhv/src/Wavelet.cpp
imagestack_nhv/src/WLS.cpp
imagestack_nhv/src/extra.cpp
main.cpp)


INCLUDE_DIRECTORIES(imagestack_nhv/src/)


TARGET_LINK_LIBRARIES(${APPNAME} fftw3f ossim)


