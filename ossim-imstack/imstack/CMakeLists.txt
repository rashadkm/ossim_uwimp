
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (ossim-gsoc)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

ADD_EXECUTABLE(ossim-app1
/code/imagestack_nhv/src/Calculus.cpp
/code/imagestack_nhv/src/Color.cpp
/code/imagestack_nhv/src/Complex.cpp
/code/imagestack_nhv/src/Control.cpp
/code/imagestack_nhv/src/Convolve.cpp
/code/imagestack_nhv/src/DFT.cpp
/code/imagestack_nhv/src/Display.cpp
/code/imagestack_nhv/src/DisplayWindow.cpp
/code/imagestack_nhv/src/Exception.cpp
/code/imagestack_nhv/src/File.cpp
/code/imagestack_nhv/src/FileEXR.cpp
/code/imagestack_nhv/src/FileCSV.cpp
/code/imagestack_nhv/src/FileHDR.cpp
/code/imagestack_nhv/src/FileJPG.cpp
/code/imagestack_nhv/src/FilePNG.cpp
/code/imagestack_nhv/src/FileTMP.cpp
/code/imagestack_nhv/src/FilePPM.cpp
/code/imagestack_nhv/src/FileTGA.cpp
/code/imagestack_nhv/src/FileTIFF.cpp
/code/imagestack_nhv/src/FileFLO.cpp
/code/imagestack_nhv/src/FileWAV.cpp
/code/imagestack_nhv/src/Filter.cpp
/code/imagestack_nhv/src/GaussTransform.cpp
/code/imagestack_nhv/src/Geometry.cpp
/code/imagestack_nhv/src/HDR.cpp
/code/imagestack_nhv/src/Image.cpp
/code/imagestack_nhv/src/LAHBPCG.cpp
/code/imagestack_nhv/src/LightField.cpp
/code/imagestack_nhv/src/LocalLaplacian.cpp
/code/imagestack_nhv/src/Arithmetic.cpp
/code/imagestack_nhv/src/Alignment.cpp
/code/imagestack_nhv/src/NetworkOps.cpp
/code/imagestack_nhv/src/Network.cpp
/code/imagestack_nhv/src/Operation.cpp
/code/imagestack_nhv/src/Paint.cpp
/code/imagestack_nhv/src/Panorama.cpp
/code/imagestack_nhv/src/PatchMatch.cpp
/code/imagestack_nhv/src/Parser.cpp
/code/imagestack_nhv/src/Prediction.cpp
/code/imagestack_nhv/src/Projection.cpp
/code/imagestack_nhv/src/Stack.cpp
/code/imagestack_nhv/src/Statistics.cpp
/code/imagestack_nhv/src/Wavelet.cpp
/code/imagestack_nhv/src/WLS.cpp
extra.cpp
main.cpp)

#INCLUDE_DIRECTORIES(/code/imagestack_nhv/include)
INCLUDE_DIRECTORIES(/code/imagestack_nhv/src)
INCLUDE_DIRECTORIES(/usr/include/SDL)
#ADD_DEFINITION(-D_GNU_SOURCE=1)
#ADD_DEPENDENCIES(ossim-app1 ${CONNECTOR})
TARGET_LINK_LIBRARIES(ossim-app1 jpeg png tiff fftw3f SDL ossim)


