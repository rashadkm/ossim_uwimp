
CMAKE_MINIMUM_REQUIRED (VERSION 2.6)
PROJECT (ossim-demosaic)

SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

ADD_EXECUTABLE(ossim-app
ossim-app.cpp)

ADD_EXECUTABLE(ossim-cc
ossim-cc.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Alignment.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Arithmetic.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Calculus.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Color.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Complex.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Control.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Convolve.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/DFT.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Exception.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/File.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FileCSV.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FileFLO.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FileHDR.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FileJPG.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FilePNG.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FilePPM.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FileTGA.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FileTIFF.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FileTMP.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/FileWAV.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Filter.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/GaussTransform.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Geometry.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/HDR.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Image.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/LAHBPCG.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/LightField.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/main.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Network.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/NetworkOps.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Operation.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Paint.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Panorama.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Parser.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/PatchMatch.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Prediction.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Projection.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Stack.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Statistics.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/Wavelet.cpp
/code/dev/gsoc2012/ossim/c++/ossim-imstack/ossim/imagestack_nhv/src/WLS.cpp
)

INCLUDE_DIRECTORIES(/usr/local/include)
INCLUDE_DIRECTORIES(imagestack_nhv/src/)

TARGET_LINK_LIBRARIES(ossim-app ossim fftw3f png tiff)


