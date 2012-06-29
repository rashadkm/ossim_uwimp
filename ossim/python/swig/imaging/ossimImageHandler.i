/*-----------------------------------------------------------------------------
Filename        : ossimImageHandler.i
Author          : Vipul Raheja
License         : See top level LICENSE.txt file.
Description     : Contains SWIG-Python of class ossimImageHandler 
-----------------------------------------------------------------------------*/

%module pyossim

%{

#include <ossim/imaging/ossimImageSource.h>
#include <ossim/imaging/ossimImageMetaData.h>
#include <ossim/base/ossimConstants.h>
#include <ossim/base/ossimNBandLutDataObject.h>
#include <ossim/base/ossimIrect.h>
#include <ossim/base/ossimFilename.h>
#include <ossim/base/ossimRefPtr.h>
#include <ossim/imaging/ossimFilterResampler.h>

#include <vector>
#include "numpy/arrayobject.h"
#include "Python.h"
%}

#ifndef TYPE_DATA
#define TYPE_DATA
#endif


%include "ossim/base/ossimConstants.h"

%include "ossim/imaging/ossimImageHandler.h"

%typemap(in,numinputs=1) (PyArrayObject *psArray)
{

  if ($input != NULL) // && PyArray_Check($input))
  {
      $1 = (PyArrayObject*)($input);
  }
  else
  {
      PyErr_SetString(PyExc_TypeError, "arg1 not a numpy array");
      SWIG_fail;
  }
}



%feature( "kwargs" ) getBufferAsNumPyArray;
%inline %{
template <typename T> void getBufferAsNumPyArray(const ossimImageData* tile,int band,ossim_uint32 num_pixels,PyArrayObject *psArray)
{      

    T* buf= (T *) tile->getBuf(band);
    for ( ossim_uint32 p=0; p<num_pixels; ++p)
    {
         //std::cout << (int)buf[p] << " ";

        psArray->data[p] =  buf[p];
  
    }

    
    return;
}

%}

%feature( "kwargs" ) ImageDataNumPy;
%inline %{
int ImageDataNumPy(const ossimImageData *tile,ossimScalarType stype, int band,ossimIpt size, PyArrayObject *psArray) 
{
    if( psArray->nd < 2 || psArray->nd > 3 )
    {
        printf("Illegal numpy array rank %d.\n", psArray->nd);
    }


    ossim_uint32 num_pixels = (ossim_uint32) (size.x * size.y);

    switch (stype)
    {
     case OSSIM_UINT8:
     case OSSIM_SINT8:
        getBufferAsNumPyArray<ossim_uint8>(tile, band, num_pixels, psArray);
        break;

     case OSSIM_UINT16:
     case OSSIM_SINT16:
     case OSSIM_USHORT11:
        getBufferAsNumPyArray<ossim_uint16>(tile, band,num_pixels, psArray);
        break;

     case OSSIM_UINT32:
     case OSSIM_SINT32:
        getBufferAsNumPyArray<ossim_uint32>(tile, band,num_pixels, psArray);
        break;

     case OSSIM_FLOAT32:
        getBufferAsNumPyArray<ossim_float32>(tile, band,num_pixels, psArray);
        break;

     case OSSIM_FLOAT64:
        getBufferAsNumPyArray<ossim_float64>(tile, band,num_pixels, psArray);
        break;

    default:
        printf("ossim datatype not supported by ossim\n");
        return -1;
    }

    return 0;
}



%}


%pythoncode %{


import numpy

dtypesmap = {  

PYOSSIM_UINT8           :   numpy.uint8,
PYOSSIM_SINT8           :   numpy.uint8,
PYOSSIM_UINT16          :   numpy.uint16,
PYOSSIM_SINT16          :   numpy.uint16,
PYOSSIM_UINT32          :   numpy.uint32,
PYOSSIM_SINT32          :   numpy.uint32,
PYOSSIM_FLOAT32         :   numpy.float32,
PYOSSIM_FLOAT64         :   numpy.float64,
PYOSSIM_CINT16          :   numpy.complex64,
PYOSSIM_CINT32          :   numpy.complex64,
PYOSSIM_CFLOAT32        :   numpy.complex64,
PYOSSIM_CFLOAT64        :   numpy.complex128

}


    
def find_dtype(dtype_):
    if isinstance(dtype_, type):
        if dtype_ == numpy.int8:
            return PYOSSIM_UINT8
        if dtype_ == numpy.complex64:
            return PYOSSIM_CFLOAT32
        
        for key, value in dtypesmap.items():
            if value == dtype_:
                return key
        return None
    else:
        try:
            return dtypesmap[dtype_]
        except KeyError:
            return None

def NumpyTypeToOssimType(numpy_dtype):
    if not isinstance(numpy_dtype, type):
        raise TypeError("Input must be a type")
    return find_dtype(numpy_dtype)

def OssimTypeToNumpyType(ossim_dtype):
    return find_dtype(ossim_dtype)
    

def ReadImageSourceNumPy(tile,ossim_dtype,band_index,numpy_dtype,size,buf_obj=None):
    buf_obj = numpy.empty([size.y,size.x], dtype = numpy_dtype)
    if ImageDataNumPy(tile.get(),ossim_dtype,band_index,size,buf_obj)!= 0:
        return None

    return buf_obj


def ossimImageSourceAsArray( handler,buf_array=None):

    extent = handler.getBoundingRect()
    tile = handler.getTile(extent)
    size = tile.getImageRectangle().size()
    nbands = tile.getNumberOfBands()
    ossim_dtype = tile.getScalarType()


    numpy_dtype = OssimTypeToNumpyType( ossim_dtype )
    if numpy_dtype == None:
        ossim_dtype = PYOSSIM_FLOAT32
        numpy_dtype = numpy.float32
    else:
        ossim_dtype = NumpyTypeToOssimType( numpy_dtype )


    array_list = []
    for band_index in range(0,nbands):
        if tile.valid():
            band_array = ReadImageSourceNumPy(tile, ossim_dtype, band_index, numpy_dtype, size)
            array_list.append( numpy.reshape( band_array, [1,size.y,size.x] ) )
        else:
            print "Invalid tile from numpy"

    return numpy.concatenate( array_list )

%}



