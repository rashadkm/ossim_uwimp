//----------------------------------------------------------------------------
//
// License:  See top level LICENSE.txt file.
//
// File: ossim-demosaic.cpp
//
// Author:  Mohammed Rashad K.M
//
// Description: 
//
//----------------------------------------------------------------------------



#include <fstream>

#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>
#include <ossim/base/ossimConstants.h>
#include <ossim/base/ossimException.h>
#include <ossim/base/ossimNotify.h>
#include <ossim/init/ossimInit.h>
#include <iostream>
#include <ossim/imaging/ossimImageHandlerRegistry.h>
#include <ossim/imaging/ossimImageSourceSequencer.h>
#include <ossim/imaging/ossimTiffTileSource.h>
#include <ossim/imaging/ossimMemoryImageSource.h>
#include <ossim/imaging/ossimTiffWriter.h>
#include <ossim/imaging/ossimU16ImageData.h>
#include <ossim/imaging/ossimImageData.h>
#include <ossim/imaging/ossimJpegWriter.h>

using namespace std;



template <typename T> bool demosaic(const ossimImageData* t)
{
    ossimIpt size = t->getImageRectangle().size();


    int num_bands= 3;
    ossimRefPtr<ossimMemoryImageSource> memSource = new ossimMemoryImageSource();
    ossimRefPtr<ossimImageData> imdata = new ossimImageData(0,OSSIM_FLOAT32,num_bands);
    int x_out = size.x/2;
    int y_out = size.x/2;
    int size_out = x_out * y_out;
    imdata->setWidth(x_out);
    imdata->setHeight(y_out);   
    imdata->initialize();


    ossim_float32* band1 = (ossim_float32*) imdata->getBuf(0);
    ossim_float32* band2 = (ossim_float32*) imdata->getBuf(1);
    ossim_float32* band3 = (ossim_float32*) imdata->getBuf(2);


    int rcount=0,bcount=0,gcount=0,g1count=0,g2count=0;
    ossim_float32 *red, *blue, *g1, *g2, *green;


    red   = new ossim_float32[size_out];
    green = new ossim_float32[size_out];
    blue  = new ossim_float32[size_out];
    g1    = new ossim_float32[size_out];
    g2    = new ossim_float32[size_out];
    int index = 0;
    
    T *buf = (T *)t->getBuf(0);

    for(int i=0;i< size.y;i++)
    {
        for(int j = 0; j<size.x;j++)
        {
            if(i%2 == 0)
            {
                if(j%2) 
                {
                    red[rcount++] = (ossim_float32)buf[index];
                }
                else
                {
                    g1[g1count++] =(ossim_float32)buf[index];
                }
            }
            else
            {
                if(j%2) 
                {
                    g2[g2count++] = (ossim_float32)buf[index];
                }
                else
                {
                    blue[bcount++] = (ossim_float32)buf[index];
                }        
          
            }
            index++;  
        }
    }
    


    for(int rx = 0; rx < size_out; rx++)
    {
        band1[rx] = red[rx];
    }

    for(int gx = 0; gx < size_out; gx++)
    {
        ossim_float32 g = (g1[gx] + g2[gx])*0.5;
        band3[gx] = g;
    }

    for(int bx = 0; bx < size_out; bx++)
    {
        band2[bx] = blue[bx];
    }




    memSource->setImage(imdata);

    ossimRefPtr<ossimTiffWriter> tiff = new ossimTiffWriter();

    tiff->setFilename("/code/data/out.tif");
    tiff->setOutputImageType("tiff_tiled_band_separate");
    tiff->connectMyInputTo(memSource.get());
    tiff->execute();
    tiff = 0;

/*
    ossimRefPtr<ossimJpegWriter> jpeg = new ossimJpegWriter();
    jpeg->setFilename("out.jpg");
    jpeg->setOutputImageType("jpg");
    jpeg->connectMyInputTo(memSource.get());
    jpeg->execute();
    jpeg = 0;
*/   




delete red;
delete green;
delete blue;
delete g1;
delete g2;






   return false;
}

int main(int argc, char *argv[])
{
   ossimArgumentParser ap(&argc, argv);
   ossimInit::instance()->addOptions(ap);
   ossimInit::instance()->initialize(ap);

   try
   {
      if (ap.argc() != 2)
      {
         cout << "\nUsage: "<<ap.getApplicationName()<<" <image>"<<endl;
         return 1;
      }
      

 
      ossimImageHandlerRegistry* registry = ossimImageHandlerRegistry::instance();
      ossimRefPtr<ossimImageHandler> handler = registry->open(ossimFilename(argv[1]));
      if (!handler.valid())
      {
         cout<<"  Could not openimage at <"<< argv[1] << ">. Aborting..."<<endl;
         return 1;
      }

      ossimRefPtr<ossimTiffTileSource> src = dynamic_cast<ossimTiffTileSource*>(handler.get());
      
      ossimIrect extent = src->getBoundingRect();
      
      ossimRefPtr<ossimImageData> tile = src->getTile(extent);

      handler->close();
      
      if (tile.valid())
         demosaic<ossim_uint16>(tile.get());


      return 0;
   }
   catch (const ossimException& e)
   {
      ossimNotify(ossimNotifyLevel_WARN) << e.what() << std::endl;
      return 1;
   }

}
