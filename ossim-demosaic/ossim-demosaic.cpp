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

using namespace std;



template <typename T> bool demosaic(const ossimImageData* t)
{
   ossimIpt size = t->getImageRectangle().size();

   int nbands = t->getNumberOfBands();


    ossim_uint32 num_pixels = (ossim_uint32) (size.x * size.y);
    
    int x = size.x;
    int y = size.y;



    
    int rcount=0,bcount=0,gcount=0,g1count=0,g2count=0;
    float *red, *blue, *g1, *g2, *green;


    red = new float[x*y];
    green = new float[x*y];
    blue = new float[x*y];
    g1 = new float[x*y];
    g2 = new float[x*y];
    int index = 0;
    
    for(int i=0;i< y;i++)
    {
        for(int j = 0; j<x;j++)
        {


    int band = 0;
    T* buf = (T*) t->getBuf(band);
          if(i%2 == 0)
          {
            if(j%2) 
            {
               red[rcount++] = (float)buf[index];
               //cout << buf[index] << " ";
               
            }
          }  
          
          if(i%2)
          {
            if(j%2 == 0) 
            {
               blue[bcount++] = (float)buf[index];
               //cout << buf[index] << " ";
            }
          }          
          if(i%2)
          {
            if(j%2)
            {
               g2[g2count++] = (float)buf[index];
               //cout << buf[index] << " ";
            }
          }
          else
          {
            if(j%2 == 0)
            {
               g1[g1count++] = (float)buf[index];
              // cout << buf[index] << " ";  
            }
          }
          
          
            index++;
        }
        
    }
    

        int gc = g1count + g2count;
        

        for(int gx = 0; gx < gc; gx++)
        {
            float g = (g1[gx] + g2[gx])*0.5;
            {
            green[gcount++] = g;
            }
        }
        ossimRefPtr<ossimImageData> redband = new ossimImageData();
        redband->copyNormalizedBufferToTile(1,red);

/*
 ossimTiffWriter *fileWriter = new ossimTiffWriter()

   fileWriter->connectMyInputTo(0, bm.get());
   

   
   fileWriter->open(output_file);
   
    if(fileWriter.valid())
    {
        try
        {
            fileWriter->execute();
        }
    }




ossimRefPtr<ossimMemoryImageSource> memSource = new ossimMemoryImageSource.h
memSource->setImage( ...... );  // pass a pointer to an ossimImageData

ossimRefPtr<ossimTiffTileSource> tiff = new ossimTiffTileSource();
tiff->setFilename( ...... ); // set your output name for file
tiff->setOutputImageType("tiff_tiled_band_separate");
tiff->connectMyInputTo(memSource.get());

tiff->execute();

tiff = 0;
*/

delete green;
delete red;
delete blue;







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


   cout << endl;
      return 0;
   }
   catch (const ossimException& e)
   {
      ossimNotify(ossimNotifyLevel_WARN) << e.what() << std::endl;
      return 1;
   }

}
