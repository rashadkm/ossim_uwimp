

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

// DUMA is a useful library that puts an electric fence at the end of
// every allocation to detect overruns.
//#include <duma.h>

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <list>

using ::std::string;
using ::std::vector;
using ::std::pair;
using ::std::make_pair;
using ::std::map;
using ::std::list;



// Some core files that everyone should include
#include "macros.h"
#include "Exception.h"
#include "Operation.h"
#include "Image.h"

#include "Geometry.h"
#include "Color.h"
#include "File.h"
#include "WLS.h"
#include "Arithmetic.h"
#include "Stack.h"


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

static Image myrgb2y(Window im);
static Image myY(Image l);



template <typename T> ossimRefPtr<ossimImageData> prepareTile(ossimImageData* im,unsigned int divisor)
{

    ossimIpt size = im->getImageRectangle().size();

    int nbands = im->getNumberOfBands();

    ossimRefPtr<ossimMemoryImageSource> memSource = new ossimMemoryImageSource();
    ossimRefPtr<ossimImageData> imdata = new ossimImageData(memSource.get(),OSSIM_FLOAT32,nbands);  
    imdata->setImageRectangle(im->getImageRectangle());
    imdata->initialize();

    ossim_uint32 num_pixels = (ossim_uint32) (size.x * size.y);

    for (int band=0; band < nbands; band++)
    {
        T* in_buffer =(T*) im->getBuf(band);
        ossim_float32* out_buffer = (ossim_float32*) imdata->getBuf(band);
    
        for (ossim_uint32 p=0; p<num_pixels; ++p)
        {
            float multiplier_c = 1.0f / divisor;
            ossim_float32 val = in_buffer[p] * multiplier_c;
            out_buffer[p] = val;
        }
    }

    return imdata;
}



void wbalance(ossimImageData* im)
{

    // auto white balance: make sure all channels have the same mean
    double sum[2][2] = {{0, 0}, {0, 0}};
    double maximum[2][2] = {{0, 0}, {0, 0}};

    ossimIpt size = im->getImageRectangle().size();

    int nbands = im->getNumberOfBands();

    if(nbands!=1)
    {
        cout << "Mosaiced images should have a single channel\n";
        return;
    }

    ossim_uint32 num_pixels = (ossim_uint32) (size.x * size.y);


   ossim_float32* buffer = im->getFloatBuf(0);

    ossim_uint32 p=0;

    p = 0;
    for (int y = 0; y < size.y; y++) {
        for (int x = 0; x < size.x; x++) {

            double val =  (double) buffer[p++];
            maximum[x & 1][y & 1] = std::max(maximum[x & 1][y & 1], val);
            sum[x & 1][y & 1] += val;
        }
    }


    double scale = sum[0][0]/maximum[0][0];
    double multiplier[2][2] = {{1.0/maximum[0][0], scale/sum[0][1]},
                               {scale/sum[1][0],   scale/sum[1][1]}};






    std::cout << "\tAWB:" 
              << (1.0/maximum[0][0]) /16 << " " 
              << (scale/sum[0][1]) /16   << " "
              << (scale/sum[1][0]) /16   << " "
              << (scale/sum[1][1]) /16   << endl;

    p = 0;

    for (int y = 0; y < size.y; y++) {
        for (int x = 0; x < size.x; x++) {
            buffer[p] = (float)(buffer[p] * multiplier[x & 1][y & 1]);
            p++;
        }
    }


}
int main(int argc, char **argv) {


    ossimArgumentParser ap(&argc, argv);
    ossimInit::instance()->addOptions(ap);
    ossimInit::instance()->initialize(ap);

    const char *f = "/code/data/in.tif";


    ossimImageHandlerRegistry* registry = ossimImageHandlerRegistry::instance();
    ossimRefPtr<ossimImageHandler> handler = registry->open(ossimFilename(f));  //ossim_uint8
    if (!handler.valid())
    {
        cout<<"  Could not openimage at <"<< argv[1] << ">. Aborting..."<<endl;
        return 1;
    }

    //cout << handler->getClassName() << endl;

    unsigned int divisor = 1;

    ossimIrect extent = handler->getBoundingRect();
      
    ossimRefPtr<ossimImageData> t = handler->getTile(extent);


    if(handler->getClassName() == "ossimTiffTileSource")
    {

        switch(t->getScalarType())
        {
            case OSSIM_SINT8:
            case OSSIM_UINT8:
                divisor = 0x000000ff;
                break;   
     
            case OSSIM_SINT16:
            case OSSIM_UINT16:
                divisor = 0x0000ffff;
                break;

            case OSSIM_SINT32: 
            case OSSIM_UINT32:
                divisor = 0xffffffff;
                break;

            case OSSIM_FLOAT32:
            case OSSIM_FLOAT64:
                divisor = 1;
                break;

            default:
            std::cout << "divisor not valid" << endl;
        }
    }




 ossimIpt size = t->getImageRectangle().size();

int ss = size.x * size.y;


   
   ossim_uint16 *buf = t->getUshortBuf(0);
cout <<  buf[1] << endl;
cout << t->getScalarTypeAsString() << endl;

int j;
//  for(j = 0; j<ss;j++)    cout <<  (float)buf[j] << " ";



int width = size.x;
int height = size.y;
int bands = 1;
//cout << width << ":" << height << endl;

Image im(width, height, 1, bands);

j = 0;

float multiplier = 1.0f / divisor;

for (int y = 0; y < im.height; y++) 
{
    for (int x = 0; x < im.width; x++) 
    {
    	for (int c = 0; c < im.channels; c++) 
        {
	        im(x, y)[c] = (float)buf[j++] * multiplier;
	    }
    }
}


//int crop_startX,crop_startY,im_width,im_height,new_width,new_height;
//Crop *crop = new Crop();
//Image cropped = crop->apply(im,0,0,1280,width,height,1024);

Save *save = new Save();


WhiteBalance *wb = new WhiteBalance();
wb->apply(im,0,0);



Demosaic *demosaic = new Demosaic();
Image de_im = demosaic->apply(im,0,0,true);

//save->apply(de_im,"/code/tmp/ossim_imstack_demosaic.png","");

ColorConvert *cc = new ColorConvert();
Image cc_im = cc->apply(de_im,"rgb","lab");


//DUP
    Image &top = cc_im;
    Image newTop(top.width, top.height, top.frames, top.channels);
    memcpy(newTop.data, top.data, top.frames * top.width * top.height * top.channels * sizeof(float));

Downsample *ds = new Downsample();
int boxWidth = 2, boxHeight = 2, boxFrames = 1;
Image ds_im = ds->apply(newTop, boxWidth, boxHeight, boxFrames);



myWLS *wls =new myWLS();
float alpha = 0, lambda = 0, tolerance = 0.005;
int max_iter = 200;

alpha = 2.5;
lambda = 0.25;
tolerance = 0.025;
max_iter = 100;

Image rgb2y = myrgb2y(ds_im);
Image wls_im = wls->apply(rgb2y, alpha, lambda, max_iter, tolerance);
Image mywls = myY(wls_im);


Upsample *us = new Upsample();
boxWidth = 2, boxHeight = 2, boxFrames = 1;
Image us_im = us->apply(mywls, boxWidth, boxHeight, boxFrames);



	//Image tmp = stack_[1]; cc_im
	//stack_[1] = stack_[0]; us_im
	//stack_[0] = tmp;

//PULL 1

Image tmp = cc_im;
cc_im = us_im;
us_im = tmp;


Add *add = new Add();
add->apply(us_im,cc_im);

	Image tmp2 =  cc_im;
	cc_im  =  us_im;
	us_im  = tmp2;
    


//stack(0) == cc_im

Image new_im = cc_im;

Crop *crop2 = new Crop();
Image crop_im = crop2->apply(new_im, 4, 4, width-8, height-8);



float hist_center = 0.5, luma = 1, chroma = 1, sig_spread = 10, sig_off = 5;
HistoAdapt *ha = new HistoAdapt();
hist_center = 0.45;
luma = 1.25;
chroma = 2.0;
sig_spread = 10;
sig_off = 4;
ha->apply(crop_im, hist_center, luma, chroma, sig_spread, sig_off);


Image cc_im2 = cc->apply(crop_im,"lab","rgb");

Clamp *clamp = new Clamp();
clamp->apply(cc_im2, 0, 1);


save->apply(cc_im2,"/code/tmp/ossim_imstack_colorcorrect.png","");




/*

bin/ImageStack -time --load '/code/data/in.tif'   --crop 0 0 1360 1024  --whitebalance --demosaic 1 0  --colorconvert rgb lab  --dup --downsample --wls2 2.5 0.25 100 0.025 --upsample --pull 1 --add --crop 4 4 width-8 height-8 --histoadapt 0.45 1.25 2.0 10 4 --colorconvert lab rgb --clamp  -save /code/data/PNG_/in.png 

*/


    return 0;

}


static Image myY(Image l)
{
    Image out(l.width, l.height, 1, 3);
    float max = 0.0;
    for (int t = 0; t < out.frames; t++) {
		for (int y = 0; y < out.height; y++) {
	    	for (int x = 0; x < out.width; x++) {
	    		float val = l(x, y, t)[0];
	    		if(val > max)
	    			max = val;
	    	}
		}
    }
    for (int t = 0; t < out.frames; t++) {
		for (int y = 0; y < out.height; y++) {
	    	for (int x = 0; x < out.width; x++) {
//				out(x, y, t)[0] = (max - 1.2*l(x, y, t)[0]);  // 1.2*(max - l(x,y,t)[0]);
				out(x, y, t)[0] = (max - l(x, y, t)[0]);
				out(x, y, t)[1] =  0; //max-l(x, y, t)[1];
				out(x, y, t)[2] =  0; //max-l(x, y, t)[2];
	    	}
		}
    }
    return out;
}


static Image myrgb2y(Window im) {
    assert(im.channels == 3, "Image does not have 3 channels\n");

    Image out(im.width, im.height, im.frames, 1);

    for (int t = 0; t < im.frames; t++) {
		for (int y = 0; y < im.height; y++) {
	    	for (int x = 0; x < im.width; x++) {
				out(x, y, t)[0] = (im(x, y, t)[0]);
	    	}
		}
    }
    return out;
}
