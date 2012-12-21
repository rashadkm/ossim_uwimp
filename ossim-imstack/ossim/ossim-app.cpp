
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




#include <tiffio.h>


#include <fstream>
#include <iostream>

#include <ossim/base/ossimArgumentParser.h>
#include <ossim/base/ossimApplicationUsage.h>
#include <ossim/base/ossimConstants.h>
#include <ossim/base/ossimException.h>
#include <ossim/base/ossimNotify.h>

#include <ossim/init/ossimInit.h>

#include <ossim/imaging/ossimImageHandlerRegistry.h>
#include <ossim/imaging/ossimImageSourceSequencer.h>
#include <ossim/imaging/ossimTiffTileSource.h>
#include <ossim/imaging/ossimMemoryImageSource.h>
#include <ossim/imaging/ossimTiffWriter.h>
#include <ossim/imaging/ossimU16ImageData.h>
#include <ossim/imaging/ossimImageData.h>
#include <ossim/imaging/ossimJpegWriter.h>






//// Some core files that everyone should include
//#include "macros.h"
//#include "Exception.h"
//#include "Operation.h"
//#include "Image.h"

//#include "Geometry.h"
//#include "Color.h"
//#include "File.h"
//#include "WLS.h"
//#include "Arithmetic.h"
//#include "Stack.h"
//using namespace ImageStack;

unsigned int divisor = 1;
//using namespace std;

template<typename T>
inline T max(const T &a, const T &b) {
    if (a > b) { return a; }
    return b;
}

template<typename T>
inline T min(const T &a, const T &b) {
    if (a < b) { return a; }
    return b;
}

void ClampSource(ossimRefPtr<ossimImageData> imdata);

//void saveOssim(ossimRefPtr<ossimImageData> imdata);
void saveIm(ossimRefPtr<ossimImageData> imdata);


ossimRefPtr<ossimImageData> createOssimImageSource(ossimImageData *data)
{
    ossimRefPtr<ossimMemoryImageSource> memSource = new ossimMemoryImageSource();
    ossimRefPtr<ossimImageData> imdata = new ossimImageData(memSource.get(),OSSIM_FLOAT32,3);  
    imdata->setWidth(data->getWidth());
    imdata->setHeight(data->getHeight());
    imdata->initialize();
    return imdata;
}
template <typename T> ossimRefPtr<ossimImageData> prepareTile(ossimImageData* im)
{

    ossimIpt size = im->getImageRectangle().size();

    int nbands = im->getNumberOfBands();
    
    

    ossimRefPtr<ossimMemoryImageSource> memSource = new ossimMemoryImageSource();
    ossimRefPtr<ossimImageData> imdata = new ossimImageData(memSource.get(),OSSIM_FLOAT32,nbands);  
    imdata->setImageRectangle(im->getImageRectangle());
    imdata->initialize();

    ossim_uint32 num_pixels = (ossim_uint32) (size.x * size.y);
    
    float multiplier = 1.0f / divisor;




 ossim_float32* out_buffer = (ossim_float32*) imdata->getBuf(0);
    //for (int band=0; band < nbands; band++)
    //{
        T* in_buffer =(T*) im->getBuf(0);
            // std::cout << "\n" << in_buffer[0] << "\n";
		    //exit(1);        

int p = 0;
    for (int y = 0; y < size.y; y++) {
        for (int x = 0; x < size.x; x++) {

                out_buffer[x +y *size.x] = (float)in_buffer[x + y *size.x] * multiplier;
                
           // if(p < 24)
                //printf("val=%lf\n",out_buffer[size.x-1 + 0]);
            //p++;    
		    //exit(1);
    
        }
    }
    
    
    //printf("val=%lf\n",out_buffer[size.x]);
//    
//    for (int y = 0; y < im.height; y++) {
//        assert(TIFFReadScanline(tiff, buffer, y, 1) != -1,
//               "Failed reading scanline\n");
//        for (int x = 0; x < im.width; x++) {
//           
//                im(x, y, c) = ((float)buffer[x * im.channels + c]) * multiplier;
//                
//             std::cout << "\n" << buffer[x * im.channels + c] << "\n";
//		    exit(1);
//            
//        }
//    }    


    return imdata;
}

void Demosaic2(ossimImageData* im, int xoff, int yoff, bool awb)
{


  ossimIpt size = im->getImageRectangle().size();

    int nbands = im->getNumberOfBands();

    if(nbands!=1)
    {
        cout << "Mosaiced images should have a single channel\n";
        return;
    }

    ossim_uint32 num_pixels = (ossim_uint32) (size.x * size.y);
    

   ossim_float32* input = im->getFloatBuf(0);




awb = false;

int width = size.x;
int height = size.y;
int ss = width;

    ossimRefPtr<ossimImageData> imdata = createOssimImageSource(im);
    ossim_float32* band1 =(ossim_float32*) imdata->getFloatBuf(0);
    ossim_float32* band2 =(ossim_float32*) imdata->getFloatBuf(1);
    ossim_float32* band3 =(ossim_float32*) imdata->getFloatBuf(2);

//std::cout << height << " h:w " << width << std::endl;

    // This algorithm is adaptive color plane interpolation (ACPI)
    // by Runs Hamilton and Adams
    // (The Adams is not me)

   // make sure the image is of even width and height
//    if(im->getWidth() & 1 || im->getHeight() & 1)
    {
//im = Window(im, 0, 0, 0, (im.width >> 1) << 1, (im.width >> 1) << 1, im.frames);

    }

    if(awb)
    {
        // Step 1
        // auto white balance: make sure all channels have the same mean
        double sum[2][2] = {{0, 0}, {0, 0}};
        double maximum[2][2] = {{0, 0}, {0, 0}};
        for (int y = 0; y < size.y; y++) {
        for (int x = 0; x < size.x; x++) {

                double val =  (double) input[x +  y * width];
                //cout << val << " ";
                maximum[x & 1][y & 1] = std::max(maximum[x & 1][y & 1], val);
                sum[x & 1][y & 1] += val;
            }
        }

        // modded scale for habcam green scale to be base value
        double scale = sum[0][0]/maximum[0][0];
        double multiplier[2][2] = {{1.0/maximum[0][0], scale/sum[0][1]},
                                   {scale/sum[1][0],   scale/sum[1][1]}};
   for (int y = 0; y < size.y; y++) {
        for (int x = 0; x < size.x; x++) {
         
                input[x + y * ss] = (float)(input[x + y * ss] * multiplier[x & 1][y & 1]);
            }
        }

	    // added NHV
        // Clamp::apply(im, 0, 1);

    }



    //assert(im.channels == 1, "Mosaiced images should have a single channel\n");

    //Image out(im.width, im.height, im.frames, 3);


    //printf("calling Demosaic::detrend(Window im, int xoff, int yoff\n");
    //Detrend::apply(im,xoff,yoff);
//pData[xloc + yloc * ss]
    // Step 2
    // Interpolate green, blending in horizontal or vertical directions depending on
    // gradient magnitudes. Add a correction factor based on the second derivative of
    // the color channel at that point.
    // Ie, calculate |dI/dx| + |d2I/dx2| and |dI/dy| + |d2I/dy2|, and interpolate
    // horizontally or vertically depending on which is smaller
    // if they're both the same, use both
    for (int y = 2; y < height-2; y++) {
	for (int x = 2; x < width-2; x++) {
	
	
		if (((x + xoff) & 1) == ((y + yoff) & 1)) {
		    // GREEN IS KNOWN
		    // we already know green here, just use it
		    //out(x, y, t)[1] = im(x, y, t)[0];
		    band2[x + y * ss] = input[x + y * ss];
		    
		    //std::cout << input[x + y * ss] << "\n";
		    //exit(1);
		} else {
		    // RED OR BLUE IS KNOWN
		    // gather neighbouring greens
		    float left1 = input[(x-1) + y * ss] , right1 = input[(x +1) + y * ss];
//		    float left1 = im(x-1, y, t)[0], right1 = im(x+1, y, t)[0];
//		    float up1 = im(x, y-1, t)[0], down1 = im(x, y+1, t)[0];
	    float up1 = input[x + (y-1) * ss] , down1 = input[x + (y+1) * ss]	;	    
		    // gather neighbouring reds or blues
		    
//		    float here = im(x, y, t)[0];
		    float here = input[x + y * ss];
//		    float left2 = im(x-2, y, t)[0], right2 = im(x+1, y, t)[0];
    	    float left2 = input[(x-2) + y * ss] , right2 = input[(x +1) + y * ss]	;	    
		    
//		    float up2 = im(x, y-2, t)[0], down2 = im(x, y+2, t)[0];
            float up2 = input[x + (y-2) * ss] , down2 = input[x + (y+2) * ss];  

		    // decide which way to interpolate
		    // (divide laplacian by two because it's across twice the baseline)
                    // the correction terms have been removed because they look retarded
                    //  nhv  added correction terms back
		    float interpHoriz = fabs(right1 - left1) + fabs(2*here - right2 - left2)/2;
		    float interpVert  = fabs(up1    - down1) + fabs(2*here - up2    - down2)/2;
		    if (interpHoriz < interpVert) { // horizontally
                        //float colAverage = (left2 + right2)/2;
			//float correction = here - colAverage;
                       float colAverage = (left2 + right2)/2;
			float correction = here - colAverage;
			// only apply half the correction, because it's across twice the baseline
				//out(x, y, t)[1] = (left1 + right1)/2;// + correction/2;
			    band2[x + y * ss] = (left1 + right1)/2 + correction/2;
		    } else if (interpVert < interpHoriz) { // vertically
                        //float colAverage = (up2 + down2)/2;
			//float correction = here - colAverage;
                        float colAverage = (up2 + down2)/2;
			 float correction = here - colAverage;
			//out(x, y, t)[1] = (up1 + down1)/2;// + correction/2;
			  band2[x + y * ss] = (up1 + down1)/2 + correction/2;
		    } else { // both
			float colAverage = (up2 + down2 + left2 + right2)/4;
			float correction = here - colAverage;
			 band2[x + y * ss] = (left1 + up1 + right1 + down1)/4 + correction/2;
		    }
		}
	    }
	}
    

    // Step 3
    // to predict blue (or red) on top of green,
    // A) average the two neighbouring blue (or red) pixels
    // B) Calculate the error you would have made if you did the same thing to predict green
    // C) Correct blue (or red) by that error

    // Step 4
    // to predict blue on red or red on blue,
    // we have 4 neighbours, diagonally around us
    // use the same approach as step 2, but take diagonal derivatives and interpolate diagonally


	for (int y = 2; y < height-2; y++) {
	    for (int x = 2;x < width-2; x++) {
		if (((x + xoff) & 1) == ((y + yoff) & 1)) {
		    // GREEN IS KNOWN (step 3)

		    // figure out which of red/blue is horizontally interpolated
		    // and which is vertically interpolated
		    int horizChannel, vertChannel;
		    if ((y + yoff) & 1) {
			horizChannel = 2;
			vertChannel = 0;
		    } else {
			horizChannel = 0;
			vertChannel = 2;
		    }

		    // do the horizontal interpolation

		    // compute an average for the color input[(x-2) + y * ss]
//		    float colLeft = im(x-1, y, t)[0], colRight = im(x+1, y, t)[0];
		    float colLeft = input[(x-1) + y * ss], colRight = input[(x+1) + y * ss];
		    float colAverage = (colLeft + colRight)/2;
		    // compute the same average for green
		    //float greenLeft = out(x-1, y, t)[1], greenRight = out(x+1, y, t)[1], greenHere = out(x, y, t)[1];
		    
		    float greenLeft = band2[(x-1) + y * ss], greenRight = band2[(x+1) + y * ss], greenHere = band2[x + y * ss];
		    
		    float greenAverage = (greenLeft + greenRight)/2;
		    // see how wrong the green average was
		    float correction = greenHere - greenAverage;
		    // set the output to the average color plus the correction factor needed for green
		    //out(x, y, t)[horizChannel] = colAverage + correction;
		    if(horizChannel == 0)
    		    band1[x + y * ss] = colAverage + correction;
		    else
       		    band2[x + y * ss]= colAverage + correction;
		    // do the vertical interpolation
		    //float colUp = im(x, y-1, t)[0], colDown = im(x, y+1, t)[0];
		    //float greenUp = out(x, y-1, t)[1], greenDown = out(x, y+1, t)[1];
		    
		    float colUp = input[x + (y-1) * ss], colDown = input[x + (y+1) * ss];
		    float greenUp = band2[x + (y-1) * ss], greenDown = band2[x + (y+1) * ss];		    
		    colAverage = (colUp + colDown)/2;
		    greenAverage = (greenUp + greenDown)/2;
		    correction = greenHere - greenAverage;
		    //out(x, y, t)[vertChannel] = colAverage + correction;
		    if(vertChannel == 0)
    		    band1[x + y * ss] = colAverage + correction;
		    else
       		    band3[x + y * ss]= colAverage + correction;		    

		} else {
		    // RED OR BLUE IS KNOWN (step 4)

		    // figure out which channel is known exactly
		    int knownChannel, unknownChannel;
		    if ((y+yoff) & 1) {
			knownChannel = 2;
			unknownChannel = 0;
		    } else {
			knownChannel = 0;
			unknownChannel = 2;
		    }
		    
		    if(knownChannel == 0)
    		    band1[x + y * ss] = input[x + y * ss];
		    else
       		    band3[x + y * ss]= input[x + y * ss];			    

		    // set the known channel to the correct value
		    //out(x, y, t)[knownChannel] = im(x, y, t)[0];

		    // for the unknown channel, do diagonal interpolation
		    // u is up left, v is down right, s is up right, t is down left
		    // p is the channel to be predicted, g is green (already interpolated)
//		    float up = im(x-1, y-1, t)[0], ug = out(x-1, y-1, t)[1];
//		    float vp = im(x+1, y+1, t)[0], vg = out(x+1, y+1, t)[1];
//		    float sp = im(x+1, y-1, t)[0], sg = out(x+1, y-1, t)[1];
//		    float tp = im(x-1, y+1, t)[0], tg = out(x-1, t+1, t)[1];
//		    float greenHere = out(x, y, t)[1];
//		    


		    float up = input[(x-1) + (y-1) * ss], ug = band2[(x-1) + (y-1) * ss];
		    
		    float vp = input[(x+1) + (y+1) * ss], vg = band2[(x+1) + (y+1) * ss];
		    float sp = input[(x+1) + (y-1) * ss], sg = band2[(x+1) + (y-1) * ss];
		    float tp = input[(x-1) + (y+1) * ss], tg = band2[(x-1) + 2 * ss]; //FIXME
		    float greenHere = band2[x + y * ss];
		    
		    		    

		    float interpUV = fabs(vp - up) + fabs(2*greenHere - vg - ug);
		    float interpST = fabs(sp - tp) + fabs(2*greenHere - sg - tg);

		    if (interpUV < interpST) {
			float greenAverage = (ug + vg)/2;
			float correction = greenHere - greenAverage;
		    
		    if(unknownChannel == 0)
    		    band1[x + y * ss] = (up + vp)/2 + correction;
		    else
       		    band3[x + y * ss]= (up + vp)/2 + correction;	
       		    			
			//out(x, y, t)[unknownChannel] = (up + vp)/2 + correction;
		    } else if (interpST < interpUV) {
			float greenAverage = (sg + tg)/2;
			float correction = greenHere - greenAverage;
			//out(x, y, t)[unknownChannel] = (sp + tp)/2 + correction;
		    if(unknownChannel == 0)
    		    band1[x + y * ss] = (sp + tp)/2 + correction;
		    else
       		    band3[x + y * ss]= (sp + tp)/2 + correction;				
		    } else {
			float greenAverage = (ug + vg + sg + tg)/4;
			float correction = greenHere - greenAverage;
		    if(unknownChannel == 0)
    		    band1[x + y * ss] = (up + vp + sp + tp)/4 + correction;
		    else
       		    band3[x + y * ss]= (up + vp + sp + tp)/4 + correction;			
//			out(x, y, t)[unknownChannel] = (up + vp + sp + tp)/4 + correction;
		    }

		}
	    }
        }
    


    // Step 5
    // zero the margins, which weren't interpolated, to avoid annoying checkerboard there
    // we could also do some more basic interpolation, but the margins don't really matter anyway

	for (int y = 0; y < height; y++) {
	    if (y == 0 || y == 1 || y == width - 2 || y == width - 1) {
		for (int x = 0; x < width; x++) {
//		    out(x, y, t)[0] = out(x, y, t)[1] = out(x, y, t)[2] = 0;
		    band1[x + y * ss] = band2[x + y * ss] = band3[x + y * ss] = 0;
		}
	    } else {
//		for (int c = 0; c < 6; c++) {
            

//		    out(0, y, t)[0] = out(0, y, t)[1] = out(0, y, t)[2] = 0;
//		    out(1, y, t)[0] = out(1, y, t)[1] = out(1, y, t)[2] = 0;
//		    out(im.width-2, y, t)[0] = out(im.width-2, y, t)[1] = out(im.width-2, y, t)[2] = 0;
//		    out(im.width-1, y, t)[0] = out(im.width-1, y, t)[1] = out(im.width-1, y, t)[2] = 0;

            
		    band1[0 + y * ss] = band2[0 + y * ss] = band3[0 + y * ss] = 0;
		    band1[1 + y * ss] = band2[1 + y * ss] = band3[1 + y * ss] = 0;		    
//		    out(1, y, t)[0] = out(1, y, t)[1] = out(1, y, t)[2] = 0;

		    band1[width-2 + y * ss] = band2[width-2 + y * ss] = band3[width-2 + y * ss] = 0;
//		    out(im.width-2, y, t)[0] = out(im.width-2, y, t)[1] = out(im.width-2, y, t)[2] = 0;
		    band1[width-1 + y * ss] = band2[width-1 + y * ss] = band3[width-1 + y * ss] = 0;
//		    out(im.width-1, y, t)[0] = out(im.width-1, y, t)[1] = out(im.width-1, y, t)[2] = 0;
//		}
	    }
        }
    

//printf("%lf\n",band1[10000]);
    // Step 6
    // clamp the output between zero and one
    // ET 2008-08-28: This works quite poorly when dealing with HDR-composited mosaiced images.  Commenting out
    //Clamp::apply(out, 0, 1);

    //return out;



ClampSource(imdata);

       

  

saveIm(imdata);
//saveOssim(imdata);
}

void ClampSource(ossimRefPtr<ossimImageData> imdata)
{
    ossimIpt size = imdata->getImageRectangle().size();
    
    ossim_float32* band1 =(ossim_float32*) imdata->getFloatBuf(0);
    ossim_float32* band2 =(ossim_float32*) imdata->getFloatBuf(1);
    ossim_float32* band3 =(ossim_float32*) imdata->getFloatBuf(2);    
    
    
    ossim_uint32 num_pixels = (ossim_uint32) (size.x * size.y);
        ossim_float32 lower = 0, upper = 1;
        for (int x = 0; x < num_pixels ; x++) 
        {
            band1[x] = std::max(lower,band1[x]);
            band1[x] = std::min(upper,band1[x]);
        }
        for (int x = 0; x < num_pixels ; x++) 
        {
            band2[x] = std::max(lower,band2[x]);
            band2[x] = std::min(upper,band2[x]);
        }
        
        for (int x = 0; x < num_pixels ; x++) 
        {
            band3[x] = std::max(lower,band3[x]);
            band3[x] = std::min(upper,band3[x]);
        }
       

}
//void saveOssim(ossimRefPtr<ossimImageData> imdata)
//{

//    ossim_float32* band1 =(ossim_float32*) imdata->getFloatBuf(0);
//    ossim_float32* band2 =(ossim_float32*) imdata->getBuf(1);
//    ossim_float32* band3 =(ossim_float32*) imdata->getBuf(2);
//    ossimIpt size = imdata->getImageRectangle().size();



//float multiplier = 65535;
//ossim_uint32 num_pixels = (ossim_uint32) (size.x * size.y);

//    
//        for (int x = 0; x < num_pixels ; x++) 
//        {
//            band1[x] = band1[x] * multiplier;

//       
//            band2[x] = band2[x] * multiplier;
//            band3[x] = band3[x] * multiplier;
//       }
//       
//                   
//ossimRefPtr<ossimMemoryImageSource> memSource = new ossimMemoryImageSource();
//memSource->setImage(imdata.get());
//   ossimRefPtr<ossimTiffWriter> tiff = new ossimTiffWriter();
//    tiff->setFilename("/code/data/out_ossim.tif");
//    tiff->setOutputImageType("tiff_tiled_band_separate");
//    tiff->connectMyInputTo(memSource.get());
//    tiff->execute();
//    tiff = 0;
//} 

void saveIm(ossimRefPtr<ossimImageData> imdata)
{

    ossim_float32* band1 =(ossim_float32*) imdata->getBuf(0);
    ossim_float32* band2 =(ossim_float32*) imdata->getBuf(1);
    ossim_float32* band3 =(ossim_float32*) imdata->getBuf(2);

ossimIpt size = imdata->getImageRectangle().size();


ossimRefPtr<ossimMemoryImageSource> memSource = new ossimMemoryImageSource();
memSource->setImage(imdata);

   ossimRefPtr<ossimTiffWriter> tiff = new ossimTiffWriter();
    tiff->setFilename("/code/data/out_demosaic.tif");
    tiff->setOutputImageType("tiff_tiled_band_separate");
    tiff->connectMyInputTo(memSource.get());
    tiff->execute();
    tiff = 0;
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


   ossim_float32* buffer = static_cast<ossim_float32*>(im->getFloatBuf(0));

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






    std::cout << "\tAWB2:" 
              << (1.0/maximum[0][0]) /16 << " " 
              << (scale/sum[0][1]) /16   << " "
              << (scale/sum[1][0]) /16   << " "
              << (scale/sum[1][1]) /16   << std::fixed << std::setprecision( 6 ) << endl;

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

  

    ossimIrect extent = handler->getBoundingRect();
      
    ossimRefPtr<ossimImageData> t = handler->getTile(extent);


    if(handler->getClassName() == "ossimTiffTileSource")
    {

        switch(t->getScalarType())
        {
            case OSSIM_SINT8:
            case OSSIM_UINT8:
                divisor = 0x000000ff;
                cout << "uint8\n";
                break;   
     
            case OSSIM_SINT16:
            case OSSIM_UINT16:
                divisor = 0x0000ffff;
                cout << "uint16\n";
                break;

            case OSSIM_SINT32: 
            case OSSIM_UINT32:
                divisor = 0xffffffff;
                cout << "uint32\n";
                break;

            case OSSIM_FLOAT32:
            case OSSIM_FLOAT64:
                divisor = 1;
                break;

            default:
            std::cout << "divisor not valid" << endl;
        }
    }

    ossimRefPtr<ossimImageData> tile = prepareTile<ossim_uint16>(t.get());
    if(tile.valid())
    {
        //cout << "valid\n";
        wbalance(tile.get());


        bool awb = false;
        int xoff = 0, yoff = 0;

//ossimRefPtr<ossimMemoryImageSource> memSource = new ossimMemoryImageSource();
//memSource->setImage(tile.get());
//   ossimRefPtr<ossimTiffWriter> tiff = new ossimTiffWriter();
//cout << "here\n\n\n";
//    tiff->setFilename("/code/data/out.tif");
//    //tiff->setOutputImageType("tiff_tiled_band_separate");
//    tiff->connectMyInputTo(memSource.get());
//    tiff->execute();
//    tiff = 0;


Demosaic2(tile.get(),xoff,yoff,awb);

////    pop();
//  //  push(im);



    }


    ossimInit::instance()->finalize();

    return 0L;
}

