

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

using namespace std;


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

void Demosaic(ossimImageData* im, int xoff, int yoff, bool awb)
{


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


  ossim_float32* buffer = im->getFloatBuf(0);

    ossim_uint32 p=0;







    // This algorithm is adaptive color plane interpolation (ACPI)
    // by Runs Hamilton and Adams
    // (The Adams is not me)

   // make sure the image is of even width and height
    if(im->getWidth() & 1 || im->getHeight() & 1)
    {
//im = Window(im, 0, 0, 0, (im.width >> 1) << 1, (im.width >> 1) << 1, im.frames);

    }

    if(awb)
    {
        // Step 1
        // auto white balance: make sure all channels have the same mean
        double sum[2][2] = {{0, 0}, {0, 0}};
        double maximum[2][2] = {{0, 0}, {0, 0}};


        p = 0;
        for (int y = 0; y < size.y; y++) {
            for (int x = 0; x < size.x; x++) {
                double val =  (double) buffer[p++];
                maximum[x & 1][y & 1] = std::max(maximum[x & 1][y & 1], val);
                sum[x & 1][y & 1] += val;
            }
        }

        // modded scale for habcam green scale to be base value
        double scale = sum[0][0]/maximum[0][0];
        double multiplier[2][2] = {{1.0/maximum[0][0], scale/sum[0][1]},
                                   {scale/sum[1][0],   scale/sum[1][1]}};

        p = 0;

        for (int y = 0; y < size.y; y++) {
            for (int x = 0; x < size.x; x++) {
                buffer[p] = (float)(buffer[p] * multiplier[x & 1][y & 1]);
                p++;
            }
        }

	    // added NHV
        // Clamp::apply(im, 0, 1);

    }








    Image out(im.width, im.height, im.frames, 3);



 


    //printf("calling Demosaic::detrend(Window im, int xoff, int yoff\n");
    //Detrend::apply(im,xoff,yoff);

    // Step 2
    // Interpolate green, blending in horizontal or vertical directions depending on
    // gradient magnitudes. Add a correction factor based on the second derivative of
    // the color channel at that point.
    // Ie, calculate |dI/dx| + |d2I/dx2| and |dI/dy| + |d2I/dy2|, and interpolate
    // horizontally or vertically depending on which is smaller
    // if they're both the same, use both

        for (int y = 2; y < size.y-2; y++) 
        {
            for (int x = 2; x < size.x-2; x++) 
            {




		        if (((x + xoff) & 1) == ((y + yoff) & 1)) 
                {
		            // GREEN IS KNOWN
		            // we already know green here, just use it
                    band2[band2_count++] = (float)buffer[p++]; 
		        } 
                else 
                {
		            // RED OR BLUE IS KNOWN
		            // gather neighbouring greens
                    float left1 = buffer[x-1 + y * size.y];
                    float right1 =buffer[x+1 + y * size.y];
                    float up1 =buffer[x + (y-1) * size.y];
                    float down1 =buffer[x + (y+1) * size.y];

		            // gather neighbouring reds or blues
		            float here = buffer[p++];

                    float left2 = buffer[x-2 + y * size.y];
                    float right2 =buffer[x+1 + y * size.y];
                    float up2 =buffer[x + (y-2) * size.y];
                    float down2 =buffer[x +(y+2) * size.y];


                    //buffer[y][x]
                    //buffer[x-1 + y * size.y]


		            // decide which way to interpolate
		            // (divide laplacian by two because it's across twice the baseline)
                    // the correction terms have been removed because they look retarded
                    //  nhv  added correction terms back
		            float interpHoriz = fabs(right1 - left1) + fabs(2*here - right2 - left2)/2;
		            float interpVert  = fabs(up1    - down1) + fabs(2*here - up2    - down2)/2;
		            if (interpHoriz < interpVert) 
                    { 
                        // horizontally
                        //float colAverage = (left2 + right2)/2;
			            //float correction = here - colAverage;
                        float colAverage = (left2 + right2)/2;
			            float correction = here - colAverage;
			            // only apply half the correction, because it's across twice the baseline
				        //out(x, y, t)[1] = (left1 + right1)/2;// + correction/2;

                        band2[band2_count++] = (left1 + right1)/2 + correction/2;
	
		            } 
                    else if (interpVert < interpHoriz) 
                    { 
                        // vertically
                        //float colAverage = (up2 + down2)/2;
			            //float correction = here - colAverage;
                        float colAverage = (up2 + down2)/2;
			            float correction = here - colAverage;
			            //out(x, y, t)[1] = (up1 + down1)/2;// + correction/2;
			            band2[band2_count++] = (up1 + down1)/2 + correction/2;
		            } 
                    else 
                    { 
                        // both
            			float colAverage = (up2 + down2 + left2 + right2)/4;
        			    float correction = here - colAverage;
        			    band2[band2_count++] = (left1 + up1 + right1 + down1)/4 + correction/2;
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

    for (int t = 0; t < im.frames; t++) {
	for (int y = 2; y < size.y-2; y++) 
    {
	    for (int x = 2; x < size.x-2; x++) 
        {
    		if (((x + xoff) & 1) == ((y + yoff) & 1)) 
            {
    		    // GREEN IS KNOWN (step 3)
    		    // figure out which of red/blue is horizontally interpolated
	    	    // and which is vertically interpolated
	    	    int horizChannel, vertChannel;
	    	    if ((y + yoff) & 1) 
                {
        			horizChannel = 2;
        			vertChannel = 0;
    		    } 
                else 
                {
        			horizChannel = 0;
        			vertChannel = 2;
    		    }

    		    // do the horizontal interpolation
    
    		    // compute an average for the color
 


    		    float colLeft = buffer[x-1 + y * size.y], colRight = buffer[x+1 + y * size.y];
    		    float colAverage = (colLeft + colRight)/2;
		        // compute the same average for green
		        float greenLeft = band2[x-1 + y * size.y], greenRight =  band2[x+1 + y * size.y], greenHere =  band2[x + y * size.y];
		        float greenAverage = (greenLeft + greenRight)/2;
		        // see how wrong the green average was
		        float correction = greenHere - greenAverage;
		        // set the output to the average color plus the correction factor needed for green
                if(horizChannel == 2)
    		        band3[x + y * size.y] = colAverage + correction;
                else if(horizChannel == 0)
    		        band1[x + y * size.y] = colAverage + correction;

    		    buffer[x + (y-1) * size.y] = colAverage + correction;
		        // do the vertical interpolation
		        float colUp =  buffer[x + (y-1) * size.y], colDown =  buffer[x + (y+1) * size.y];
 
 
		        float greenUp = band2[x + (y-1) * size.y], greenDown = band2[x + (y+1) * size.y];
		        colAverage = (colUp + colDown)/2;
		        greenAverage = (greenUp + greenDown)/2;
		        correction = greenHere - greenAverage;
                if(vertChannel == 2)
    		        band3[x + y * size.y] = colAverage + correction;
                else if(vertChannel == 0)
    		        band1[x + y * size.y] = colAverage + correction;

		    } 
            else 
            {
		        // RED OR BLUE IS KNOWN (step 4)

		        // figure out which channel is known exactly
		        int knownChannel, unknownChannel;
		        if ((y+yoff) & 1) 
                {
        			knownChannel = 2;
        			unknownChannel = 0;
    		    } 
                else 
                {
        			knownChannel = 0;
        			unknownChannel = 2;
    		    }

                if(knownChannel == 2)
    		        band3[x + y * size.y] = buffer[x + y * size.y];
                else if(knownChannel == 0)
    		        band1[x + y * size.y] = buffer[x + y * size.y];

		        // set the known channel to the correct value
	

    		    // for the unknown channel, do diagonal interpolation
    		    // u is up left, v is down right, s is up right, t is down left
    		    // p is the channel to be predicted, g is green (already interpolated)
    		    float up = im(x-1, y-1, t)[0], ug = out(x-1, y-1, t)[1];
    		    float vp = im(x+1, y+1, t)[0], vg = out(x+1, y+1, t)[1];
    		    float sp = im(x+1, y-1, t)[0], sg = out(x+1, y-1, t)[1];
    		    float tp = im(x-1, y+1, t)[0], tg = out(x-1, t+1, t)[1];
    		    float greenHere = out(x, y, t)[1];

		    float interpUV = fabs(vp - up) + fabs(2*greenHere - vg - ug);
		    float interpST = fabs(sp - tp) + fabs(2*greenHere - sg - tg);

		    if (interpUV < interpST) {
			float greenAverage = (ug + vg)/2;
			float correction = greenHere - greenAverage;
			out(x, y, t)[unknownChannel] = (up + vp)/2 + correction;
		    } else if (interpST < interpUV) {
			float greenAverage = (sg + tg)/2;
			float correction = greenHere - greenAverage;
			out(x, y, t)[unknownChannel] = (sp + tp)/2 + correction;
		    } else {
			float greenAverage = (ug + vg + sg + tg)/4;
			float correction = greenHere - greenAverage;
			out(x, y, t)[unknownChannel] = (up + vp + sp + tp)/4 + correction;
		    }

		}
	    }
        }
    }


    // Step 5
    // zero the margins, which weren't interpolated, to avoid annoying checkerboard there
    // we could also do some more basic interpolation, but the margins don't really matter anyway
    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    if (y == 0 || y == 1 || y == im.height - 2 || y == im.height - 1) {
		for (int x = 0; x < im.width; x++) {
		    out(x, y, t)[0] = out(x, y, t)[1] = out(x, y, t)[2] = 0;
		}
	    } else {
//		for (int c = 0; c < 6; c++) {
		    out(0, y, t)[0] = out(0, y, t)[1] = out(0, y, t)[2] = 0;
		    out(1, y, t)[0] = out(1, y, t)[1] = out(1, y, t)[2] = 0;
		    out(im.width-2, y, t)[0] = out(im.width-2, y, t)[1] = out(im.width-2, y, t)[2] = 0;
		    out(im.width-1, y, t)[0] = out(im.width-1, y, t)[1] = out(im.width-1, y, t)[2] = 0;
//		}
	    }
        }
    }


    // Step 6
    // clamp the output between zero and one
    // ET 2008-08-28: This works quite poorly when dealing with HDR-composited mosaiced images.  Commenting out
    Clamp::apply(out, 0, 1);

    return out;










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

    ossimRefPtr<ossimImageData> tile = prepareTile<ossim_uint16>(t.get(),divisor);
    if(tile.valid())
    {
        //cout << "valid\n";
        wbalance(tile.get());


        bool awb = false;
        int xoff = 0, yoff = 0;


        Demosaic(tile.get(),xoff,yoff,awb);

//    pop();
  //  push(im);



    }


    ossimInit::instance()->finalize();

    return 0L;
}

