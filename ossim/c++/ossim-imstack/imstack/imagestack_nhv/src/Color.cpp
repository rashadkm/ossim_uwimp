#include "main.h"
#include "Color.h"
#include "Arithmetic.h"
#include "WLS.h"
#include "GaussTransform.h"
#include "Geometry.h"
#include "NetworkOps.h"
#include "Statistics.h"
#include "LocalLaplacian.h"

void ColorMatrix::help() {
    printf("\n-colormatrix treats each pixel as a vector over its channels and multiplies\n"
	   "the vector by the given matrix. The matrix size and shape is deduced from the\n"
	   "number of arguments. The matrix is specified in column major order.\n\n"
	   "Converting rgb to grayscale:\n"
	   "  ImageStack -load color.tga -colormatrix 1 1 1 -scale 0.33333 -save gray.tga\n\n"
	   "Making orange noise:\n"
	   "  ImageStack -push 100 100 1 1 -noise -colormatrix 1 0.3 0 -save noise.tga\n\n"
	   "Making noise that varies between orange and blue:\n"
	   "  ImageStack -push 100 100 1 2 -noise -colormatrix 1 0.3 0 0 0 1 -save noise.tga\n\n");
}

void ColorMatrix::parse(vector<string> args) {
    assert(args.size() > 0, "-colormatrix requires arguments\n");

    vector<float> matrix(args.size());
    for (size_t i = 0; i < args.size(); i++) {
	matrix[i] = readFloat(args[i]);
    }

    Image im = apply(stack(0), matrix);
    pop();
    push(im);
}

Image ColorMatrix::apply(Window im, vector<float> matrix) {
    assert(matrix.size() % im.channels == 0,
	   "-colormatrix requires a number of arguments that is a multiple of the number of\n"
	   "channels of the current image\n");

    int inChannels = im.channels;
    int outChannels = (int)matrix.size() / inChannels;

    Image out(im.width, im.height, im.frames, outChannels);

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		for (int c = 0; c < im.channels; c++) {
		    for (int i = 0; i < out.channels; i++) {
			out(x, y, t)[i] += im(x, y, t)[c] * matrix[c * outChannels + i];
		    }
		}
	    }
	}
    }

    return out;
}




void ColorConvert::help() {
    printf("\n-colorconvert converts from one colorspace to another. It is called with two\n"
	   "arguments representing these colorspaces.\n\n"
	   "Allowable colorspaces are rgb, yuv, hsv, xyz, lab and y (luminance alone). grayscale,\n"
	   "gray, and luminance are synonyms for y, and hsb and hsl are synonyms for hsv.\n\n"
	   "Usage: ImageStack -load a.tga -colorconvert rgb hsv -scale 0.1 1 1\n"
	   "                  -colorconvert hsv rgb -save out.tga\n\n");
}

void ColorConvert::parse(vector<string> args) {
    assert(args.size() == 2, "-colorconvert requires two arguments\n");
    Image im = apply(stack(0), args[0], args[1]);
    pop();
    push(im);
}

Image ColorConvert::apply(Window im, string from, string to) {
    // check for the trivial case
    assert(from != to, "color conversion from %s to %s is pointless\n", from.c_str(), to.c_str());

    // unsupported destination color spaces
    if (to == "yuyv" ||
	to == "uyvy") {
	panic("Unsupported destination color space: %s\n", to.c_str());
    }

    // direct conversions that don't have to go via rgb
    if (from == "yuyv" && to == "yuv") {
	return yuyv2yuv(im);
    } else if (from == "uyvy" && to == "yuv") {
	return uyvy2yuv(im);
    } else if (from == "xyz" && to == "lab") {
	return xyz2lab(im);
    } else if (from == "lab" && to == "xyz") {
	return lab2xyz(im);
    } else if (from != "rgb" && to != "rgb") {
	// conversions that go through rgb
	Image halfway = apply(im, from, "rgb");
	return apply(halfway, "rgb", to);
    } else if (from == "rgb") { // from rgb
	if (to == "hsv" || to == "hsl" || to == "hsb") {
	    return rgb2hsv(im);
	} else if (to == "yuv") {
	    return rgb2yuv(im);
	} else if (to == "xyz") {
	    return rgb2xyz(im);
	} else if (to == "y" || to == "gray" ||
		   to == "grayscale" || to == "luminance") {
	    return rgb2y(im);
	} else if (to == "lab"){
	    return rgb2lab(im);
	} else {
	    panic("Unknown color space %s\n", to.c_str());
	}
    } else { //(to == "rgb")
	if (from == "hsv" || from == "hsl" || from == "hsb") {
	    return hsv2rgb(im);
	} else if (from == "yuv") {
	    return yuv2rgb(im);
	} else if (from == "xyz") {
	    return xyz2rgb(im);
	} else if (from == "y" || from == "gray" ||
		   from == "grayscale" || from == "luminance") {
	    return y2rgb(im);
	} else if (from == "lab") {
	    return lab2rgb(im);
	} else if (from == "uyvy") {
	    return uyvy2rgb(im);
	} else if (from == "yuyv") {
	    return yuyv2rgb(im);
	} else {
	    panic("Unknown color space %s\n", from.c_str());
	}
    }

    // keep the compiler happy
    return Image();

}

//conversions to and from lab inspired by CImg (http://cimg.sourceforge.net/)
Image ColorConvert::xyz2lab(Window im){
    #define labf(x)  ((x)>=0.008856?(powf(x,1/3.0)):(7.787*(x)+16.0/116.0))

    assert(im.channels == 3, "Image does not have 3 channels\n");

    Image out(im.width, im.height, im.frames, im.channels);

    //left in this form to allow for changes/fine-tuning
    float Xn = 1.0f/(0.412453 + 0.357580 + 0.180423);
    float Yn = 1.0f/(0.212671 + 0.715160 + 0.072169);
    float Zn = 1.0f/(0.019334 + 0.119193 + 0.950227);

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		float L, a, b;
		float X = im(x, y, t)[0];
		float Y = im(x, y, t)[1];
		float Z = im(x, y, t)[2];
		L = 1.16f * labf(Y*Yn) - 0.16f;
		a = 5.0f * (labf(X*Xn) - labf(Y*Yn));
		b = 2.0f * (labf(Y*Yn) - labf(Z*Zn));
		out(x, y, t)[0] = L;
		out(x, y, t)[1] = a;
		out(x, y, t)[2] = b;
	    }
	}
    }
    return out;

}

Image ColorConvert::lab2xyz(Window im){
    assert(im.channels == 3, "Image does not have 3 channels\n");
    Image out(im.width, im.height, im.frames, im.channels);

    float s = 6.0/29;

    float Xn = 0.412453 + 0.357580 + 0.180423;
    float Yn = 0.212671 + 0.715160 + 0.072169;
    float Zn = 0.019334 + 0.119193 + 0.950227;

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		float X, Y, Z, fy, fx, fz;
		float L = im(x, y, t)[0];
		float a = im(x, y, t)[1];
		float b = im(x, y, t)[2];
		fy = (L + 0.16f)/1.16f;
		fx = fy + a/5.0f;
		fz = fy - b/2.0f;

		if (fy > s) {
		    Y = Yn*(fy * fy * fy);
		} else {
		    Y = (fy - 16.0/116)*3*s*s*Yn;
		}
		if (fx > s) {
		    X = Xn*(fx * fx * fx);
		} else {
		    X = (fx - 16.0/116)*3*s*s*Xn;
		}
		if (fz > s) {
		    Z = Zn*(fz * fz * fz);
		} else {
		    Z = (fz - 16.0/116)*3*s*s*Zn;
		}

		out(x, y, t)[0] = X;
		out(x, y, t)[1] = Y;
		out(x, y, t)[2] = Z;
	    }
	}
    }

    return out;
}

Image ColorConvert::rgb2lab(Window im){
    assert(im.channels == 3, "Image does not have 3 channels\n");

    return xyz2lab(rgb2xyz(im));
}

Image ColorConvert::lab2rgb(Window im){
    assert(im.channels == 3, "Image does not have 3 channels\n");
    return xyz2rgb(lab2xyz(im));
}


Image ColorConvert::rgb2hsv(Window im) {
    assert(im.channels == 3, "Image does not have 3 channels\n");

    float mult = 1.0f / 6;

    Image out(im.width, im.height, im.frames, im.channels);

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		float minV, maxV, delta;
		float h, s, v;
		float r = im(x, y, t)[0];
		float g = im(x, y, t)[1];
		float b = im(x, y, t)[2];

		minV = min(r, g, b);
		maxV = max(r, g, b);
		v = maxV;

		delta = maxV - minV;

		if (delta != 0) {
		    s = delta / maxV;
		    if (r == maxV)      h = 0 + (g - b) / delta; // between yellow & magenta
		    else if (g == maxV) h = 2 + (b - r) / delta; // between cyan & yellow
		    else                h = 4 + (r - g) / delta; // between magenta & cyan
		    h *= mult;
		    if (h < 0) h++;
		} else {
		    // r = g = b = 0 so s = 0, h is undefined
		    s = 0;
		    h = 0;
		}

		out(x, y, t)[0] = h;
		out(x, y, t)[1] = s;
		out(x, y, t)[2] = v;
	    }
	}
    }

    return out;
}

Image ColorConvert::hsv2rgb(Window im) {
    assert(im.channels == 3, "Image does not have 3 channels\n");

    Image out(im.width, im.height, im.frames, im.channels);

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		float h = im(x, y, t)[0];
		float s = im(x, y, t)[1];
		float v = im(x, y, t)[2];

		if (s == 0) {
		    // achromatic (grey)
		    im(x, y, t)[0] = im(x, y, t)[1] = im(x, y, t)[2] = v;
		} else {

		    h *= 6;	// sector 0 to 5
		    int i = (int)h;
		    if (i == 6) i = 5;
		    float f = h - i;
		    float p = v * (1 - s);
		    float q = v * (1 - s * f);
		    float u = v * (1 - s * (1 - f));

		    switch (i) {
		    case 0:
			out(x, y, t)[0] = v;
			out(x, y, t)[1] = u;
			out(x, y, t)[2] = p;
			break;
		    case 1:
			out(x, y, t)[0] = q;
			out(x, y, t)[1] = v;
			out(x, y, t)[2] = p;
			break;
		    case 2:
			out(x, y, t)[0] = p;
			out(x, y, t)[1] = v;
			out(x, y, t)[2] = u;
			break;
		    case 3:
			out(x, y, t)[0] = p;
			out(x, y, t)[1] = q;
			out(x, y, t)[2] = v;
			break;
		    case 4:
			out(x, y, t)[0] = u;
			out(x, y, t)[1] = p;
			out(x, y, t)[2] = v;
			break;
		    default:  // case 5:
			out(x, y, t)[0] = v;
			out(x, y, t)[1] = p;
			out(x, y, t)[2] = q;
			break;
		    }
		}
	    }
	}
    }

    return out;
}

Image ColorConvert::rgb2y(Window im) {
    assert(im.channels == 3, "Image does not have 3 channels\n");

    Image out(im.width, im.height, im.frames, 1);

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		out(x, y, t)[0] = (im(x, y, t)[0] * 0.299f +
				   im(x, y, t)[1] * 0.587f +
				   im(x, y, t)[2] * 0.114f);
	    }
	}
    }

    return out;

}

Image ColorConvert::y2rgb(Window im) {
    assert(im.channels == 1, "Image does not have one channel\n");

    Image out(im.width, im.height, im.frames, 3);

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		out(x, y, t)[2] = out(x, y, t)[1] = out(x, y, t)[0] = im(x, y, t)[0];
	    }
	}
    }

    return out;
}

Image ColorConvert::rgb2yuv(Window im) {
    assert(im.channels == 3, "Image does not have 3 channels\n");
    Image out(im.width, im.height, im.frames, 3);

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		float r = im(x, y, t)[0];
		float g = im(x, y, t)[1];
		float b = im(x, y, t)[2];
		out(x, y, t)[0] = r *  0.299f + g *  0.587f + b *  0.114f;
		out(x, y, t)[1] = r * -0.169f + g * -0.332f + b *  0.500f + 0.5f;
		out(x, y, t)[2] = r *  0.500f + g * -0.419f + b * -0.0813f + 0.5f;
	    }
	}
    }

    return out;

}

Image ColorConvert::yuv2rgb(Window im) {
    assert(im.channels == 3, "Image does not have 3 channels\n");

    Image out(im.width, im.height, im.frames, 3);

    for (int t = 0; t < im.frames; t++) {
	for (int y = 0; y < im.height; y++) {
	    for (int x = 0; x < im.width; x++) {
		float Y = im(x, y, t)[0];
		float U = im(x, y, t)[1] - 0.5f;
		float V = im(x, y, t)[2] - 0.5f;
		out(x, y, t)[0] = Y + 1.4075f * V;
		out(x, y, t)[1] = Y - 0.3455f * U - 0.7169f * V;
		out(x, y, t)[2] = Y + 1.7790f * U;
	    }
	}
    }

    return out;
}

Image ColorConvert::rgb2xyz(Window im) {
    assert(im.channels == 3, "Image does not have 3 channels\n");

    float mat[9] = {0.412453, 0.212671, 0.019334,
		    0.357580, 0.715160, 0.119193,
		    0.180423, 0.072169, 0.950227};
    vector<float> matrix(mat, mat+9);

    return ColorMatrix::apply(im, matrix);

}

Image ColorConvert::xyz2rgb(Window im) {
    assert(im.channels == 3, "Image does not have 3 channels\n");

    float mat[9] = {3.240479, -0.969256, 0.055648,
		    -1.537150, 1.875992, -0.204043,
		    -0.498535, 0.041556, 1.057311};
    vector<float> matrix(mat, mat+9);

    return ColorMatrix::apply(im, matrix);

}

Image ColorConvert::uyvy2yuv(Window im) {
    assert(im.channels == 2,
	   "uyvy images should be stored as a two channel image where the second"
	   " channel represents luminance (y), and the first channel alternates"
	   " between u and v.\n");
    assert((im.width & 1) == 0,
	   "uyvy images must have an even width\n");

    Image out(im.width, im.height, im.frames, 3);
    float *outPtr = out(0, 0, 0);
    for (int t = 0; t < out.frames; t++) {
	for (int y = 0; y < out.height; y++) {
	    float *imPtr = im(0, y, t);
	    for (int x = 0; x < out.width; x+=2) {
		*outPtr++ = imPtr[1]; // Y
		*outPtr++ = imPtr[0]; // U
		*outPtr++ = imPtr[2]; // V
		*outPtr++ = imPtr[3]; // Y
		*outPtr++ = imPtr[0]; // U
		*outPtr++ = imPtr[2]; // V
		imPtr += 4;
	    }
	}
    }

    return out;
}

Image ColorConvert::yuyv2yuv(Window im) {
    assert(im.channels == 2,
	   "yuyv images should be stored as a two channel image where the first"
	   " channel represents luminance (y), and the second channel alternates"
	   " between u and v.\n");
    assert((im.width & 1) == 0,
	   "uyvy images must have an even width\n");

    Image out(im.width, im.height, im.frames, 3);
    float *outPtr = out(0, 0, 0);
    for (int t = 0; t < out.frames; t++) {
	for (int y = 0; y < out.height; y++) {
	    float *imPtr = im(0, y, t);
	    for (int x = 0; x < out.width; x+=2) {
		*outPtr++ = imPtr[0]; // Y
		*outPtr++ = imPtr[1]; // U
		*outPtr++ = imPtr[3]; // V
		*outPtr++ = imPtr[2]; // Y
		*outPtr++ = imPtr[1]; // U
		*outPtr++ = imPtr[3]; // V
		imPtr += 4;
	    }
	}
    }

    return out;
}

Image ColorConvert::uyvy2rgb(Window im) {
    return yuv2rgb(uyvy2yuv(im));
}

Image ColorConvert::yuyv2rgb(Window im) {
    return yuv2rgb(yuyv2yuv(im));
}

//TODO: rgb2lab
//lab2rgb
//xyz2lab
//lab2xyz




void Demosaic::help() {
    printf("\n-demosaic demosaics a raw bayer mosaiced image camera. It should be a one\n"
           "channel image. The algorithm used is adaptive color plane interpolation (ACPI).\n"
           "Demosaic optionally takes two or three arguments. Two arguments specify an offset\n"
           "of the standard bayer pattern in x and y. The presence of a third argument\n"
           "indicates that auto-white-balancing should be performed.\n\n"
           "Usage: ImageStack -load photo.dng -demosaic -save out.png\n"
           "       ImageStack -load raw.yuv -demosaic 0 1 awb -save out.png\n");
}

void Demosaic::parse(vector<string> args) {
    bool awb = false;
    int xoff = 0, yoff = 0;
    if (args.size() == 0) {
        awb = false;
    } else if (args.size() == 2) {
        xoff = readInt(args[0]);
        yoff = readInt(args[1]);
    } else if (args.size() == 3) {
        xoff = readInt(args[0]);
        yoff = readInt(args[1]);
        awb = true;
    } else {
        panic("-demosaic takes zero, two, or three arguments");
    }
    Image im = apply(stack(0), xoff, yoff, awb);
    pop();
    push(im);
}


Image Demosaic::apply(Window im, int xoff, int yoff, bool awb) {

    assert(im.channels == 1, "Mosaiced images should have a single channel\n");

    Image out(im.width, im.height, im.frames, 3);

    // This algorithm is adaptive color plane interpolation (ACPI)
    // by Runs Hamilton and Adams
    // (The Adams is not me)

    // make sure the image is of even width and height
    if (im.width & 1 || im.height & 1) {
	im = Window(im, 0, 0, 0, (im.width >> 1) << 1, (im.width >> 1) << 1, im.frames);
    }



   if (awb) {

        // Step 1
        // auto white balance: make sure all channels have the same mean
        double sum[2][2] = {{0, 0}, {0, 0}};
        double maximum[2][2] = {{0, 0}, {0, 0}};
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
#if 0
                	// nhv  make log of image
		    		double val = log(im(x, y, t)[0]+1);
                	im(x, y, t)[0] = val;
#else
		    		double val = im(x, y, t)[0];
#endif
                    maximum[x & 1][y & 1] = max(maximum[x & 1][y & 1], val);
                    sum[x & 1][y & 1] += val;
                }
            }
        }

        // modded scale for habcam green scale to be base value
        double scale = sum[0][0]/maximum[0][0];
        double multiplier[2][2] = {{1.0/maximum[0][0], scale/sum[0][1]},
                                   {scale/sum[1][0],   scale/sum[1][1]}};

//        printf("\tAWB: %.4f  %.4f  %.4f  %.4f\n", (1.0/maximum[0][0]) /16, (scale/sum[0][1]) /16,  (scale/sum[1][0]) /16,   (scale/sum[1][1]) /16 );
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    im(x, y, t)[0] = (float)(im(x, y, t)[0] * multiplier[x & 1][y & 1]);
                }
            }
        }
	    // added NHV
    	Clamp::apply(im, 0, 1);
    }

    //printf("calling Demosaic::detrend(Window im, int xoff, int yoff\n");
    //Detrend::apply(im,xoff,yoff);

    // Step 2
    // Interpolate green, blending in horizontal or vertical directions depending on
    // gradient magnitudes. Add a correction factor based on the second derivative of
    // the color channel at that point.
    // Ie, calculate |dI/dx| + |d2I/dx2| and |dI/dy| + |d2I/dy2|, and interpolate
    // horizontally or vertically depending on which is smaller
    // if they're both the same, use both
    for (int t = 0; t < im.frames; t++) {
	for (int y = 2; y < im.height-2; y++) {
	    for (int x = 2; x < im.width-2; x++) {
		if (((x + xoff) & 1) == ((y + yoff) & 1)) {
		    // GREEN IS KNOWN
		    // we already know green here, just use it
		    out(x, y, t)[1] = im(x, y, t)[0];
		} else {
		    // RED OR BLUE IS KNOWN
		    // gather neighbouring greens
		    float left1 = im(x-1, y, t)[0], right1 = im(x+1, y, t)[0];
		    float up1 = im(x, y-1, t)[0], down1 = im(x, y+1, t)[0];
		    // gather neighbouring reds or blues
		    float here = im(x, y, t)[0];
		    float left2 = im(x-2, y, t)[0], right2 = im(x+1, y, t)[0];
		    float up2 = im(x, y-2, t)[0], down2 = im(x, y+2, t)[0];

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
			out(x, y, t)[1] = (left1 + right1)/2 + correction/2;
		    } else if (interpVert < interpHoriz) { // vertically
                        //float colAverage = (up2 + down2)/2;
			//float correction = here - colAverage;
                        float colAverage = (up2 + down2)/2;
			 float correction = here - colAverage;
			//out(x, y, t)[1] = (up1 + down1)/2;// + correction/2;
			 out(x, y, t)[1] = (up1 + down1)/2 + correction/2;
		    } else { // both
			float colAverage = (up2 + down2 + left2 + right2)/4;
			float correction = here - colAverage;
			out(x, y, t)[1] = (left1 + up1 + right1 + down1)/4 + correction/2;
		    }
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
	for (int y = 2; y < im.height-2; y++) {
	    for (int x = 2; x < im.width-2; x++) {
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

		    // compute an average for the color
		    float colLeft = im(x-1, y, t)[0], colRight = im(x+1, y, t)[0];
		    float colAverage = (colLeft + colRight)/2;
		    // compute the same average for green
		    float greenLeft = out(x-1, y, t)[1], greenRight = out(x+1, y, t)[1], greenHere = out(x, y, t)[1];
		    float greenAverage = (greenLeft + greenRight)/2;
		    // see how wrong the green average was
		    float correction = greenHere - greenAverage;
		    // set the output to the average color plus the correction factor needed for green
		    out(x, y, t)[horizChannel] = colAverage + correction;

		    // do the vertical interpolation
		    float colUp = im(x, y-1, t)[0], colDown = im(x, y+1, t)[0];
		    float greenUp = out(x, y-1, t)[1], greenDown = out(x, y+1, t)[1];
		    colAverage = (colUp + colDown)/2;
		    greenAverage = (greenUp + greenDown)/2;
		    correction = greenHere - greenAverage;
		    out(x, y, t)[vertChannel] = colAverage + correction;

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

		    // set the known channel to the correct value
		    out(x, y, t)[knownChannel] = im(x, y, t)[0];

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



void Detrend::help() {
	printf("\n Detrend a Bayer Image before it is deconvoluted\n"
           "Detrend optionally takes two arguments specifying an offset\n"
           "of the standard bayer pattern in x and y.\n\n"
           "Usage: ImageStack -load photo.raw -detrend -demosaic -save out.png\n");
 }

void Detrend::parse(vector<string> args) {
    int xoff = 0, yoff = 0;
    if (args.size() == 0) {
        xoff = yoff = 0;
    } else if (args.size() == 2) {
        xoff = readInt(args[0]);
        yoff = readInt(args[1]);
    } else {
        panic("-detrend takes zero, two");
    }
    apply(stack(0), xoff, yoff);
}


// modifies in place passed image
// xoff, yoff  offsets from standard Bayer pattern
// these aren't important in this application but may be
// hence including them in the interface
void Detrend::apply(Window im, int xoff, int yoff) {

	 assert(im.channels == 1, "Mosaiced images should have a single channel\n");

	// should expose these to the parser
   	int resamp = 0;
   	float alpha = 2.50;                       // 5     .125   for 1/2 size image
   	float lambda = 0.125;

   //	printf("\nperforming operation detrend %d %.3f %.3f\n", resamp,alpha,lambda);

	/* remove low frequency illumination artifacts due to non uniform lighting */
	/* this is done on a per channel basis */

        //printf("Detrend::detrend(Window im, int xoff, int yoff\t");

	int width = im.width >> 1;
	int height = im.height >> 1;

	// storage for the demuxed Bayer Image channels
#if 1
   	Image tmp(width, height, 4, 1);
   	//    Window(Window im, int minx_, int miny_, int mint_, int width_, int height_, int frames_) {
   	Window red(tmp,0,0,0,width,height,1);
   	Window green1(tmp,0,0,1,width,height,1);
   	Window green2(tmp,0,0,2,width,height,1);
   	Window blue(tmp,0,0,3,width,height,1);
#else
   	Image red(width,height,1,1);
   	Image green1(width,height,1,1);
   	Image green2(width,height,1,1);
   	Image blue(width,height,1,1);
 #endif

   	// hardwired for 12 bit
   	//float scale = 1.0f/4095;
   	//float scale = 1.0f;

        //printf("disassemble\t");
	// demux bayer image into single color images
	for (int t = 0; t < im.frames; t++) {
        	for (int y = 0; y < im.height; y++) {
        	        int y2 = y >> 1;
                	for (int x = 0; x < im.width; x++) {
				if (((x + xoff) & 1) == ((y + yoff) & 1)) {
					if(x&1) {
						green1(x>>1,y2,t)[0] =  log(im(x, y, t)[0] + 1); // * scale;
					} else {
						green2(x>>1,y2,t)[0] =  log(im(x, y, t)[0] + 1); // * scale;
					}
				} else {
					if(x&1) {
						red(x>>1,y2,t)[0] =  log(im(x, y, t)[0] + 1); // * scale;
					} else {
						blue(x>>1,y2,t)[0] =  log(im(x, y, t)[0] + 1); // * scale;
                			}
            			}
            		}
            	}

                //printf("calling _detrend()\t");
		_detrend(red, resamp, alpha, lambda);
		_detrend(blue, resamp, alpha, lambda);
		//_detrend(green1,green2, resamp, alpha, lambda);
		_detrend(green1, resamp, alpha, lambda);
		_detrend(green2, resamp, alpha, lambda);

                //printf("reassemble\t");
	        // encode Bayer Image from individual color images
	        //scale = 4095;
	        for (int y = 0; y < im.height; y++) {
	                int y2 = y >> 1;
                	for (int x = 0; x < im.width; x++) {
				if (((x + xoff) & 1) == ((y + yoff) & 1)) {
					if(x&1) {
						im(x, y, t)[0] = expf(green1(x>>1,y2,t)[0]) - 1; // * scale;
					} else {
						im(x, y, t)[0] = expf(green2(x>>1,y2,t)[0]) - 1; // * scale;
					}
				} else {
					if(x&1) {
						im(x, y, t)[0] = expf(red(x>>1,y2,t)[0]) - 1; // * scale;
					} else {
						im(x, y, t)[0] = expf(blue(x>>1,y2,t)[0]) - 1; // * scale;
					}
                		}
                	}
            	}
        }
        //printf("detrend returns\n");
}

void Detrend::_detrend(Window im, int resamp,float alpha,float lambda) {
        //printf("Demosaic::_detrend(Window im)\n");
       	float max = 0.0;
        //int t = 0;

   	//Image im2 = im.copy();

   	Image im2(im.width,im.height,im.channels,im.frames,im.data);

   	// hardwired for now  really need a better filter as this can still cause halos
   	// maybe a thin plate spline with a good 'stiffness'  ??

   	// quick linear downsampling
   	//Downsample::apply(im2,(im2.width>>resamp),(im2.height>>resamp),1);


       	Image tmp = myWLS::apply(im2, alpha, lambda, 100, 0.01);

        //printf("\tbefore max\n");
        for (int y = 0; y < tmp.height; y++) {
                for (int x = 0; x < tmp.width; x++) {
                        float val = tmp(x,y)[0];
               		if (val > max) {
               			max = val;
               		}
               	}
          }
          //printf("\tafter max\n");

          // upsample using  high quality lanczos method back to original size
          //Resample::apply(tmp,im.width,im.height,1);
          //Upsample::apply(tmp,4,4,1);

          // add the inverse of the low frequency field to 'flatten' the low frequency field to
          // to the 'found maximum' in field
          for (int y = 0; y < tmp.height; y++) {
                for (int x = 0; x < tmp.width; x++) {
               		float val = tmp(x,y)[0] - max;
               		im(x, y)[0] -= val;
               	}
        }
        //printf("\treturns\n");
}


void Detrend::_detrend(Window im, Window im2, int resamp, float alpha, float lambda) {
        //printf("Demosaic::_detrend(Window im)\n");
       	float max = 0.0;
        //int t = 0;

   	//Image im2 = im.copy();
   	Image _im(im.width,im.height,im.channels,im.frames,im.data);

   	// hardwired for now  really need a better filter as this can still cause halos
   	// maybe a thin plate spline with a good 'stiffness'  ??

   	// could use a gaussian pyramid here
   	//Downsample::apply(_im,(_im.width>>resamp),(_im.height>>resamp),1);

       	Image tmp = myWLS::apply(_im, alpha, lambda, 100, 0.01);

        //printf("\tbefore max\n");
        for (int y = 0; y < tmp.height; y++) {
                for (int x = 0; x < tmp.width; x++) {
                        float val = tmp(x,y)[0];
               			if (val > max) {
               				max = val;
               			}
               	}
          }
          //printf("\tafter max\n");

          // resample to original size
          //Resample::apply(tmp,im.width,im.height,1);
          //Upsample::apply(tmp,4,4,1);

          // add the inverse of the low frequency field to 'flatten' the low frequency field to
          // to the 'found maximum' in field
          for (int y = 0; y < tmp.height; y++) {
             for (int x = 0; x < tmp.width; x++) {
              	float val = tmp(x,y)[0] - max;
               	im(x, y)[0] -= val;
            }
       	}

         // add the inverse of the low frequency field to 'flatten' the low frequency field to
         // to the 'found maximum' in field
       	for (int y = 0; y < tmp.height; y++) {
        	for (int x = 0; x < tmp.width; x++) {
              	float val = tmp(x,y)[0] - max;
              	im2(x, y)[0] -= val;
              }
        }
        //printf("\treturns\n");
}


void Detrend2::help() {
	printf("\n Detrend a Bayer Image before it is deconvoluted\n"
           "Detrend optionally takes two arguments specifying an offset\n"
           "of the standard bayer pattern in x and y.\n\n"
           "Usage: ImageStack -load photo.raw -detrend -demosaic -save out.png\n");
 }

void Detrend2::parse(vector<string> args) {
    int xoff = 0, yoff = 0, log_space = 0;
   	float alpha = 2.5; //5.0;
   	float lambda = 0.25; //0.25
   	float scale = 1.0;
    if (args.size() == 0) {
        xoff = yoff = log_space = 0;
    } else if (args.size() == 1) {
        log_space = 1;
    } else if (args.size() == 2) {
        xoff = readInt(args[0]);
        yoff = readInt(args[1]);
    } else if (args.size() == 3) {
        xoff = readInt(args[0]);
        yoff = readInt(args[1]);
        log_space = readInt(args[2]) != 0;
    } else {
        panic("-detrend takes zero, two");
    }
    // !!!  HARDWIRED FOR NOW  !!!
    //log_space = 0;
   	scale = 1.0;
    apply(stack(0), xoff, yoff, log_space, alpha, lambda, scale);
}


// modifies in place passed image
// xoff, yoff  offsets from standard Bayer pattern
// these aren't important in this application but may be
// hence including them in the interface
void Detrend2::apply(Window im, int xoff, int yoff, int log_space, float alpha, float lambda, float scale) {

	 assert(im.channels == 1, "Mosaiced images should have a single channel\n");

	// should expose these to the parser
   	int resamp = 0;

   //	printf("\nperforming operation detrend %d %.3f %.3f\n", resamp,alpha,lambda);

	/* remove low frequency illumination artifacts due to non uniform lighting */
	/* this is done on a per channel basis */

    printf("Detrend2::detrend( log_space %d, alpha %5.3f, lambda %5.3f, scale %5.3f )\n", log_space, alpha, lambda, scale);

	int width = im.width >> 1;
	int height = im.height >> 1;

	// storage for the demuxed Bayer Image channels
#if 1
   	Image tmp(width, height, 4, 1);
   	//    Window(Window im, int minx_, int miny_, int mint_, int width_, int height_, int frames_) {
   	Window red(tmp,0,0,0,width,height,1);
   	Window green1(tmp,0,0,1,width,height,1);
   	Window green2(tmp,0,0,2,width,height,1);
   	Window blue(tmp,0,0,3,width,height,1);
#else
   	Image red(width,height,1,1);
   	Image green1(width,height,1,1);
   	Image green2(width,height,1,1);
   	Image blue(width,height,1,1);
 #endif

   	if (log_space) {
	   	// Offset::apply(im,1.0f);
   		// Log::apply(im);
     	for (int t = 0; t < im.frames; t++) {
			for (int y = 0; y < im.height; y++) {
	    		for (int x = 0; x < im.width; x++) {
		    		im(x, y, t)[0] = log( im(x, y, t)[0] + 1.0f );
	    		}
			}
    	}
    }

    //printf("disassemble\t");
	// demux bayer image into single color images
	for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < im.height; y++) {
            int y2 = y >> 1;
            for (int x = 0; x < im.width; x++) {
				if (((x + xoff) & 1) == ((y + yoff) & 1)) {
					if(x&1) {
						green1(x>>1,y2)[0] =  im(x, y)[0];
					} else {
						green2(x>>1,y2)[0] =  im(x, y)[0];
					}
				} else {
					if(x&1) {
						red(x>>1,y2)[0] =  im(x, y)[0];
					} else {
						blue(x>>1,y2)[0] =  im(x, y)[0];
            		}
        		}
			}
    	}
    }

        //printf("calling _detrend()\t");
	_detrend(red, resamp, alpha, lambda, scale);
	//_detrend(green1,green2, resamp, alpha, lambda);
	_detrend(green1, resamp, alpha, lambda, scale);
	_detrend(green2, resamp, alpha, lambda, scale);
	_detrend(blue, resamp, alpha, lambda, scale);

        //printf("reassemble\t");
	 	 // encode Bayer Image from individual color images
	for (int t = 0; t < im.frames; t++) {
	    for (int y = 0; y < im.height; y++) {
	        int y2 = y >> 1;
            for (int x = 0; x < im.width; x++) {
				if (((x + xoff) & 1) == ((y + yoff) & 1)) {
					if(x&1) {
						im(x, y)[0] = green1(x>>1,y2)[0];
					} else {
						im(x, y)[0] = green2(x>>1,y2)[0];
					}
				} else {
					if(x&1) {
						im(x, y)[0] = red(x>>1,y2)[0];
					} else {
						im(x, y)[0] = blue(x>>1,y2)[0];
					}
                }
            }
        }
    }
    /*
    printf("Green1\n");
	Statistics::apply(green1);
    printf("Green2\n");
	Statistics::apply(green2);
    printf("Red\n");
	Statistics::apply(red);
    printf("Blue\n");
	Statistics::apply(blue);
    */
    if( log_space ) {
	   	// Exp::apply(im);
   		// Offset::apply(im, -1.0f);
     	for (int t = 0; t < im.frames; t++) {
			for (int y = 0; y < im.height; y++) {
	    		for (int x = 0; x < im.width; x++) {
		    		im(x, y, t)[0] = exp(im(x, y, t)[0]) - 1.0f;
	    		}
			}
    	}
    }

	//printf("detrend returns\n");
}

void Detrend2::_detrend(Window im, int resamp,float alpha,float lambda, float scale) {
//	return;
        //printf("Demosaic::_detrend(Window im)\n");
#define TESTING 1
#if TESTING ////  testing LL
    float max_ = 0.0;
#endif
	//int t = 0;

   	//Image im2 = im.copy();

	Image im2(im);

#if 0
	// HotPixel supression
    for (int t = 0; t < im.frames; t++) {
 	        for (int y = 1; y < im.height-1; y++) {
 	            for (int x = 1; x < im.width-1; x++) {
 	                for (int c = 0; c < im.channels; c++) {
 	                    float n1 = im(x-1, y, t)[c];
 	                    float n2 = im(x+1, y, t)[c];
 	                    float n3 = im(x, y-1, t)[c];
 	                    float n4 = im(x, y+1, t)[c];
 	                    float here = im(x, y, t)[c];
 	                    float maxn = max(max(n1, n2), max(n3, n4));
 	                    float minn = min(min(n1, n2), min(n3, n4));
 	                    if (here > maxn) here = maxn;
 	                    if (here < minn) here = minn;
 	                    im2(x, y, t)[c] = here;
 	                }
 	            }
 	        }
 	    }
#endif        //printf("\treturns\n");

   		// hardwired for now  really need a better filter as this can still cause halos
   		// maybe a thin plate spline with a good 'stiffness'  ??

   		// quick linear downsampling
   		if (resamp)
   			Downsample::apply(im2,(im2.width>>resamp),(im2.height>>resamp),1);

   	#if TESTING
       	Image tmp = WLS::apply(im2, 2.5, 0.25, 100, 0.025);
    #else
       	//-bilateral 0.25 25. 25. 0 grid
//       Window tmp = LocalLaplacian::apply(im2, 1.0,  2.0);
       Window tmp = LocalLaplacian::apply(im2, -1.0,  -2.0);
//       Window tmp = LocalLaplacian::apply(im2, 1.0,  -1.0);
//       Window tmp = LocalLaplacian::apply(im2, 1.0,  -0.5);
//       Bilateral::apply(im2, 25, 25, 1, 0.25, GaussTransform::AUTO);
///       Bilateral::apply(im2, 50, 50, 1, 1.1, GaussTransform::AUTO);
       //Window tmp(im2);
       // -nlmeans 1.0 6 50 0.02
	#endif

#if TESTING ////  testing LL
        //printf("\tbefore max\n");
        for (int y = 0; y < tmp.height; y++)
        {
        	for (int x = 0; x < tmp.width; x++)
            {
	            float val = tmp(x,y)[0];
   		      		if (val > max_)	max_ = val;
            }
        }
#endif
          //printf("\tafter max\n");
        if (resamp)
      	// upsample using  high quality lanczos method back to original size
	      	Resample::apply(tmp,im.width,im.height,1);
          	//Upsample::apply(tmp,4,4,1);

        // add the inverse of the low frequency field to 'flatten' the low frequency field to
    	// to the 'found maximum' in field
        for (int y = 0; y < tmp.height; y++) {
           for (int x = 0; x < tmp.width; x++) {
#if TESTING ////  testing LL
           		float val = max_ - tmp(x,y)[0];
          		im(x, y)[0] += (val);
#else
				im(x, y)[0] = tmp(x,y)[0];
#endif
}
        }

}


void Detrend2::_detrend(Window im, Window im2, int resamp, float alpha, float lambda, float scale) {
    //printf("Demosaic::_detrend(Window im)\n");
    float max = 0.0;
    //int t = 0;

   	//Image im2 = im.copy();
   	Image _im(im.width,im.height,im.channels,im.frames,im.data);

   	// hardwired for now  really need a better filter as this can still cause halos
   	// maybe a thin plate spline with a good 'stiffness'  ??

   	// could use a gaussian pyramid here
   	if (resamp)
   		Downsample::apply(_im,(_im.width>>resamp),(_im.height>>resamp),1);

    Image tmp = WLS::apply(_im, alpha, lambda, 100, 0.01);

    //printf("\tbefore max\n");
    for (int y = 0; y < tmp.height; y++) {
        for (int x = 0; x < tmp.width; x++) {
            float val = tmp(x,y)[0] * scale;
            tmp(x,y)[0] = val;
       		if (val > max) {
        		max = val;
            }
        }
    }
    //printf("\tafter max\n");

    if (resamp)
        Resample::apply(tmp,im.width,im.height,1);
        //Upsample::apply(tmp,4,4,1);

    // add the inverse of the low frequency field to 'flatten' the low frequency field to
    // to the 'found maximum' in field
    for (int y = 0; y < tmp.height; y++) {
        for (int x = 0; x < tmp.width; x++) {
            float val = tmp(x,y)[0] - max;
            im(x, y)[0] -= val;
        }
    }

    // add the inverse of the low frequency field to 'flatten' the low frequency field to
    // to the 'found maximum' in field
    for (int y = 0; y < tmp.height; y++) {
    	for (int x = 0; x < tmp.width; x++) {
        	float val = tmp(x,y)[0] - max;
            im2(x, y)[0] -= val;
        }
    }
    //printf("\treturns\n");
}


void WhiteBalance::help() {
	printf("\n WhiteBalance a Bayer Image before it is deconvoluted\n"
           "WhiteBalance optionally takes two arguments specifying an offset\n"
           "of the standard bayer pattern in x and y.\n\n"
           "Usage: ImageStack -load photo.raw -detrend -whitebalance -demosaic -save out.png\n");
 }

void WhiteBalance::parse(vector<string> args) {
    int xoff = 0, yoff = 0;
    if (args.size() == 0) {
        xoff = yoff = 0;
    } else if (args.size() == 2) {
        xoff = readInt(args[0]);
        yoff = readInt(args[1]);
    } else {
        panic("-whitebalance takes zero, two");
    }
    apply(stack(0), xoff, yoff);
}


// modifies in place passed image
// xoff, yoff  offsets from standard Bayer pattern
// these aren't important in this application but may be
// hence including them in the interface
void WhiteBalance::apply(Window im, int xoff, int yoff) {

	 assert(im.channels == 1, "Mosaiced images should have a single channel\n");

   //	printf("\nperforming operation detrend %d %.3f %.3f\n", resamp,alpha,lambda);

        // auto white balance: make sure all channels have the same mean
        double sum[2][2] = {{0, 0}, {0, 0}};
        double maximum[2][2] = {{0, 0}, {0, 0}};
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
#if 0
                	// nhv  make log of image
//		    		double val = log(log(im(x, y, t)[0]+1)+1);
//					printf("WhiteBalance using log values\n");
		    		double val = log(im(x, y, t)[0]+1);
                	im(x, y, t)[0] = val;
#else
		    		double val = im(x, y, t)[0];
#endif
                    maximum[x & 1][y & 1] = max(maximum[x & 1][y & 1], val);
                    sum[x & 1][y & 1] += val;
                }
            }
        }
#if 1  // calculate scale from first channel
        double scale = sum[0][0]/maximum[0][0];
        double multiplier[2][2] = {{1.0/maximum[0][0], scale/sum[0][1]},
                                   {scale/sum[1][0],   scale/sum[1][1]}};

/*
        double scale = sum[1][1]/maximum[1][1];
        double multiplier[2][2] = {{scale/sum[0][0], scale/sum[0][1]},
                                   {scale/sum[1][0],   1.0/maximum[1][1]}};
        // scale to avoid oversaturation  NHV
		multiplier[0][0] *= 0.9;
		multiplier[0][1] *= 0.9;
		multiplier[1][0] *= 0.9;
		multiplier[1][1] *= 0.9;
*/
#else
        // modded scale for habcam green scale to be base value
        double green_sum = 0.5 * (sum[0][0] + sum[1][1]);
        double green_max = 0.5 * (maximum[0][0] + maximum[1][1]);
        double scale = green_sum/green_max;
        // green seems to need a little boost  NHV
		green_max *= 0.975;
        double multiplier[2][2] = {{1.0/green_max, 0.95*scale/sum[0][1]},
                                   {scale/sum[1][0],   1.0/green_max}};
#endif
        printf("\tAWB: %.4f  %.4f  %.4f  %.4f\n", (1.0/maximum[0][0]) /16, (scale/sum[0][1]) /16,  (scale/sum[1][0]) /16,   (scale/sum[1][1]) /16 );
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                	// Grey World assumption
                    im(x, y, t)[0] = (float)(im(x, y, t)[0] * multiplier[x & 1][y & 1]);
                    // White Patch assumption
//                    im(x, y, t)[0] = (float)(im(x, y, t)[0] / maximum[x & 1][y & 1]);
                }
            }
        }

    	Clamp::apply(im, 0, 1);
}





////added by rashad
Image WhiteBalance::apply(Window im) {

	 assert(im.channels == 1, "Mosaiced images should have a single channel\n");

 
  int xoff = 0; //added by rashad
    int yoff  = 0; //added by rashad

   //	printf("\nperforming operation detrend %d %.3f %.3f\n", resamp,alpha,lambda);

        // auto white balance: make sure all channels have the same mean
        double sum[2][2] = {{0, 0}, {0, 0}};
        double maximum[2][2] = {{0, 0}, {0, 0}};
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
#if 0
                	// nhv  make log of image
//		    		double val = log(log(im(x, y, t)[0]+1)+1);
//					printf("WhiteBalance using log values\n");
		    		double val = log(im(x, y, t)[0]+1);
                	im(x, y, t)[0] = val;
#else
		    		double val = im(x, y, t)[0];
#endif
                    maximum[x & 1][y & 1] = max(maximum[x & 1][y & 1], val);
                    sum[x & 1][y & 1] += val;
                }
            }
        }
#if 1  // calculate scale from first channel
        double scale = sum[0][0]/maximum[0][0];
        double multiplier[2][2] = {{1.0/maximum[0][0], scale/sum[0][1]},
                                   {scale/sum[1][0],   scale/sum[1][1]}};

/*
        double scale = sum[1][1]/maximum[1][1];
        double multiplier[2][2] = {{scale/sum[0][0], scale/sum[0][1]},
                                   {scale/sum[1][0],   1.0/maximum[1][1]}};
        // scale to avoid oversaturation  NHV
		multiplier[0][0] *= 0.9;
		multiplier[0][1] *= 0.9;
		multiplier[1][0] *= 0.9;
		multiplier[1][1] *= 0.9;
*/
#else
        // modded scale for habcam green scale to be base value
        double green_sum = 0.5 * (sum[0][0] + sum[1][1]);
        double green_max = 0.5 * (maximum[0][0] + maximum[1][1]);
        double scale = green_sum/green_max;
        // green seems to need a little boost  NHV
		green_max *= 0.975;
        double multiplier[2][2] = {{1.0/green_max, 0.95*scale/sum[0][1]},
                                   {scale/sum[1][0],   1.0/green_max}};
#endif
        printf("\tAWB: %.4f  %.4f  %.4f  %.4f\n", (1.0/maximum[0][0]) /16, (scale/sum[0][1]) /16,  (scale/sum[1][0]) /16,   (scale/sum[1][1]) /16 );
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                	// Grey World assumption
                    im(x, y, t)[0] = (float)(im(x, y, t)[0] * multiplier[x & 1][y & 1]);
                    // White Patch assumption
//                    im(x, y, t)[0] = (float)(im(x, y, t)[0] / maximum[x & 1][y & 1]);
                }
            }
        }

    	return im;
}






void HistoAdapt::help() {
	printf("\n HistoAdapt an Image\n"
           "HistoAdapt optionally takes three arguments specifying hist center a Luma and a Chroma scale factor sigmoid spread sigmoid offset\n\n"
           "Usage: ImageStack -load photo.raw -detrend -whitebalance -demosaic -histoadapt 0.5 1.3 1.2 10 4 -save out.png\n");
 }

void HistoAdapt::parse(vector<string> args) {
    float hist_center = 0.5, luma = 1, chroma = 1, sig_spread = 10, sig_off = 5;
    if (args.size() == 1) {
        hist_center   = readFloat(args[0]);
    } else if (args.size() == 3) {
        hist_center   = readFloat(args[0]);
        luma             = readFloat(args[1]);
        chroma         = readFloat(args[2]);
     } else if (args.size() == 5) {
        hist_center   = readFloat(args[0]);
        luma             = readFloat(args[1]);
        chroma         = readFloat(args[2]);
        sig_spread   = readFloat(args[3]);
        sig_off          = readFloat(args[4]);
   } else {
        panic("-histoadapt takes one three or five");
    }
    apply(stack(0), hist_center, luma, chroma, sig_spread, sig_off);
}


    // im(x, y, t)[0] = 1.0f/(1.0f+expf(-(((val<0.0f?0.0f:(val>1.0f?1.0f:val))*10.0f)-4.0f)));
inline float my_sigmoid( float val, float sig_spread, float sig_off ) {
	val = val < 0.0f ? 0.0f : val > 1.0f ? 1.0f : val;
	val *= sig_spread;
	val -= sig_off;
	return( 1.0f / (1.0f + expf(-val)));
}

// modifies in place passed image
void HistoAdapt::apply(Window im, float hist_center, float luma, float chroma, float sig_spread, float sig_off) {

	 assert(im.channels >= 3, "Mosaiced images should have at least 3 channels\n");

   //	printf("\nperforming operation histoadapt %f %.3f %.3f\n", center, luma, chroma);
/*
//			' -evalchannels "[0]*(0.5/mean(0))" "[1]" "[2]" ' \
//			' -evalchannels "(([0]-mean(0))*1.3)+mean(0)" "1.2*[1]" "1.2*[2]" ' \
//			' -evalchannels "1.0/(1.0+exp(-(((val<0?0:(val>1?1:val))*10)-4.0)))" "[1]" "[2]"' \
*/
   	//Image im2 = im.copy();

   //im = ColorConvert::rgb2lab(im);



        double sum = 0;
        for (int t = 0; t < im.frames; t++)
            for (int y = 0; y < im.height; y++)
                for (int x = 0; x < im.width; x++)
                    sum += im(x, y, t)[0];

        float mean =  (sum / (im.height * im.width));
        float scale = hist_center / mean;
        printf("mean %f   scale %f   mn*scl %f\n", mean, scale, mean*scale);
        mean *= scale;

        for (int t = 0; t < im.frames; t++)
            for (int y = 0; y < im.height; y++)
                for (int x = 0; x < im.width; x++)
                {
                	float val = (((im(x, y, t)[0]*scale) - mean) * luma) + mean;
                    im(x, y, t)[0] = my_sigmoid(val, sig_spread, sig_off);
                    im(x, y, t)[1] = im(x, y, t)[1] * chroma;
                    im(x, y, t)[2] = im(x, y, t)[2] * chroma;
                }

//	Image im2 = ColorConvert::lab2rgb(im);
//	pop();
//	push(im2);
}

