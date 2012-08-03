#ifndef IMAGESTACK_COLOR_H
#define IMAGESTACK_COLOR_H

class ColorMatrix : public Operation {
  public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, vector<float> matrix);
};

class ColorConvert : public Operation {
  public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, string from, string to);
    static Image rgb2hsv(Window im);
    static Image hsv2rgb(Window im);
    static Image rgb2y(Window im);
    static Image y2rgb(Window im);	
    static Image rgb2yuv(Window im);
    static Image yuv2rgb(Window im);
    static Image rgb2xyz(Window im);
    static Image xyz2rgb(Window im);
    static Image lab2xyz(Window im);
    static Image xyz2lab(Window im);
    static Image rgb2lab(Window im);
    static Image lab2rgb(Window im);

    static Image uyvy2yuv(Window im);
    static Image yuyv2yuv(Window im);

    static Image uyvy2rgb(Window im);
    static Image yuyv2rgb(Window im);
};

class Demosaic : public Operation {
public:
    void help();
    void parse(vector<string> args);
     static Image apply(Window win, int xoff, int yoff, bool awb);
};

class WhiteBalance : public Operation {
public:
        void help();
        void parse(vector<string> args);
        static void apply(Window im, int xoff, int yoff);
        static Image apply(Window im); //added by rashad
};

class HistoAdapt : public Operation {
public:
        void help();
        void parse(vector<string> args);
        static void apply(Window im, float hist_center, float luma, float chroma, float sig_spread, float sig_off);
};

class Detrend : public Operation {
public:
        void help();
        void parse(vector<string> args);
        static void apply(Window im, int xoff, int yoff);
        static void _detrend(Window im, int resample, float alpha, float lambda);
        static void _detrend(Window im, Window im2, int resample, float alpha, float lambda);    
};

class Detrend2 : public Operation {
public:
        void help();
        void parse(vector<string> args);
        static void apply(Window im, int xoff, int yoff, int log_space, float alpha, float lambda, float scale);
        static void _detrend(Window im, int resample, float alpha, float lambda, float scale);
        static void _detrend(Window im, Window im2, int resample, float alpha, float lambda, float scale);    
};

#endif
