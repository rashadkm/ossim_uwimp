#ifndef IMAGESTACK_LOCAL_LAPLACIAN_H
#define IMAGESTACK_LOCAL_LAPLACIAN_H
//#include "header.h"

class LocalLaplacian : public Operation {
public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, float alpha, float beta);
 protected:
    static Image pyramidDown(Window im);
    static Image pyramidUp(Window im, int w, int h, int f);
};

class myLocalLaplacian : public LocalLaplacian {
public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, float alpha, float beta);
};

//#include "footer.h"
#endif
