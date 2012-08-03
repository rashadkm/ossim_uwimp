#ifndef WLS_H
#define WLS_H

class WLS : public Operation {
public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, float alpha, float lambda, int iter, float tolerance);
};


class myWLS : public WLS {
public:
    void help();
    void parse(vector<string> args);
//    static Image apply(Window im, float alpha, float lambda, float tolerance);
};

#endif
