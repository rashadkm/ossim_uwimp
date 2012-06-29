import sys
sys.path.append("../lib")
from pyossim import *
import numpy as np
def pyossimImageCompare(argc,argv):
    if argc !=2:
        print "Usage: python ossim-numpy-test.py <image1>";
        sys.exit(0)

    imgname = ossimFilename(argv[1]) 

    

    
    registry = ossimImageHandlerRegistry.instance()
    handler = registry.open(imgname)    
    if not handler.isOpen():
        print "Could not open first image at <" + imgname.c_str() +  ">. Aborting..."
        return 1
    
    #np.set_printoptions(threshold='nan')

    imdata = ossimImageSourceAsArray(handler)

    print imdata[0]
    print type(imdata)


    return 0   

       



if __name__ == "__main__":
    init = ossimInit.instance()
    init.initialize()
    ret = pyossimImageCompare(len(sys.argv),sys.argv)


