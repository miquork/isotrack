#! /usr/bin/python
import os

##################
# Run 3 IOV list #
##################
#IOV_list= ['24C','24D','24E','24F']
#IOV_list= ['24CDEF']
IOV_list= ['24CDEFGHI']
#IOV_list= ['24C','24D','24E','24F','24G','24H','24I','24CDEFGHI']
version = 'lxplus_v20'

#os.system("rm *.so *.d *.pcm")
os.system("root -l -b -q mk_compile.C")
for iov in IOV_list:
    print("Process mk_IsoTrack.C for IOV "+iov+" and version "+version)
    #os.system("nohup root -l -b -q 'mk_IsoTrack.C(\""+iov+"\",\""+version+"\")' > logs/log_"+iov+"_"+version+".txt &")
    #os.system("nohup root -l -b -q 'mk_IsoTrack.C(\""+iov+"\",\""+version+"\")' > logs/log_"+iov+"_"+version+".txt")
    os.system("nohup root -l -b -q 'mk_IsoTrack.C(\"{}\", \"{}\")' > logs/log_{}_{}.txt".format(iov, version, iov, version))

#    os.system("fs flush")
#    wait()
#    time.sleep(sleep_time)
