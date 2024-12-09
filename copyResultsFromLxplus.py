#! /usr/bin/python
import os

##################
# Run 3 IOV list #
##################
#IOV_list= ['24C','24D','24E','24F']
#IOV_list= ['24CDEF']
#IOV_list= ['24CDEFGHI']
#IOV_list= ['24C','24D','24E','24F','24G','24H','24I','24CDEFGHI']
IOV_list= ['24CDEF','24GHI','24CDEFGHI']
version = 'lxplus_v29'

#os.system("rm *.so *.d *.pcm")
for iov in IOV_list:
    print("Copy files from lxplus for IOV "+iov+" and version "+version)
    print("First *.root to .");
    print("rsync -rutP voutila@lxplus.cern.ch:~/scratch0/isotrack/*{}_{}.root .".format(version,iov))
    os.system("rsync -rutP voutila@lxplus.cern.ch:~/scratch0/isotrack/*{}_{}.root .".format(version,iov))
    print("Then rootfiles/*.root to rootfiles/");
    print("rsync -rutP voutila@lxplus.cern.ch:~/scratch0/isotrack/rootfiles/*{}_{}*.root rootfiles/".format(version,iov))
    os.system("rsync -rutP voutila@lxplus.cern.ch:~/scratch0/isotrack/rootfiles/*{}_{}*.root rootfiles/".format(version,iov))
    print("pdf/*.pdf to pdf/");
    print("rsync -rutP voutila@lxplus.cern.ch:~/scratch0/isotrack/pdf/*{}_{}.pdf pdf/".format(version,iov))
    os.system("rsync -rutP voutila@lxplus.cern.ch:~/scratch0/isotrack/pdf/*{}_{}.pdf pdf/".format(version,iov))
    print("logs/*.txt to logs/");
    print("rsync -rutP voutila@lxplus.cern.ch:~/scratch0/isotrack/logs/*{}_{}.txt logs/".format(iov,version))
    os.system("rsync -rutP voutila@lxplus.cern.ch:~/scratch0/isotrack/logs/*{}_{}.txt logs/".format(iov,version))
#    os.system("fs flush")
#    wait()
#    time.sleep(sleep_time)
