#! /usr/bin/python
import os

print("Copy code to lxplus");
print("rsync -rutP runAllIOVs.py mk_compile.C mk_IsoTrack.C IsoTrack.C drawIsoTrack.C compareCorrFact.C drawPileUp.C rebinProfiles.C drawCovariance.C hybridCorrFact.C voutila@lxplus.cern.ch:~/scratch0/isotrack/")
os.system("rsync -rutP runAllIOVs.py mk_compile.C mk_IsoTrack.C IsoTrack.C drawIsoTrack.C compareCorrFact.C drawPileUp.C rebinProfiles.C drawCovariance.C hybridCorrFact.C voutila@lxplus.cern.ch:~/scratch0/isotrack/")

#    os.system("fs flush")
#    wait()
#    time.sleep(sleep_time)
