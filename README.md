# isotrack
IsoTrack calibration for HCAL with optional depth-dependence

How to run this (e.g. on lxplus):
root -l -b -q mk_compile
root -l -b -q mk_IsoTrack.C

Can then also repeat plotting:
root -l -b -q drawIsoTrack.C
root -l -b -q drawCorrFact.C
root -l -b -q drawCovariance.C