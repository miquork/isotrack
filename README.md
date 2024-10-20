# isotrack
IsoTrack calibration for HCAL with optional depth-dependence

How to run this (e.g. on lxplus):
root -l -b -q mk_compile
root -l -b -q mk_IsoTrack.C

Can then also rebin and repeat plotting:
root -l -b -q rebinProfiles.C
root -l -b -q drawIsoTrack.C
root -l -b -q drawCorrFact.C
root -l -b -q drawCovariance.C

To-do:
- add analysis of residuals and pulls into drawCorrFact.C
- add timer to IsoTrack::Loop() for bigger runs
- add HO?
- improve pileup subtraction: cone->strip for HE?