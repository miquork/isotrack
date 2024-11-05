# isotrack
IsoTrack calibration for HCAL with optional depth-dependence

How to run this (e.g. on lxplus):  
`git clone https://github.com/miquork/isotrack.git`  
`cd isotrack`  
`root -l -b -q mk_compile.C`  
`root -l -b -q mk_IsoTrack.C`  

Can then also rebin and repeat plotting:
`root -l -b -q rebinProfiles.C`  
`root -l -b -q drawIsoTrack.C`  
`root -l -b -q drawCorrFact.C`  
`root -l -b -q drawCovariance.C`  

To-do:
- add analysis of residuals and pulls into drawCorrFact.C
- add timer to IsoTrack::Loop() for bigger runs
- add jet veto map filtering (for ECAL holes, pixel holes)?
- add HO?
- improve pileup subtraction: cone->strip for HE?

HCAL-2024-10-13:
- document method and variables in slides
- merge depths 1+2 as suggested by Salavat, and to compare to Sunanda
   - S.B. changing AlCaRaw cuts from ECAL<2 GeV to 10 GeV to get some HE hadrons
   - direct ratio to Sunanda's results once available
- add GHI as soon as available

Versions
- v6: git code test, same as v5
- v7: switch off all quality cuts (trk, MIP, dR)
- v8: quality cuts back on
- v9: merge depths 1+2