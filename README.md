# isotrack
IsoTrack calibration for HCAL with optional depth-dependence

How to run this on lxplus:  
`git clone https://github.com/miquork/isotrack.git`  
`cd isotrack`  
`root -l -b -q mk_compile.C`
`ls /eos/cms/store/group/dpg_hcal/my_file_directory/*.root > myeraEA.txt`
`nohup root -l -b -q 'mk_IsoTrack.C("myera","myversion")' > logs/log_myera_myversion.txt &`
`tail -f logs/log_myera_myversion.txt`

The code will typically run from 5 minutes to 1 hour depending on input file sizes. You are encouraged to split the processing into smaller pieces and simply `hadd` them together afterwards. This will also make it simpler to study time stability later. The files are only ~1 MB each.

The output file will be copied to
`rootfiles/IsoTrack_version_era_orig.root`

Copy this file from lxplus to local `rootfiles/` and run locally
`root -l -b -q 'mk_IsoTrack.C("myera","myversion")'` or
`root -l -b -q 'drawIsoTrack.C("myera", "myversion")'` and
`root -l -b -q 'drawCorrFact("","myera","myversion")'`
`root -l -b -q compareCorrFact.C` [edit file to add myera_myversion]

Can also rebin and then repeat plotting:  
`root -l -b -q rebinProfiles.C`  
`root -l -b -q 'drawIsoTrack.C("24GHI","25June06",false)'` // era,version,stat  
`root -l -b -q 'drawCorrFact.C("24GHI","25Jun06",false)'`   // era,version,stat
`root -l -b -q drawCovariance.C`  

To-do:
- add 2025 HCalPFCuts (https://indico.cern.ch/event/1475489/contributions/6216571/attachments/2961918/5209794/Run3Winter25_MC_Nov6_2024.pdf?#page=13)
- add explicit analysis of pileup biases
- add analysis of residuals and pulls into drawCorrFact.C
- add jet veto map filtering (for ECAL holes, pixel holes)?
- add HO?
- improve pileup subtraction: cone->strip for HE?

Versions
- v33: reference for IsoTrackV2 with data, noPU MC, PU MC
  - why CorrFact.root is missing depths? => prevcorr?
- v29: Sunanda-style, full corrs, iterated
- v24: Sunanda-style, full corrs
- v18: useSunanda with CalibCorr.C gains + phi
- v17: useSunanda with my gains
- v9: merge depths 1+2
- v8: quality cuts back on
- v7: switch off all quality cuts (trk, MIP, dR)
- v6: git code test, same as v5