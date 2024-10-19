// Run this script to compile CondFormats libraries. After this can easily run 
// root -l -b -q mk_GamHistosFill.C
// using R__LOAD_LIBRARY to load *.so
{

  // For IsoTrack code (v6.30/04)
  gROOT->ProcessLine(".L IsoTrack.C+g");

}
