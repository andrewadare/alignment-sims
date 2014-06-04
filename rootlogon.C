/*
   rootlogon.C
   This script is automatically invoked whenever ROOT is started.
   Add session-level configurations as needed.
*/

void rootlogon()
{
  Printf("Starting ROOT version %s.", gROOT->GetVersion());
  Printf("Running %s/rootlogon.C on %s.",
         gSystem->Getenv("PWD"),
         gSystem->HostName());

  if (gSystem->Getenv("TMPDIR"))
    gSystem->SetBuildDir(gSystem->Getenv("TMPDIR"));

  // Include paths
  // =============
  // There are 2 separate paths: one for ACLiC, and one for CINT or CLING.
  // 1. ACLiC include path
  const char *pwd = gSystem->WorkingDirectory();
  const char *inc = pwd; //gSystem->DirName(pwd); // Parent directory
  const char *geo = "/Users/adare/phenix/svxgeo";
  const char *mpd = "/Users/adare/millepede";

  gSystem->AddIncludePath(Form("-I%s ", inc));
  gSystem->AddIncludePath(Form("-I%s ", geo));
  gSystem->AddIncludePath(Form("-I%s ", mpd));

  // 2. Interpreter include path
  // Type .include (ROOT 5) or .I (ROOT 6) at the ROOT REPL to see a listing
#ifdef __CINT__
  gROOT->ProcessLine(Form(".include %s", inc));
  gROOT->ProcessLine(Form(".include %s", geo));
  gROOT->ProcessLine(Form(".include %s", mpd));
#endif
#ifdef __CLING__
  gROOT->ProcessLine(Form(".I %s", inc));
  gROOT->ProcessLine(Form(".I %s", geo));
  gROOT->ProcessLine(Form(".I %s", mpd));
#endif

  // Build shared libs using ACLiC
  gROOT->LoadMacro(Form("%s/Mille.cc+", mpd));
  gROOT->LoadMacro(Form("%s/SvxTGeo.C+", geo));
  gROOT->LoadMacro(Form("%s/SvxGeoTrack.C+", geo));
  gROOT->LoadMacro(Form("%s/SvxProj.C+", geo));

  return;
}
