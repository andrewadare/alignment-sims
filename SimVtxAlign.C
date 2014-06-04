
#include "VtxAlignmentUtils.h"
#include "SvxTGeo.h"
#include "SvxGeoTrack.h"
#include "SvxProj.h"
#include "Mille.h"

#include <TFile.h>
#include <TCanvas.h>
#include <TArrow.h>
#include <TRandom3.h>
#include <TGeoManager.h>
#include <TGeoTrack.h>
#include <TNtuple.h>
#include <TLatex.h>
#include <TArrow.h>
#include <TMarker.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TPolyLine.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <string>
#include <map>

using namespace std;

typedef vector<double> vecd;
typedef vector<SvxGeoTrack> geoTracks;

// Globals
const double BField = 0.0;
const int nTracks   = (int)1e5;
const char *pedeInputFile = "vtx-pede-input.bin";
const char *pedeSteerFile = "vtx-pede-steer.txt";
const char *pedeConstFile = "vtx-pede-const.txt";

void GenTrack(SvxTGeo *geo, const double BField, SvxGeoTrack &t);
void AddHitNoise(SvxGeoTrack &t, double xsigma, double zsigma);
void TrackLoop(geoTracks &tracks, TNtuple *hitTree = 0, TNtuple *trkTree = 0);
bool TrackOk(SvxGeoTrack &t);
void WriteConstFile(const char *filename, SvxTGeo *geo);
void WriteSteerFile(const char *filename, const char *binfile, const char *constfile);
int Label(int layer, int ladder, string coord);
void ParInfo(int label, int &layer, int &ladder, string &coord);
void GetLadderXYZ(SvxTGeo *tgeo, vecd &x, vecd &y, vecd &z);
int GetCorrections(const char *resFile, std::map<int, double> &mpc);

// Call DrawDiffs() after DrawXY()
TCanvas *DrawXY(SvxTGeo *geo, const char *name, const char *title, TString opt);
void DrawDiffs(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2);






void SimVtxAlign(int iter = 1)
{
  // No point in continuing if Millepede II is not installed...
  if (TString(gSystem->GetFromPipe("which pede")).IsNull())
  {
    Printf("\"which pede\" returns nothing. Exiting.");
    gSystem->Exit(-1);
  }

  TString pisaFileIn  = (iter==0) ? "geom/svxPISA.par.ideal" :
                        Form("geom/svxPISA.par.ideal.%d", iter - 1);
  TString pisaFileOut = Form("geom/svxPISA.par.ideal.%d", iter);
  TString inFileName  = Form("rootfiles/vtx-fake.%d.root", iter - 1);
  TString outFileName = Form("rootfiles/vtx-fake.%d.root", iter);


  TFile *inFile = 0;
  TFile *outFile = new TFile(outFileName.Data(), "recreate");
  const char *hitvars = "layer:ladder:sensor:xs:ys:zs:x:y:z:"
                        "xsigma:zsigma:dz:ds:trkid";
  TNtuple *ht0 = new TNtuple("ht0", "SvxGeoHit variables in ideal geometry",
                             hitvars);
  TNtuple *ht1 = new TNtuple("ht1", "Pre-alignment SvxGeoHit variables",
                             hitvars);
  TNtuple *ht2 = new TNtuple("ht2", "Post-alignment SvxGeoHit variables",
                             hitvars);
  TNtuple *trktree = new TNtuple("trktree", "Tracks from SimVtxAlign.C",
                                 "id:y0:z0:phi:theta");

  SvxTGeo *tgeo = new SvxTGeo;
  tgeo->ReadParFile(pisaFileIn.Data());
  tgeo->MakeTopVolume(100, 100, 100);
  tgeo->AddSensors();

  TGeoManager *mgr = tgeo->GeoManager();
  TGeoVolume *top = mgr->GetTopVolume();
  mgr->CloseGeometry();

  WriteConstFile(pedeConstFile, tgeo);
  WriteSteerFile(pedeSteerFile, pedeInputFile, pedeConstFile);

  // Original ladder positions
  vecd x0; vecd y0; vecd z0;
  GetLadderXYZ(tgeo, x0, y0, z0);

  Printf("Generating %d tracks...", nTracks);
  geoTracks tracks;
  while ((int)tracks.size() < nTracks)
  {
    // Generate a new track
    SvxGeoTrack track;
    GenTrack(tgeo, BField, track);
    AddHitNoise(track, 50e-4, 80e-4);

    // Ensure track has four hits, one per layer.
    if (TrackOk(track))
    {
      tracks.push_back(track);
      FillHitNTuple(track, ht0);
    }

  }

  if (iter==0)
  {
    // Misalignment: move B2L3 by +1 mm in z
    tgeo->TranslateLadder(2, 3, 0.0, 0.0, 0.1);

    // Update global hit positions to reflect misalignments.
    // Local hit positions (x,z on sensor) remain unchanged.
    for (unsigned int i=0; i<tracks.size(); i++)
      tracks[i].UpdateHits();
  }

  // Record ladder positions after applied misalignment
  vecd x1; vecd y1; vecd z1;
  GetLadderXYZ(tgeo, x1, y1, z1);

  TrackLoop(tracks, ht1, trktree);

  // Shell out to pede executable
  gSystem->Exec(Form("pede %s", pedeSteerFile));

  // Retrieve Millepede's corrections to the global parameters
  // Key is label, value is correction.
  map<int, double> mpc;
  GetCorrections("millepede.res", mpc);

  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<tgeo->GetNLadders(i); j++)
    {
      // Phi correction from ds
      tgeo->RotateLadderRPhi(i, j, mpc[Label(i,j,"s")]);
      // Longitudinal (z) correction
      tgeo->TranslateLadder(i, j, 0. ,0., mpc[Label(i,j,"z")]);
    }

  for (unsigned int i=0; i<tracks.size(); i++)
    tracks[i].UpdateHits();

  // Record positions after alignment
  vecd x2; vecd y2; vecd z2;
  GetLadderXYZ(tgeo, x2, y2, z2);

  Printf("Post-alignment refit...");
  for (unsigned int i=0; i<tracks.size(); i++)
  {
    double pars[4] = {0}; /* y0, z0, phi, theta */
    ZeroFieldResiduals(tracks[i], pars);
    FillHitNTuple(tracks[i], ht2);
  }

  // Draw changes in ladder positions
  const int NC = 4;
  TCanvas *c[NC];
  c[0] = DrawXY(tgeo, "vtx_xy", "VTX ladders", "L");
  c[1] = DrawXY(tgeo, "misalign", "Applied misalignment", "L,faint");
  DrawDiffs(x0,y0,z0,x1,y1,z1);
  c[2] = DrawXY(tgeo, "millepede_dp", "Alignment corrections", "L,faint");
  DrawDiffs(x1,y1,z1,x2,y2,z2);
  c[3] = DrawXY(tgeo, "align_error", "Alignment errors", "L,faint");
  DrawDiffs(x0,y0,z0,x2,y2,z2);

  Printf("Writing %s", pisaFileOut.Data());
  tgeo->WriteParFile(pisaFileOut.Data());

  Printf("Writing %s", outFileName.Data());
  ht0->Write("ht0");
  ht1->Write("ht1");
  ht2->Write("ht2");
  trktree->Write();
  for (int i=0; i<NC; i++)
    c[i]->Write();

  Printf("Done.");
  return;
}

bool
TrackOk(SvxGeoTrack &t)
{
  if (t.nhits == 4)
  {
    if (t.hits[0].layer==0 &&
        t.hits[1].layer==1 &&
        t.hits[2].layer==2 &&
        t.hits[3].layer==3)
    {
      return true;
    }
  }
  return false;
}

void
GenTrack(SvxTGeo *geo, const double BField, SvxGeoTrack &t)
{
  static TRandom3 ran3;
  SvxProj proj;
  proj.SetVerbosity(0);
  static long trackid = 0;

  t.nhits  = 0;
  t.charge = 0;
  t.vx     = 0.; //0.01*ran3.Gaus();
  t.vy     = 0.; //0.01*ran3.Gaus();
  t.vz     = 0.; //0.05*ran3.Gaus();
  t.mom    = ran3.Uniform(0.5, 5.0);
  t.phi0   = ran3.Uniform(0,TMath::TwoPi());
  t.the0   = ran3.Uniform(TMath::PiOver4(), 3*TMath::PiOver4());
  t.bfield = BField;

  if (false)
    Printf("Initialized track at (vx,vy,vz) %.2f,%.2f,%.2f with "
           "(phi0,the0) %.2f, %.2f charge %d, magField %.2f, mom %.2f",
           t.vx, t.vy, t.vz, t.phi0, t.the0, t.charge, t.bfield, t.mom);

  // Add hits to track
  proj.FindHitsFromVertex(t, geo);

  for (int j=0; j<t.nhits; j++)
    t.hits[j].trkid = trackid;
  trackid++;

  return;
}

void
AddHitNoise(SvxGeoTrack &t, double xsigma, double zsigma)
{
  static TRandom3 ran3;
  for (int j=0; j<t.nhits; j++)
  {
    t.hits[j].xsigma = xsigma;
    t.hits[j].zsigma = zsigma;

    double dx = t.hits[j].xsigma*ran3.Gaus();
    double dz = t.hits[j].zsigma*ran3.Gaus();

    t.hits[j].xs += dx;
    t.hits[j].zs += dz;

    t.hits[j].ds = dx;
    t.hits[j].dz = dz;

    SvxGeoHit hit = t.hits[j];
    double lxyz[3] = {hit.xs, hit.ys, hit.zs};
    double gxyz[3] = {0};
    TGeoNode *s = hit.node;
    if (s)
      s->GetMatrix()->LocalToMaster(lxyz, gxyz);

    t.hits[j].x = gxyz[0];
    t.hits[j].y = gxyz[1];
    t.hits[j].z = gxyz[2];
  }
  return;
}

void
TrackLoop(geoTracks &tracks, TNtuple *hitTree, TNtuple *trkTree)
{
  Printf("Track loop: fit, compute residuals, write to file...");
  Mille m(pedeInputFile);
  for (unsigned int i=0; i<nTracks; i++)
  {
    // Perform straight-line fit --> residuals, track parameters
    double pars[4] = {0}; /* y0, z0, phi, theta */
    ZeroFieldResiduals(tracks[i], pars);

    if (hitTree)
      FillHitNTuple(tracks[i], hitTree);

    if (trkTree)
    {
      // Track id and fit parameters for TTree
      float trkvars[5] = {0};
      trkvars[0] = tracks[i].hits[0].trkid;
      for (int j=0; j<4; j++)
        trkvars[j+1] = pars[j];
      trkTree->Fill(trkvars);
    }

    // Hit loop
    float dergl[1] = {1.0}; // Global derivatives
    int slabel[1] = {0};
    int zlabel[1] = {0};
    float sigma_s[4] = {32e-4, 39e-4, 42e-4, 25e-4};
    float sigma_z[4] = {52e-4, 63e-4, 66e-4, 40e-4};
    for (int j=0; j<tracks[i].nhits; j++)
    {
      SvxGeoHit hit = tracks[i].GetHit(j);
      float r = hit.x*hit.x + hit.y*hit.y;
      float sderlc[4] = {1.0,   r, 0.0, 0.0}; // Local derivatives
      float zderlc[4] = {0.0, 0.0, 1.0,   r}; // Local derivatives
      // float derlc[2] = {1.0, r}; // Local derivatives
      slabel[0] = Label(hit.layer, hit.ladder, "s");
      zlabel[0] = Label(hit.layer, hit.ladder, "z");

      m.mille(4, sderlc, 1, dergl, slabel, hit.ds, sigma_s[j]);
      m.mille(4, zderlc, 1, dergl, zlabel, hit.dz, sigma_z[j]);
    }
    m.end(); // Write residuals for this track & reset for next one
  }

  // Mille object must go out of scope for output file to close properly.
  return;
}

void
WriteConstFile(const char *filename, SvxTGeo *geo)
{
  Printf("Writing %s", filename);
  ofstream fs(filename);

  // Write linear constraints such that \sum_label w_label * p_label = c
  // where w is the weight to be applied for each term.
  //
  // The format is as follows:
  // Constraint c
  // label w_label
  // label w_label
  // ...

  // Parameter
  // 104 0.0 -1 !

  // Translation constraints for W then E arm
  fs << "Constraint 0.0" << endl;
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=0; ldr<geo->GetNLadders(lyr)/2; ldr++)
    {
      int label = Label(lyr, ldr, "z");
      float val = label==104 ? 0.0 : 1.0;
      fs << label << " " << val << endl;
    }
  fs << "Constraint 0.0" << endl;
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=geo->GetNLadders(lyr)/2; ldr < geo->GetNLadders(lyr); ldr++)
    {
      int label = Label(lyr, ldr, "z");
      float val = label==104 ? 0.0 : 1.0;
      fs << label << " " << val << endl;
    }

  // Shear constraints for W then E arm
  fs << "Constraint 0.0" << endl;
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=0; ldr<geo->GetNLadders(lyr)/2; ldr++)
    {
      int label = Label(lyr, ldr, "z");
      float val = label==104 ? 0.0 : geo->SensorRadius(lyr, ldr, 0);
      fs << label << " " << val << endl;
    }
  fs << "Constraint 0.0" << endl;
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=geo->GetNLadders(lyr)/2; ldr < geo->GetNLadders(lyr); ldr++)
    {
      int label = Label(lyr, ldr, "z");
      // float val = geo->SensorRadius(lyr, ldr, 0);
      float val = label==104 ? 0.0 : geo->SensorRadius(lyr, ldr, 0);
      fs << label << " " << val << endl;
    }

  // // Phi sum constraints for W then E arm
  // fs << "Constraint 0.0" << endl;
  // for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
  //   for (int ldr=0; ldr<geo->GetNLadders(lyr)/2; ldr++)
  //   {
  //     int label = Label(lyr, ldr, "z");
  //     float val = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0), TMath::TwoPi());
  //     if (label==104)
  //       val = 0;
  //     fs << label << " " << val << endl;
  //   }

  // fs << "Constraint 0.0" << endl;
  // for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
  //   for (int ldr=geo->GetNLadders(lyr)/2; ldr < geo->GetNLadders(lyr); ldr++)
  //   {
  //     int label = Label(lyr, ldr, "z");
  //     float val = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0), TMath::TwoPi());
  //     if (label==104)
  //       val = 0;
  //     fs << label << " " << val << endl;
  //   }

  // r*phi sum constraints for W then E arm
  fs << "Constraint 0.0" << endl;
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=0; ldr<geo->GetNLadders(lyr)/2; ldr++)
    {
      int label = Label(lyr, ldr, "z");
      float val = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0), TMath::TwoPi());
      val *= geo->SensorRadius(lyr, ldr, 0);
      if (label==104)
        val = 0;
      fs << label << " " << val << endl;
    }
  fs << "Constraint 0.0" << endl;
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=geo->GetNLadders(lyr)/2; ldr < geo->GetNLadders(lyr); ldr++)
    {
      int label = Label(lyr, ldr, "z");
      float val = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0), TMath::TwoPi());
      val *= geo->SensorRadius(lyr, ldr, 0);
      if (label==104)
        val = 0;
      fs << label << " " << val << endl;
    }

  // r*phi sum constraints for W then E arm
  fs << "Constraint 0.0" << endl;
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=0; ldr<geo->GetNLadders(lyr)/2; ldr++)
    {
      int label = Label(lyr, ldr, "s");
      float val = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0), TMath::TwoPi());
      val *= geo->SensorRadius(lyr, ldr, 0);
      if (label==104)
        val = 0;
      fs << label << " " << val << endl;
    }
  fs << "Constraint 0.0" << endl;
  for (int lyr=0; lyr<geo->GetNLayers(); lyr++)
    for (int ldr=geo->GetNLadders(lyr)/2; ldr < geo->GetNLadders(lyr); ldr++)
    {
      int label = Label(lyr, ldr, "s");
      float val = TMath::PiOver2() + fmod(geo->SensorPhiRad(lyr, ldr, 0), TMath::TwoPi());
      val *= geo->SensorRadius(lyr, ldr, 0);
      if (label==104)
        val = 0;
      fs << label << " " << val << endl;
    }

  fs.close();
  Printf("Done.");
  return;
}

void
WriteSteerFile(const char *filename, const char *binfile, const char *constfile)
{
  Printf("Writing %s", filename);

  ofstream fs(filename);

  fs << Form("* This is %s created by %s", filename, "SimVtxAlign.C") << endl;
  fs << Form("* Pass this file to pede: pede %s", filename) << endl;
  fs << endl;

  fs << "Fortranfiles  ! Fortran/text inputs listed here:" << endl;
  fs << Form("%s  ! constraints text file", constfile) << endl;
  fs << endl;

  fs << "Cfiles  ! c/c++ binary input files listed here:" << endl;
  fs << Form("%s  ! binary data file", binfile) << endl;
  fs << endl;

  fs << "method inversion 3 0.001  ! Gauss. elim., #iterations, tol." << endl;
  fs << "end" << endl;

  fs.close();

  Printf("Done.");

  return;
}

int
Label(int layer, int ladder, string coord)
{
  // Return a global parameter label for this layer, ladder, and coordinate.
  // In Millepede II, any unique integer > 0 will do.
  // Labels are not required to be sequential.
  int start[4] = {0,10,30,46}; // # ladders within & excluding layer 0,1,2,3
  int ic = 99;

  if (coord.compare("s") == 0)
    ic = 0;
  if (coord.compare("z") == 0)
    ic = 1;
  // Add other coordinates / degrees of freedom here as needed
  // Then ParInfo() would need corresponding modification.

  return 70*ic + start[layer] + ladder + 1;
}

void
ParInfo(int label, int &layer, int &ladder, string &coord)
{
  // Get layer, ladder, and coordinate string from global parameter label.
  // Inverse of Label() function.

  int start[4] = {0,10,30,46}; // # ladders within & excluding layer 0,1,2,3
  int l = 0;

  if (label > 70)
  {
    coord = "z";
    l = label - 70;
  }
  else
  {
    coord = "s";
    l = label;
  }

  layer = 3;
  for (int i=3; i>=0; i--)
    if (l < start[i] + 1)
      layer = i - 1;

  ladder = l - start[layer] - 1;

  return;
}

TCanvas *
DrawXY(SvxTGeo *geo, const char *name, const char *title, TString opt)
{
  TCanvas *c = new TCanvas(name, title, 900, 900);
  TH1F *hf = c->DrawFrame(-20, -20, 20, 20, title);
  hf->SetXTitle("East                 x [cm]                 West");
  hf->GetXaxis()->CenterTitle();
  hf->SetYTitle("y [cm]");
  hf->GetYaxis()->CenterTitle();

  if (opt.Contains("L"))
  {
    TLatex label;
    label.SetTextSize(0.02);
    double xyz[3] = {0};
    for (int i=0; i<geo->GetNLayers(); i++)
      for (int j=0; j<geo->GetNLadders(i); j++)
      {
        TPolyLine *s = geo->LadderOutlineXY(i,j);
        s->SetLineColor(opt.Contains("faint") ? kGray : kGray+2);
        s->SetLineWidth(opt.Contains("faint") ? 1 : 2);
        s->Draw();

        if (opt.Contains("faint"))
          continue;

        geo->GetSensorXYZ(i, j, 0, xyz);
        int horz = xyz[0] > 0 ? 1 : 3;
        int vert = xyz[1] > 0 ? 1 : 3;
        label.SetTextAlign(10*horz + vert);
        label.DrawLatex(xyz[0], xyz[1], Form(" %d ", j));
      }
  }

  return c;
}

void
GetLadderXYZ(SvxTGeo *tgeo, vecd &x, vecd &y, vecd &z)
{
  int ipt = 0;
  for (int i=0; i<tgeo->GetNLayers(); i++)
    for (int j=0; j<tgeo->GetNLadders(i); j++)
    {
      double xyz[3] = {0};
      tgeo->GetSensorXYZ(i,j,0,xyz);
      x.push_back(xyz[0]);
      y.push_back(xyz[1]);
      z.push_back(xyz[2]);

      if (false)
        Printf("xyz %f %f %f", x[ipt], y[ipt], z[ipt]);

      ipt++;
    }

  return;
}

void
DrawDiffs(vecd &x1, vecd &y1, vecd &z1, vecd &x2, vecd &y2, vecd &z2)
{
  int n = (int)x1.size();

  for (int i=0; i<n; i++)
  {
    double f = 20;
    double x = x1[i], y = y1[i];
    double dx = x2[i] - x;
    double dy = y2[i] - y;
    double dz = z2[i] - z1[i];

    // Draw points showing displacement in z coordinate
    TMarker mkr;
    mkr.SetMarkerStyle(kOpenCircle);
    mkr.SetMarkerColor(kGray+1);
    mkr.SetMarkerSize(0.5);
    //    mkr.DrawMarker(x,y);

    if (TMath::Abs(dz) > 5e-4) // Only label changes > 5 um
    {
      mkr.SetMarkerStyle(kOpenCircle);
      mkr.SetMarkerColor(dz>0 ? kRed-4 : kAzure-4); // Red/blue shift mnemonic
      mkr.SetMarkerSize(100*TMath::Abs(dz)); // 1 = 8px diam, 2 = 16px, ...
      mkr.DrawMarker(x,y);

      TLatex ltx;
      ltx.SetTextSize(0.018);
      ltx.SetTextAlign(22);
      ltx.SetTextColor(dz>0 ? kRed+2 : kAzure+3);
      ltx.DrawLatex(x, y, Form("%.0f", 1e4*dz));
    }

    // Draw arrows showing (significant) displacements in xy plane
    if (dx*dx + dy*dy > 5e-4) // Only label changes > 5 um
    {
      TArrow a;
      a.SetLineWidth(2);
      a.DrawArrow(x, y, x + f*dx, y + f*dy, 0.005);

      TLatex ltx;
      ltx.SetTextSize(0.01);
      ltx.DrawLatex(x + f*dx, y + f*dy, Form("(%.0f, %.0f)", 1e4*dx, 1e4*dy));
    }
  }

  return;
}

int
GetCorrections(const char *resFile, std::map<int, double> &mpc)
{
  std::ifstream filein(resFile);
  int label;
  double p, col3, col4, col5;

  if (!filein)
  {
    Printf("Error opening %s", resFile);
    return -1;
  }
  else
    for (std::string line; std::getline(filein, line);)
    {
      if (filein >> label >> p >> col3 >> col4 >> col5)
        mpc[label] = p;
    }

  return (int)mpc.size();
}
