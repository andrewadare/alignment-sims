#include "UtilFns.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TNtuple.h"

TH1D *Resid(int lyr, int ldr, int stage, TNtuple *t, const char *var, const double xmax);

void DrawFakeVtxAlignResults()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  int nLadders[4] = {10,20,16,24}; // ladders/layer

  TFile *inFile = new TFile("rootfiles/vtx-fake.0.root", "read");
  TObjArray *cList = new TObjArray();
  TNtuple *ht[3] = {0};
  ht[0] = (TNtuple *) inFile->Get("ht0");
  ht[1] = (TNtuple *) inFile->Get("ht1");
  ht[2] = (TNtuple *) inFile->Get("ht2");

  TLatex ltx;
  ltx.SetNDC();
  TH1D *hs[4][24][3]; // Layer, ladder, {0=ideal, 1=misaligned, 2=corrected}
  TH1D *hz[4][24][3];
  TH1D *hds[4];
  TH1D *hdz[4];
  for (int lyr=0; lyr<4; lyr++)
    for (int ldr=0; ldr<nLadders[lyr]; ldr++)
      for (int stage=0; stage<3; stage++)
      {
        hs[lyr][ldr][stage] = Resid(lyr, ldr, stage, ht[stage], "ds", 0.12);
        hz[lyr][ldr][stage] = Resid(lyr, ldr, stage, ht[stage], "dz", 0.12);
      }

  // Draw individual residual distributions on sub-pads
  for (int lyr=0; lyr<4; lyr++)
  {
    int nl = nLadders[lyr];
    int nx = lyr ? nl/4 : nl/2; // Number of pads along x
    int ny = lyr ? 4 : 2;       // Number of pads along y
    int ph = 200;               // pad height in px
    int pw = 200;               // pad width in px
    TCanvas *cs = new TCanvas(Form("cs%d", lyr), Form("ds layer %d", lyr),
                              pw*nx, ph*ny);
    TCanvas *cz = new TCanvas(Form("cz%d", lyr), Form("dz layer %d", lyr),
                              pw*nx, ph*ny);
    cs->Divide(nx, ny, 0.001, 0.001);
    cz->Divide(nx, ny, 0.001, 0.001);

    for (int ldr=0; ldr<nl; ldr++)
    {
      cs->cd(ldr+1);
      hs[lyr][ldr][0]->Draw("");
      hs[lyr][ldr][1]->Draw("same");
      hs[lyr][ldr][2]->Draw("same");
      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.2);
      gPad->SetRightMargin(0.01);
      gPad->SetTopMargin(0.01);
      ltx.SetTextSize(0.07);
      ltx.DrawLatex(0.23, 0.92, hs[lyr][ldr][1]->GetTitle());
      ltx.DrawLatex(0.23, 0.85, "ds [cm]");
      ltx.SetTextSize(0.06);
      ltx.DrawLatex(0.62,0.93, Form("Mean %.0f #mum", 1e4*hs[lyr][ldr][1]->GetMean()));
      ltx.DrawLatex(0.62,0.87, Form("RMS  %.0f #mum", 1e4*hs[lyr][ldr][1]->GetRMS()));


      cz->cd(ldr+1);
      hz[lyr][ldr][0]->Draw("");
      hz[lyr][ldr][1]->Draw("same");
      hz[lyr][ldr][2]->Draw("same");
      gPad->SetBottomMargin(0.15);
      gPad->SetLeftMargin(0.2);
      gPad->SetRightMargin(0.01);
      gPad->SetTopMargin(0.01);
      ltx.SetTextSize(0.07);
      ltx.DrawLatex(0.23, 0.92, hz[lyr][ldr][1]->GetTitle());
      ltx.DrawLatex(0.23, 0.85, "dz [cm]");
      ltx.SetTextSize(0.06);
      ltx.DrawLatex(0.62,0.93, Form("Mean %.0f #mum", 1e4*hz[lyr][ldr][1]->GetMean()));
      ltx.DrawLatex(0.62,0.87, Form("RMS  %.0f #mum", 1e4*hz[lyr][ldr][1]->GetRMS()));

      if (false)
        DrawObject(hz[lyr][ldr][1], "", Form("resid_lyr%d_ldr%d",lyr,ldr), cList);
    }
    cList->Add((TCanvas *)cs);
    cList->Add((TCanvas *)cz);
  }


  for (int lyr=0; lyr<4; lyr++)
  {
    hdz[lyr] = new TH1D(Form("hdz%d",lyr), Form("z-resid layer %d", lyr), nLadders[lyr], 0, float(nLadders[lyr]));
    hds[lyr] = new TH1D(Form("hds%d",lyr), Form("s-resid layer %d", lyr), nLadders[lyr], 0, float(nLadders[lyr]));
    SetHistProps(hds[lyr], kBlack, kNone, kRed, kFullCircle, 1.5);
    SetHistProps(hdz[lyr], kBlack, kNone, kBlue, kFullCircle, 1.5);
    hds[lyr]->SetXTitle("Ladder");
    hdz[lyr]->SetXTitle("Ladder");

    hds[lyr]->SetYTitle("#Deltas [cm]");
    hdz[lyr]->SetYTitle("#Deltaz [cm]");

    hds[lyr]->GetYaxis()->SetRangeUser(-0.09, 0.09);
    hdz[lyr]->GetYaxis()->SetRangeUser(-0.09, 0.09);

    hds[lyr]->GetXaxis()->SetNdivisions(210);
    hdz[lyr]->GetXaxis()->SetNdivisions(210);

    for (int ldr=0; ldr<nLadders[lyr]; ldr++)
    {
      hds[lyr]->SetBinContent(ldr+1, hs[lyr][ldr][1]->GetMean());
      hdz[lyr]->SetBinContent(ldr+1, hz[lyr][ldr][1]->GetMean());
      hds[lyr]->SetBinError(ldr+1, hs[lyr][ldr][1]->GetRMS());
      hdz[lyr]->SetBinError(ldr+1, hz[lyr][ldr][1]->GetRMS());
    }

    ltx.SetTextSize(0.07);
    DrawObject(hds[lyr], "e0p", Form("ds_lyr%d",lyr), cList);
    ltx.DrawLatex(0.23, 0.92, hds[lyr]->GetTitle());
    DrawObject(hdz[lyr], "e0p", Form("dz_lyr%d",lyr), cList);
    ltx.DrawLatex(0.23, 0.92, hdz[lyr]->GetTitle());
  }

  ltx.SetTextSize(0.06);
  ltx.SetTextFont(42);

  const char *xyplots[4] = {"vtx_xy","misalign","millepede_dp","align_error"};
  for (int i=0; i<4; i++)
  {
    TCanvas *c = (TCanvas *) inFile->Get(xyplots[i]);
    c->Draw();
    ltx.DrawLatex(0.15, 0.92, gPad->GetTitle());
    cList->Add(c);
  }

  PrintPDFs(cList, "pdfs", "");
  PrintPDF(cList, "pdfs/fake-align");
}

TH1D *
Resid(int lyr, int ldr, int stage, TNtuple *t, const char *var, const double xmax)
{
  TH1D *h = new TH1D(Form("%s%d%d%d", var, lyr, ldr, stage),
                     Form("B%dL%d", lyr, ldr),
                     100, -xmax, xmax);
  //  h->GetXaxis()->SetTitle(Form("%s [cm]"))
  t->Draw(Form("%s >> %s", var, h->GetName()),
          Form("layer==%d && ladder==%d", lyr, ldr),
          "goff");
  if (strcmp(var, "ds")==0)
  {
    //   h->SetFillColor(kGreen-7);
    h->SetLineColor(kGreen+3);
    if (stage==0)
    {
      h->SetLineColor(kGray+1);
      h->SetFillColor(kGray);
    }
    if (stage==1)
    {
      h->SetLineStyle(kDashed);
    }
    if (stage==2)
    {
      ;//      h->SetLineStyle(kDashed);
    }

  }
  if (strcmp(var, "dz")==0)
  {
    //    h->SetFillColor(kAzure-4);
    h->SetLineColor(kAzure+2);
    if (stage==0)
    {
      h->SetLineColor(kGray+1);
      h->SetFillColor(kGray);
    }
    if (stage==1)
    {
      h->SetLineStyle(kDashed);
    }
    if (stage==2)
    {
      ;//      h->SetLineStyle(kDashed);
    }
  }

  if (true)
  {
    TAxis *ax = h->GetXaxis();
    TAxis *ay = h->GetYaxis();
    ay->SetRangeUser(0, 1.3*h->GetMaximum());
    ax->SetNdivisions(205);
    ax->SetLabelSize(0.06);
    ay->SetNdivisions(205);
    ay->SetLabelSize(0.05);
    ax->SetTitleSize(0);
  }
  return h;
}
