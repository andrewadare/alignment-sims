// $Id: VtxAlignmentUtils.h,v 1.12 2014/05/08 00:02:14 adare Exp $

#include "SvxGeoTrack.h"
#include <TNtuple.h>
#include <TDecompSVD.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMath.h>
#include <TGeoNode.h>
#include <TGeoMatrix.h>

using namespace std;

TVectorD SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L);
void FillHitNTuple(SvxGeoTrack &gt, TNtuple *ntuple);
void TrackFitZResid(SvxGeoTrack &gt, double *pars = 0);
void TrackFitSResid(SvxGeoTrack &gt, double *pars = 0);
void ZeroFieldResiduals(SvxGeoTrack &gt, double *pars /* y0, z0, phi, theta */);
void Residuals(SvxGeoTrack &tt, SvxGeoTrack &mt, TNtuple *t);


TVectorD
SolveGLS(TMatrixD &X, TVectorD &y, TMatrixD &L)
{
  // Simple generalized linear least-squares fitter.
  // Solve y(X) = beta'*X + error(L) for beta (vector of regression coefs.)
  // L is the inverse covariance matrix for y.
  // Least-squares solution involves solving X' L X beta = X' L y

  TMatrixD XT(X); XT.T();
  TMatrixD A = XT*L*X;
  TVectorD b = XT*L*y;

  // Now solve A*beta = b using SVD. Decompose A = U S V'.
  TDecompSVD svd(A);
  TMatrixD UT = svd.GetU(); UT.T();

  // Construct Moore-Penrose pseudoinverse of S
  TVectorD s = svd.GetSig();
  TMatrixD Sd(s.GetNrows(), s.GetNrows());
  for (int i=0; i<s.GetNrows(); i++)
    Sd(i,i) = s(i)>0 ? 1./s(i) : 0.;

  TVectorD beta = svd.GetV() * Sd * UT * b;

  return beta;
}

void
ZeroFieldResiduals(SvxGeoTrack &gt, double *pars /* y0, z0, yslope, zslope */)
{
  if (gt.nhits < 1)
  {
    Printf("ZeroFieldResiduals(): No hits in track. Skipping.");
    return;
  }

  double zpars[2] = {0}, spars[2] = {0};
  TrackFitZResid(gt, zpars);
  TrackFitSResid(gt, spars);

  pars[0] = spars[0]; // y0 (y-intercept)
  pars[1] = zpars[0]; // z0 (z-intercept)
  pars[2] = spars[1]; // slope in transverse (y vs x) plane
  pars[3] = zpars[1]; // slope in longitudinal (z vs r) plane

  return;
}

void
TrackFitZResid(SvxGeoTrack &gt, double *pars)
{
  // Perform straight-line fit z(r) = z0 + c*r.
  // Assign residuals to gt.hits[i].dz. Put [z0, theta] from fit in pars.
  int m = gt.nhits, n = 2;
  TMatrixD X(m, n);
  TMatrixD Cinv(m, m);
  TVectorD y(m);

  for (int ihit=0; ihit < m; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    X(ihit,0) = 1;
    X(ihit,1) = TMath::Sqrt(hit.x*hit.x + hit.y*hit.y);
    Cinv(ihit,ihit) = (hit.zsigma > 0) ? 1./hit.zsigma : 1.;
    y(ihit) = hit.z;
  }

  TVectorD beta = SolveGLS(X, y, Cinv);
  double z0 = beta(0), c = beta(1);

  for (int ihit=0; ihit<m; ihit++)
  {
    SvxGeoHit hit = gt.GetHit(ihit);
    double zproj = z0 + c*TMath::Sqrt(hit.x*hit.x + hit.y*hit.y);

    // double gxyz[3] = {hit.x, hit.y, zproj};
    // double lxyz[3] = {0};
    // assert(hit.node);
    // hit.node->GetMatrix()->MasterToLocal(gxyz, lxyz);

    // gt.hits[ihit].dz = hit.zs - lxyz[2];
    //    gt.hits[ihit].dz = lxyz[2] - hit.zs;
    gt.hits[ihit].dz = zproj - hit.z;
  }

  if (pars)
  {
    pars[0] = z0;
    pars[1] = TMath::ATan2(1.0, c);
  }

  return;
}

void
TrackFitSResid(SvxGeoTrack &gt, double *pars)
{
  // Perform straight-line fit y' = m'*x' + b' after rotating points
  // to approximately lie along the x axis. Then by construction,
  // x' ~ r, m' ~ 0, and the error is (mostly) in the y' = phi direction.
  // Assign s = r*phi residual to gt.hits[i].ds.
  // Optionally put [y0, phi] from fit in pars array.

  int m = gt.nhits, n = 2;
  TMatrixD points(n, m); // Columns are x',y' pairs
  TMatrixD X(m, n);      // Column 0 is 1's, column 1 is x' (~r) hit coords.
  TMatrixD Cinv(m, m);   // Inverse covariance matrix. Currently diagonal.
  TVectorD y(m);         // Dependent variables y'.
  double phirot = 0;

  for (int ihit=0; ihit<m; ihit++)
  {
    SvxGeoHit t = gt.GetHit(ihit);
    points(0,ihit) = t.x;
    points(1,ihit) = t.y;
    X(ihit,0) = 1;
    Cinv(ihit,ihit) = (t.xsigma > 0) ? 1./t.xsigma : 1.;
    phirot += 1./m * TMath::ATan2(t.y, t.x);
  }

  // Rotate x, y by -phirot so error is approximately in y' direction only.
  // In rotated frame (xp, yp) = (cx + sy, cy - sx)
  TMatrixD R(2,2);
  double c = TMath::Cos(phirot), s = TMath::Sin(phirot);
  R(0,0) =  c; R(0,1) = s;
  R(1,0) = -s; R(1,1) = c;

  points = R * points;

  TMatrixDColumn(X, 1) = TMatrixDRow(points, 0);
  y = TMatrixDRow(points, 1);

  // Fit track to get [b', m'] in rotated system.
  TVectorD beta = SolveGLS(X,y,Cinv);
  double bp = beta(0), mp = beta(1);

  // Rotate back to get y-intercept b of track in x=0 plane
  // and phi angle of track.
  double y0  = bp / (c - mp*s);
  double phi = phirot + TMath::ATan(mp);
  if (phi < 0)
    phi += TMath::TwoPi();
  if (pars)
  {
    pars[0] = y0;
    pars[1] = phi;
  }

  for (int ihit=0; ihit<m; ihit++)
  {
    double x  = points(0, ihit);
    double ds = mp*x + bp - points(1, ihit); // projected - measured y'
    gt.hits[ihit].ds = (x < 0) ? -ds : +ds;
  }

  /*
    // Compute track projections at hit radii --> ds residuals
    for (int ihit=0; ihit<m; ihit++)
    {
      SvxGeoHit hit  = gt.GetHit(ihit);
      double rho     = TMath::Sqrt(hit.x*hit.x + hit.y*hit.y);
      double xproj   = rho*TMath::Cos(phi);
      double yproj   = rho*TMath::Sin(phi) + y0;
      double gxyz[3] = {xproj, yproj, hit.z};
      double lxyz[3] = {0}; // Projection wrt sensor center (only x is valid)
      hit.node->GetMatrix()->MasterToLocal(gxyz, lxyz);

      double dphi = atan2(lxyz[0] - hit.xs, rho);

      gt.hits[ihit].ds = hit.x > 0 ? rho*dphi : -rho*dphi;
    }
  */

  return;
}

void
FillHitNTuple(SvxGeoTrack &gt, TNtuple *ntuple)
{
  // TNtuple *ntuple = new TNtuple("t", "SvxGeoHit variables",
  // "layer:ladder:sensor:xs:ys:zs:x:y:z:xsigma:zsigma:dz:ds:trkid");

  assert(ntuple);

  if (false)
    Printf("ds (%8.4f,%8.4f,%8.4f,%8.4f),  "
           "dz (%8.4f,%8.4f,%8.4f,%8.4f)",
           gt.hits[0].ds, gt.hits[1].ds, gt.hits[2].ds, gt.hits[3].ds,
           gt.hits[0].dz, gt.hits[1].dz, gt.hits[2].dz, gt.hits[3].dz);

  for (int ihit=0; ihit<gt.nhits; ihit++)
  {
    int nj = 14;
    std::vector<float> vars(nj, 0.);
    int j = 0;
    SvxGeoHit hit = gt.GetHit(ihit);

    vars[j++] = hit.layer  ;
    vars[j++] = hit.ladder ;
    vars[j++] = hit.sensor ;
    vars[j++] = hit.xs     ;
    vars[j++] = hit.ys     ;
    vars[j++] = hit.zs     ;
    vars[j++] = hit.x      ;
    vars[j++] = hit.y      ;
    vars[j++] = hit.z      ;
    vars[j++] = hit.xsigma ;
    vars[j++] = hit.zsigma ;
    vars[j++] = hit.dz     ;
    vars[j++] = hit.ds     ;
    vars[j++] = hit.trkid  ;
    ntuple->Fill(&vars[0]);
  }

  return;
}
