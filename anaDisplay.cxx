

// ROOTANA header files
#include "TRootanaDisplay.hxx"
#include "TV792Data.hxx"
#include "TV1190Data.hxx"

// Conditional ROOTANA header files
//#undef USE_V792
//#define USE_V1190
//#define USE_L2249
//#define USE_AGILENT

#ifdef  USE_V792
#include "TV792Histogram.h"
#endif 
#ifdef  USE_V1190
#include "TV1190Histogram.h"
#endif 
#ifdef  USE_L2249
#include "TL2249Histogram.h"
#endif 
#ifdef  USE_AGILENT
#include "TAgilentHistogram.h"
#endif 

// ROOT header files
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TPolyLine3D.h"
#include "TMath.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TFancyHistogramCanvas.hxx"

// C++ header files
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <iomanip>
#include <fstream>
#include <math.h>

const int Nchan = 14;
const int Ncounts = 5;
const Float_t Tcal = 5.0;
int maxcount = 0;
const int Nray = 10;
int nevent = 0;
const bool docal=1;//turn calibration on
Bool_t dogate=1;//turn gating on
const bool free_range=false;//set histogram to auto-range
const bool match_sim=true;//match histogram ranges and binning to simulation
Bool_t bdiag=kTRUE;
Bool_t set_showtest=kFALSE;
//Bool_t set_showoutput=kFALSE;
const Int_t nevents[84]={15460,0,0,0,4570762,8407363,2765657,0,3769,10783,8011,8181,0,13614,0,0,1748,20100,31849,20702,3486088,20412,12227,14655,236743,10033212,7141975,6967745,0,8441,5733,5625,37194,6737045,1587827,20628,0,0,0,0,2375383,32159,16836,15420,3983463,909,0,24941,19733,38811,4916,4911678,8095,39730,28987,27590,5175814,38737,35818,25818,499468,3072779,0,11207650,4980179,828280,33211,0,22758,0,22083,11334,0,5817,40467,31042,20511,8141,12470,29756,2679,29251,2073,28465};

class MyTestLoop: public TRootanaDisplay {

#ifdef  USE_V792
  TV792Histograms *v792_histos;
#endif 
#ifdef  USE_V1190
  TV1190Histograms *v1190_histos;
#endif 

  

  // *********************************************************
  // GLOBAL VARIABLES
  // *********************************************************

  // Define histograms to be used
  
  // A1 - A2
  // [0] = A1B - A2B
  // [1] = A1B - A2M
  // [2] = A1B - A2T
  // [3] = A1M - A2B
  // [4] = A1M - A2M
  // [5] = A1M - A2T
  // [6] = A1T - A2B
  // [7] = A1T - A2M
  // [8] = A1T - A2T
  // [9] = A1 - A2 (All Straight)
  // [10] = A1 - A2 (All Crossing)
  // [11] = A1B - A2B (Copy for plotting purpose)
  // [12] = A1M - A2M (Copy for plotting purpose)
  // [13] = A1T - A2T (Copy for plotting purpose)
  TH1F *a1_a2_diff[14];
  // The X and Y differences
  // [0] = X1R - X1L
  // [1] = Y1B - Y1T
  // [2] = X2R - X2L
  // [3] = Y2B - Y2T
  TH1F *x_y_diff[4];
  // The X and Y sums
  // [0] = X1R + X1L - 1AM
  // [1] = Y1B + Y1T - 1AM
  // [2] = X2R + X2L - 2AM
  // [3] = Y2B + Y2T - 2AM
  TH1F *x_y_sum[4];

  // Debugging histogram
  TH1F *m_counts[12];
  TH1F *anode_pathway_multiplicity;

  // T2 histograms for diff vs sum
  // TH2F *x_y_diff_vs_sum[4];
  //TH2F *x_y_diff_vs_sum_large[4];
  TH2F *x_y_diff_vs_sum_silly[4];

  int ct;
  //  ofstream raw_output;

public:

  //Calibration constants
  Float_t hxxcal[Nchan/2][2];//linear calibration constants
  Float_t hxcal[4][2];
  Float_t hasumcal[4][3];
  Float_t hsumzcal[Nchan/2][3];
  Float_t hdiffzcal[Nchan/2][3];

  MyTestLoop() {
    SetOutputFilename("output_display/emma_anaDisplay_");
    DisableRootOutput(false);

#ifdef  USE_V792
    v792_histos = new TV792Histograms();
#endif 
#ifdef  USE_V1190
    v1190_histos = new TV1190Histograms();
#endif 
  }

  void AddAllCanvases() {

    SetNumberSkipEvent(100);
    // Set up tabbed canvases
#ifdef  USE_V792
    TFancyHistogramCanvas* v792_all = new TFancyHistogramCanvas(new TV792Histograms(),"V792");
    AddSingleCanvas(v792_all);
#endif 
#ifdef  USE_V1190
    TFancyHistogramCanvas* v1190_all = new TFancyHistogramCanvas(new TV1190Histograms(),"V1190");
    AddSingleCanvas(v1190_all);
#endif 
#ifdef  USE_L2249
    TFancyHistogramCanvas* l2249_all = new TFancyHistogramCanvas(new TL2249Histograms(),"L2249");
    AddSingleCanvas(l2249_all);
#endif 
#ifdef  USE_AGILENT
    TFancyHistogramCanvas* agilent_all = new TFancyHistogramCanvas(new TAgilentHistograms(),"AGILENT");
    AddSingleCanvas(agilent_all);
#endif 

    // ****************************************
    // ADD CANVAS
    // ****************************************

    AddSingleCanvas("A1-A2 Summary"); 
    AddSingleCanvas("A1-A2 Detail"); 
    AddSingleCanvas("M_Counts"); 

    AddSingleCanvas("XY Diffs"); 
    AddSingleCanvas("XY Sums"); 

    AddSingleCanvas("XY Diff vs Sum (Full)");
   
    //    AddSingleCanvas("XY Positions"); 
    SetDisplayName("Online Display");
    gStyle->SetOptStat("neMi"); //sets what info is disp. in histograms
    gStyle->SetOptFit(0111); //sets info disp. in fits
  };

  virtual ~MyTestLoop() {};

  int v1742_newrun;

  void BeginRun(int transition,int run,int time) {

    v1742_newrun = 1;

#ifdef  USE_V792
    v792_histos->BeginRun(transition,run,time);
#endif 
#ifdef  USE_V1190
    v1190_histos->BeginRun(transition,run,time);
#endif 
  
    // **************************************
    // INITIALIZATION 
    // **************************************
    ct = 0;
    // initialize histograms
    char name[100];
    Float_t xdiffmin=-90;
    Float_t ydiffmin=-150;
    Float_t xdiffrange=470;
    Float_t ydiffrange=310;
    Float_t xsummin=310;
    Float_t ysummin=90;
    Float_t xsumrange=80;
    Float_t ysumrange=xsumrange;
    Float_t diffbin=512;
    Float_t sumbin=100;

    // initialize anode_pathway_multiplicity histogram
    anode_pathway_multiplicity = new TH1F("APM","Anode Pathway Multiplicity Index",511,0.5,511.5);
    anode_pathway_multiplicity->SetXTitle("Index");

    // initialize a1_a2_diff histogram
    char *title[] = {"A1B - A2B", "A1B - A2M", "A1B - A2T",
                     "A1M - A2B", "A1M - A2M", "A1M - A2T",
                     "A1T - A2B", "A1T - A2M", "A1T - A2T",
                     "A1 - A2 (Straight)", "A1 - A2 (Crossing)",
                     "A1B - A2B", "A1M - A2M", "A1T - A2T"};
    for(int i = 0; i < 14; i++){   
      sprintf(name,"a1_a2_diff_%i",i);
      a1_a2_diff[i] = new TH1F(name,title[i],150,-14.5,15.5);
      a1_a2_diff[i]->SetXTitle("A1-A2 time diff (ns)");
    }
    a1_a2_diff[0]->SetLineColor(2);
    a1_a2_diff[4]->SetLineColor(1);
    a1_a2_diff[8]->SetLineColor(4);
    a1_a2_diff[11]->SetLineColor(2);
    a1_a2_diff[12]->SetLineColor(1);
    a1_a2_diff[13]->SetLineColor(4);

    // initialize m_count histogram
    title = {"Anode Multiplicity", "A2M", "X1L", "X1R", "Y1B", "Y1T", "X2L", "X2R", "Y2B", "Y2T", "Anode & Cathode Measurement Total", "64-Ch Measurement Total"};
    for(int i = 0; i < 12; i++){
      sprintf(name,"m_count_%i",i);
      if(i<10)
        m_counts[i] = new TH1F(name,title[i],5,0,4);
      else
        m_counts[i] = new TH1F(name,title[i],100,0,100);
    }

    // initialize x_y_diff histogram
    title = {"X1R - X1L", "Y1B - Y1T", "X2R - X2L", "Y2B - Y2T"};
    for(int i = 0; i < 4; i++){
      sprintf(name,"x_y_diff_%i",i);
      if(i == 0)
	x_y_diff[i] = new TH1F(name,title[i],diffbin,xdiffmin,xdiffmin+xdiffrange);
      if(i == 1)
	x_y_diff[i] = new TH1F(name,title[i],diffbin,ydiffmin,ydiffmin+ydiffrange);
      if(i==2)
	x_y_diff[i] = new TH1F(name,title[i],diffbin,xdiffmin,xdiffmin+xdiffrange);     
      if(i==3)
	x_y_diff[i] = new TH1F(name,title[i],diffbin,ydiffmin,ydiffmin+ydiffrange);	
      x_y_diff[i]->SetXTitle("Time Diff (ns)");
    }

    // initialize x_y_sum histogram
    title = {"X1R + X1L", "Y1T + Y1B", "X2R + X2L", "Y2T + Y2B"};
    for(int i = 0; i < 4; i++){
      sprintf(name,"x_y_sum_%i",i);
      if(i == 0)
	x_y_sum[i] = new TH1F(name,title[i],sumbin,xsummin,xsummin+xsumrange);
      if(i == 1)
	x_y_sum[i] = new TH1F(name,title[i],sumbin,ysummin,ysummin+ysumrange);
      if(i==2)
	x_y_sum[i] = new TH1F(name,title[i],sumbin,xsummin,xsummin+xsumrange);
      if(i==3)
	x_y_sum[i] = new TH1F(name,title[i],sumbin,ysummin,ysummin+ysumrange);
      x_y_sum[i]->SetXTitle("Time Sum (ns)");
    }

    // initialize x_y_diff_vs_sum histogram
    title = {"X1R/X1L Diff vs Sum", "Y1B/Y1T Diff vs Sum", "X2R/X2L Diff vs Sum", "Y2B/Y2T Diff vs Sum"};
    for(int i = 0; i < 4; i++){
      sprintf(name,"x_y_diff_vs_sum_silly_%i",i);
      if(i == 0)
	x_y_diff_vs_sum_silly[i] = new TH2F(name,title[i],sumbin,xsummin,xsummin+xsumrange,diffbin,xdiffmin,xdiffmin+xdiffrange);
      if(i == 1)
	x_y_diff_vs_sum_silly[i] = new TH2F(name,title[i],sumbin,ysummin,ysummin+ysumrange,diffbin,ydiffmin,ydiffmin+ydiffrange);
      if(i==2)
	x_y_diff_vs_sum_silly[i] = new TH2F(name,title[i],sumbin,xsummin,xsummin+xsumrange,diffbin,xdiffmin,xdiffmin+xdiffrange);
      if(i==3)
	x_y_diff_vs_sum_silly[i] = new TH2F(name,title[i],sumbin,ysummin,ysummin+ysumrange,diffbin,ydiffmin,ydiffmin+ydiffrange);
      x_y_diff_vs_sum_silly[i]->SetXTitle("Cathode Time Sum (ns)");      
      x_y_diff_vs_sum_silly[i]->SetYTitle("Cathode Time Differecnce (ns)");
    }
  }

  void EndRun(int transition,int run,int time){
  }

  void ResetHistograms() {
    // ***********************************
    // RESET HISTOGRAMS
    // ***********************************
    for(int i =0; i < 14; i++)
      a1_a2_diff[i]->Reset();
    for(int i =0; i < 12; i++)
      m_counts[i]->Reset();
    for(int i =0; i < 4; i++) {
      x_y_diff[i]->Reset();
      x_y_sum[i]->Reset();
      x_y_diff_vs_sum_silly[i]->Reset();
    }
  }

  int adc_vmax[32];
  int adc_tmax[32];
  int adc_tthr[32];

  int adc_ttri[4];

  TH2D*hhxx1;
  TH2D*hhxx2;

  TH2D*hhxx1aa;
  TH2D*hhxx1tt;

  TH2D*hhxx2aa;
  TH2D*hhxx2tt;

  double t6;
  double t7;

  void doV1742(TDataContainer& dataContainer)
  {
    memset(adc_vmax, 0, sizeof(adc_vmax));
    memset(adc_tmax, 0, sizeof(adc_tmax));
    memset(adc_tthr, 0, sizeof(adc_tthr));

    const TMidasEvent& e = dataContainer.GetMidasEvent();

    void *ptr;
    int bklen,bktype;
    int status = e.FindBank("VADC", &bklen, &bktype, &ptr);

    /// If we couldn't find bank, return null.
    if (status != 1)
      return;

    uint32_t *data32 = (uint32_t*)ptr;
    //printf("V1742 pointer: %p, status %d, len %d, type %d\n", data32, status, bklen, bktype);

    int verbose_unpack = true;

    static int once = 1;
    if (once) {
      once = 0;
    } else {
      verbose_unpack = false;
    }

    // print header

    if (verbose_unpack) {
      printf("Header:\n");
      printf("  total event size: 0x%08x (%d)\n", data32[0], data32[0]&0x0FFFFFFF);
      printf("  board id, pattern, gr mask: 0x%08x\n", data32[1]);
      printf("  event counter: 0x%08x (%d)\n", data32[2], data32[2]);
      printf("  event time tag: 0x%08x (%d)\n", data32[3], data32[3]);
    }

    int grpMask = data32[1] & 0xF;

    int adc[32][1024];
    int adc_tr[4][1024];

    memset(adc, 0, sizeof(adc));
    memset(adc_tr, 0, sizeof(adc_tr));

    uint32_t *g = data32 + 4;
    for (int i=0; i<4; i++) {

      if (((1<<i)&grpMask)==0)
	continue;

      int len = g[0] & 0xfff;
      int tr  = (g[0]>>12)&1;
      int freq = (g[0]>>16)&3;
      int cell = (g[0]>>20)&0x3ff;

      if (verbose_unpack) {
	printf("Group %d:\n", i);
	printf("  group description: 0x%08x, cell %4d, freq %d, tr %d, size %5d\n", g[0], cell, freq, tr, len);
      }

      g += 1;
      // g points to the data
      
      //for (int k=0; k<10; k++)
      //	printf("  adc data[k]: 0x%08x\n", g[k]);

      int k=0;

      const uint8_t* p = (const uint8_t*)g;
      int x = 0;
      for (int s=0; s<1024; s++)
	for (int a=0; a<8; a++) {
	  int v = 0;
	  if (x==0) {
	    v = (p[0]) | ((p[1]&0xF)<<8);
	    p += 1;
	    x = 1;
	  } else {
	    v = ((p[0]&0xF0)>>4) | ((p[1]&0xFF)<<4);
	    p += 2;
	    x = 0;
	  }

	  cell = 0;
	  int xs = (s + cell)%1024;

	  //printf("group %d, channel %d, sample %d: value %6d (0x%03x)\n", i, a, s, v, v);

	  adc[i*8+a][xs] = v;

	  k++;
	  //if (k > 10)
	  //abort();
	}

      g += len;

      if (tr) {
	int trlen = len/8;

	const uint8_t* p = (const uint8_t*)g;
	int x = 0;
	for (int s=0; s<1024; s++) {
	  int v = 0;
	  if (x==0) {
	    v = (p[0]) | ((p[1]&0xF)<<8);
	    p += 1;
	    x = 1;
	  } else {
	    v = ((p[0]&0xF0)>>4) | ((p[1]&0xFF)<<4);
	    p += 2;
	    x = 0;
	  }
	    
	  cell = 0;
	  int xs = (s + cell)%1024;
	    
	  //printf("group %d, channel %d, sample %d: value %6d (0x%03x)\n", i, a, s, v, v);
	    
	  adc_tr[i][xs] = v;
	    
	  k++;
	  //if (k > 10)
	  //abort();
	}

	g += trlen;
      }

      // g points to the time tag
      if (verbose_unpack) {
	printf("  group trigger time tag: 0x%08x\n", g[0]);
      }
      g += 1;
      // g point s to the next group
    }

#if 0
    static TCanvas* gX = NULL;
    if (!gX) {
      gX = new TCanvas("V1742 ADC");
      gX->cd();
      gX->Divide(1,10);
    }
#endif

    static TH1D*hhxt[32];
    static TH1D*hhxv[32];

#define MEMZERO(a) (memset((a), 0, sizeof(a)))

    if (v1742_newrun) {
      MEMZERO(hhxt);
      MEMZERO(hhxv);
    }

    //static int hxmap[8] = { 0, 3, 2, 5, 6, 18, 7, 19 };

    for (int c=0; c<32; c++) {

      if (!hhxv[c]) {
	char buf[256];
	sprintf(buf, "hhxv%d", c);
	hhxv[c] = new TH1D(buf, buf, 100, 0, 4200);
      }

      if (!hhxt[c]) {
	char buf[256];
	sprintf(buf, "hhxt%d, adc < 3550", c);
	hhxt[c] = new TH1D(buf, buf, 1024, 0, 1024);
      }

      int tmax = 0;
      int vmax = adc[c][tmax];
      for (int s=0; s<1024; s++) {
	if (adc[c][s] < vmax) {
	  tmax = s;
	  vmax = adc[c][s];
	}
      }

      adc_vmax[c] = vmax;
      adc_tmax[c] = tmax;

      hhxv[c]->Fill(vmax);
      if (vmax < 3550)
	hhxt[c]->Fill(tmax);

      int ped = adc[c][0];
      int thr = (ped+vmax)/2;
      int tthr = 0;
      for (int s=0; s<1024; s++) {
	if (adc[c][s] < thr) {
	  tthr = s;
	  break;
	}
      }

      //if (tthr > tmax)
      //tthr = 0;

      adc_tthr[c] = tthr;
    }


    static TH1D*htri[4];

    if (v1742_newrun)
      MEMZERO(htri);

    for (int c=0; c<4; c++) {

      if (!htri[c]) {
	char buf[256];
	sprintf(buf, "htri%d", c);
	htri[c] = new TH1D(buf, buf, 1024, 0, 1024);
      }

      int tmax = 0;
      int vmax = adc_tr[c][tmax];
      for (int s=0; s<1024; s++) {
	if (adc_tr[c][s] < vmax) {
	  tmax = s;
	  vmax = adc_tr[c][s];
	}
      }

      int ped = adc_tr[c][0];
      int thr = (ped+vmax)/2;
      int tthr = 0;
      for (int s=0; s<1024; s++) {
	if (adc_tr[c][s] < thr) {
	  tthr = s;
	  break;
	}
      }

      adc_ttri[c] = tthr;

      htri[c]->Fill(tthr);
    }

    static TH1D*hhxa12a;
    static TH1D*hhxa12b;
    static TH1D*hhxa12c;

    if (v1742_newrun) {
      hhxa12a = NULL;
      hhxa12b = NULL;
      hhxa12c = NULL;
      hhxx1 = NULL;
      hhxx2 = NULL;
      hhxx1aa = NULL;
      hhxx1tt = NULL;
      hhxx2aa = NULL;
      hhxx2tt = NULL;
    }

    //static int hxmap[8] = { 0, 3, 2, 5, 6, 18, 7, 19 };

    if (!hhxa12a)
      hhxa12a = new TH1D("A1-A2 section A", "A1-A2 section A, ADC time bins", 101, -50, 50);

    if (!hhxa12b)
      hhxa12b = new TH1D("A1-A2 section B", "A1-A2 section B, ADC time bins", 101, -50, 50);

    if (!hhxa12c)
      hhxa12c = new TH1D("A1-A2 section C", "A1-A2 section C, ADC time bins", 101, -50, 50);

    int thr = 3550;
    
    if (adc_vmax[0]<thr && adc_vmax[3] < thr)
      hhxa12a->Fill(adc_tthr[0] - adc_tthr[3]);

    if (adc_vmax[1]<thr && adc_vmax[4] < thr)
      hhxa12b->Fill(adc_tthr[1] - adc_tthr[4]);

    if (adc_vmax[2]<thr && adc_vmax[5] < thr)
      hhxa12c->Fill(adc_tthr[2] - adc_tthr[5]);

    if (!hhxx1aa)
      hhxx1aa = new TH2D("X1 ADC vs ADC", "X1 ADC vs ADC (time bin)", 100, 0, 1024, 100, 0, 1024);

    if (!hhxx1tt)
      hhxx1tt = new TH2D("X1 TDC vs TDC", "X1 TDC vs TDC (ns)", 100, -100, 500, 100, -100, 500);

    if (!hhxx1)
      hhxx1 = new TH2D("X1 ADC vs TDC", "X1 ADC (time bin) vs TDC (ns)", 100, 0, 1024, 100, -100, 500);

    if (!hhxx2)
      hhxx2 = new TH2D("X2 ADC vs TDC", "X2 ADC (time bin) vs TDC (ns)", 100, 0, 1024, 100, -100, 500);

    static TH1D* hhxx1a;
    static TH1D* hhxx2a;

    if (v1742_newrun) {
      hhxx1a = NULL;
      hhxx2a = NULL;
    }

    if (!hhxx1a)
      hhxx1a = new TH1D("X1 ADC - TDC", "X1 ADC - TDC (ns)", 101, -100, 100);

    if (!hhxx2a)
      hhxx2a = new TH1D("X2 ADC - TDC", "X2 ADC - TDC (ns)", 101, -100, 100);

    static TH1D* hhxxsum;

    if (v1742_newrun) {
      hhxxsum = NULL;
    }

    if (!hhxxsum)
      hhxxsum = new TH1D("X ADC SUM", "X ADC SUM (ns)", 200, 200, 400);

    int x1a = 8;
    int x1b = 9;

    int a6 = adc_tthr[x1a];
    int a7 = adc_tthr[x1b];

    double abin = 1.0; // ns per adc time bin
    double a6a = abin*(a6 - adc_ttri[1]);
    double a7a = abin*(a7 - adc_ttri[1]);

    double asum = a6a + a7a;

    if ((adc_vmax[x1a] < thr) && (adc_vmax[x1b] < thr)) {
      hhxx1tt->Fill(t6, t7);
      hhxx1aa->Fill(a6, a7);
    }

    if (adc_vmax[x1a] < thr) {
      hhxx1->Fill(a6, t6);
      hhxx1a->Fill(a6a - t6);
    }

    if (adc_vmax[x1b] < thr) {
      hhxx2->Fill(a7, t7);
      hhxx2a->Fill(a7a - t7);
    }

    hhxxsum->Fill(asum);

    /*    printf("TDC6 %6.1f, TDC7 %6.1f, TSUM %6.1f, ATRI %4d, ADC%d %4d %4d %4d %6.1f, ADC%d %4d %4d %4d %6.1f, ASUM %6.1f\n",
	  t6, 
	  t7,
	  t6+t7,
	  adc_ttri[0],
	  x1a, adc_vmax[x1a], adc_tmax[x1a], adc_tthr[x1a], a6a,
	  x1b, adc_vmax[x1b], adc_tmax[x1b], adc_tthr[x1b], a7a,
	  asum);
    */

#if 0
    static int scaledown = 0;
    scaledown++;
    if ((scaledown%10) != 1)
      return;
#endif

    static time_t tx = 0;
    time_t now = time(NULL);
    //printf("%d %d %d\n", tx, now, now-tx);
    if (tx==0 || now>=tx) {
      // continue
      tx = now + 2;
    } else {
      if (!v1742_newrun) // kludge alert!
	return;
    }

    printf("plot v1742 waveforms!\n");

    static TCanvas* gC = NULL;
    if (!gC) {
      //gC = new TCanvas("V1742 ADC","V1742 ADC",2000,1200);//online
      gC = new TCanvas("V1742 ADC","V1742 ADC",1400,700);
      if(!(gC->GetShowEventStatus()))gC->ToggleEventStatus();
      if(!(gC->GetShowToolBar()))gC->ToggleToolBar();
      gStyle->SetOptStat("neMi"); //sets what info is disp. in histograms
      gStyle->SetOptFit(0111); //sets info disp. in fits
      gC->cd();      
      gC->Divide(3,3);
    }

    static TH1D*ha = NULL;
    static TH1D*hb = NULL;
    static TH1D*hc = NULL;
    static TH1D*hd = NULL;

    static TH1D*hh[14];

    if (v1742_newrun) {
      ha = NULL;
      hb = NULL;
      hc = NULL;
      hd = NULL;
      MEMZERO(hh);
    }

    static int hmap[14] = { 0, 3, 1, 4, 2, 5,   8, 9, 10, 11, 12, 13, 14, 15 };
    TString title;    
    TString titles[14]={"A1B","A1M","A1T","A2B","A2M","A2T","X1L","X1R","Y1B","Y1T","X2L","X2R","Y2B","Y2T"};

    for (int h=0; h<14; h++) {

      int a = hmap[h];

      gC->cd((int)(h/2)+1);
      
      if (!hh[h]) {
	char buf[256];
	sprintf(buf, "ch%d", a);
	title=titles[a];
	title+=" V1742 Ch.";
	title+=a;
	hh[h] = new TH1D(buf, title, 1024, 0, 1024);
      }

      hh[h]->Clear();
      hh[h]->SetMinimum(-100);
      //hh[h]->SetMaximum(0xFFF);
      hh[h]->SetMaximum(4300);


      for (int i=0; i<1024; i++)
	hh[h]->SetBinContent(i+1, adc[a][i]);
      //h->Modified();
      if(h%2==0)
	hh[h]->Draw();
      else {
	hh[h]->SetLineColor(2);
	hh[h]->Draw("same");
      }
    }
#if 0
    {
      gC->cd(1);
      if (!ha) {
	ha = new TH1D("cha", "cha", 1024, 0, 1024);
      }

      ha->Clear();
      ha->SetMinimum(0);
      ha->SetMaximum(0xFFF);

      int a = 0;
      for (int i=0; i<1024; i++)
	ha->SetBinContent(i+1, adc[a][i]);
      //h->Modified();
      ha->Draw();
    }
    {
      gC->cd(2);
      if (!hb) {
	hb = new TH1D("chb", "chb", 1024, 0, 1024);
      }

      hb->Clear();
      hb->SetMinimum(0);
      hb->SetMaximum(0xFFF);

      int a = 1;
      for (int i=0; i<1024; i++)
	hb->SetBinContent(i+1, adc[a][i]);
      //h->Modified();
      hb->Draw();
    }
#endif

    {
      gC->cd(8);
      if (!hc) {
	hc = new TH1D("chc", "chc", 1024, 0, 1024);
      }

      hc->Clear();
      hc->SetMinimum(0);
      hc->SetMaximum(0xFFF);

      int a = 0;
      for (int i=0; i<1024; i++)
	hc->SetBinContent(i+1, adc_tr[a][i]);
      //h->Modified();
      hc->Draw();
    }

    {
      gC->cd(8);
      if (!hd) {
	hd = new TH1D("chd", "chd", 1024, 0, 1024);
      }

      hd->Clear();
      hd->SetMinimum(0);
      hd->SetMaximum(0xFFF);

      int a = 1;
      for (int i=0; i<1024; i++)
	hd->SetBinContent(i+1, adc_tr[a][i]);
      //h->Modified();
      hd->SetLineColor(2);

      hd->Draw("same");
    }

    gC->Modified();
    gC->Draw();
    gC->Update();

    now = time(NULL);
    tx = now + 2;

    v1742_newrun = 0;
  }//end doV1742

  void UpdateHistograms(TDataContainer& dataContainer){

    TV1190Data *data = dataContainer.GetEventData<TV1190Data>("EMMT");
    if(!data) return;
    std::vector<TDCMeasurement> measurements = data->GetMeasurements();
    int chan = -1;
    // Seems to be some noise in the measurements.  In the case of multiple
    // measurements for the same channel, get the earliest measurement.
    //    double a1m_earliest = 9999999.0, a2m_earliest = 9999999.0;
    // Vector of the earliest TDC time for each channel
    std::vector<int> anode_pathway(9,0);
    std::vector<double> earliest_times(64,999999);
    std::vector<int> counts(64,0);
    //    std::vector<double> a1_pulse_time(10,-1.0);
    //    int a1_counter = 0;
    for(unsigned int i = 0; i < measurements.size(); i++){ // loop over measurements
      chan = measurements[i].GetChannel();
      if((GetCurrentRunNumber()>480)&&(GetCurrentRunNumber()<612)) {      
	if((chan>5)&&(chan<10)) {
	  chan=15-chan;
	}
      }
      double t = measurements[i].GetMeasurement()/5.0; // convert to nsec
      counts[chan] = counts[chan] + 1;
      if(t < earliest_times[chan])
	earliest_times[chan] = t;
    }

    // Output a text file of the raw data for offline debugging, note, can generate really big file
    /*    for(int i = 0; i < 14; i++) {
	  if(earliest_times[i] < 999999)
	  raw_output << std::setw(8) << earliest_times[i];
	  else
	  raw_output << std::setw(8) << 0;
	  }
	  raw_output << "\n";*/
    

    // ***************************************
    // RAW DATA
    // ***************************************
    // Compute and fill anode pathway multiplicity data
    for (int i=0; i < 3; i++) {
      if(earliest_times[i] < 999990) {
        if(earliest_times[3] < 999990) {
          anode_pathway[0+3*i] = 1;
        }
        if(earliest_times[4] < 999990) {
          anode_pathway[1+3*i] = 1;
        }
        if(earliest_times[5] < 999990) {
          anode_pathway[2+3*i] = 1;
        }
      }
    }
    int apm = 0;
    for (int i=0; i<9; i++)
      apm = apm + anode_pathway[i]*int(pow(2,i));
    // multiplicity of 1 means APM = 1, 2, 4, 8, 16, 32, 64, 128,or 256
    anode_pathway_multiplicity->Fill(apm);

    // This is used elsewhere
    if (counts[6] > 0)
      t6 = earliest_times[6] - earliest_times[1];
    else
      t6 = 0;

    if (counts[7] > 0)
      t7 = earliest_times[7] - earliest_times[1];
    else
      t7 = 0;

    // Fill m_count debugger data
    int multiplicity = 0;
    for (int i=0; i<9; i++)
      multiplicity = multiplicity + anode_pathway[i];

    m_counts[0]->Fill(multiplicity);
    m_counts[1]->Fill(counts[4]);
    m_counts[2]->Fill(counts[6]);
    m_counts[3]->Fill(counts[7]);
    m_counts[4]->Fill(counts[8]);
    m_counts[5]->Fill(counts[9]);
    m_counts[6]->Fill(counts[10]);
    m_counts[7]->Fill(counts[11]);
    m_counts[8]->Fill(counts[12]);
    m_counts[9]->Fill(counts[13]);
    m_counts[10]->Fill(counts[1]+counts[4]+counts[6]+counts[7]+counts[8]+counts[9]+counts[10]+counts[11]+counts[12]+counts[13]);
    int total = 0;
    for(int i=0; i<64; i++)
      total = total + counts[i];
    m_counts[11]->Fill(total);

    // Fill anode relevant data
    // A1 - A2
    // [0] = A1B - A2B
    // [1] = A1B - A2M
    // [2] = A1B - A2T
    // [3] = A1M - A2B
    // [4] = A1M - A2M
    // [5] = A1M - A2T
    // [6] = A1T - A2B
    // [7] = A1T - A2M
    // [8] = A1T - A2T
    if (multiplicity > 0) {
      for (int i=0; i < 3; i++) {
        if(earliest_times[i] < 999990) {
          if(earliest_times[3] < 999990)
            a1_a2_diff[0+3*i]->Fill(earliest_times[i]-earliest_times[3]);
          if(earliest_times[4] < 999990)
            a1_a2_diff[1+3*i]->Fill(earliest_times[i]-earliest_times[4]);
          if(earliest_times[5] < 999990)
            a1_a2_diff[2+3*i]->Fill(earliest_times[i]-earliest_times[5]);
        }
      }
    }

    // Fill cathode relevant data
    for(int i = 0; i < 4; i++) {
      int index = 6+i*2;
      // Only fill data for straight paths, multiplicity of 1
      //      if (apm == 1 || apm == 16 || apm == 256) {
      if (multiplicity > 0) {
        // First find a valid anode time for corresponding anode
        // This should be the smallest time between the 3 earliest_time,
        // since multiplicity = 1, the other 2 would be 999999
        double anode_time = 999999.0;
        for (int j=0; j<3; j++) {
          if (i == 0 || i == 1) {
            if (earliest_times[j] < anode_time)
              anode_time = earliest_times[j];
          }
          else {
            if (earliest_times[3+j] < anode_time)
              anode_time = earliest_times[3+j];
          }
        }
        // Check both cathode data are non-zero
        if(earliest_times[index] < 999999.0 && earliest_times[index+1] < 999999.0) {
          // Fill x/y_sum - 2*anode_time
          x_y_sum[i]->Fill(earliest_times[index+1]+earliest_times[index] - 2*anode_time);
          // R-L or B-T
          if(i==0 || i==2) {
            x_y_diff[i]->Fill(earliest_times[index+1]-earliest_times[index]);
            x_y_diff_vs_sum_silly[i]->Fill(earliest_times[index+1]+earliest_times[index] - 2*anode_time, earliest_times[index+1]-earliest_times[index]);
          }
          else {
            x_y_diff[i]->Fill(earliest_times[index]-earliest_times[index+1]);
            x_y_diff_vs_sum_silly[i]->Fill(earliest_times[index+1]+earliest_times[index] - 2*anode_time, earliest_times[index]-earliest_times[index+1]);
          }
        }
      }
    }
    doV1742(dataContainer);
  } //end UpdateHistograms

  void PlotCanvas(TDataContainer& dataContainer){

    // ********************************
    // DRAW HISTOGRAMS
    // ********************************
    // A1-A2 Summary tab
    //   4 pads:    [3 straight paths overplot]    [straight combined]
    //              [crossing combined]            [anode pathway multiplicity]
    //
    // A1 - A2
    // [0] = A1B - A2B
    // [1] = A1B - A2M
    // [2] = A1B - A2T
    // [3] = A1M - A2B
    // [4] = A1M - A2M
    // [5] = A1M - A2T
    // [6] = A1T - A2B
    // [7] = A1T - A2M
    // [8] = A1T - A2T
    // [9] = A1 - A2 (All Straight)
    // [10] = A1 - A2 (All Crossing)
    // [11] = A1B - A2B (Copy for plotting purpose)
    // [12] = A1M - A2M (Copy for plotting purpose)
    // [13] = A1T - A2T (Copy for plotting purpose)
 
    if(GetDisplayWindow()->GetCurrentTabName().compare("A1-A2 Summary") == 0) {       
      TCanvas* c1 = GetDisplayWindow()->GetCanvas("A1-A2 Summary");
      c1->Clear();
      c1->Divide(2,2);
      // dumb/quick way of getting sum and getting duplicate copy
      for(int i=0; i<152; i++) {/*
				  a1_a2_diff[9]->SetBinContent(i,a1_a2_diff[0]->GetBinContent(i) + a1_a2_diff[4]->GetBinContent(i) + a1_a2_diff[8]->GetBinContent(i));
				  a1_a2_diff[10]->SetBinContent(i,a1_a2_diff[1]->GetBinContent(i) + a1_a2_diff[2]->GetBinContent(i) + a1_a2_diff[3]->GetBinContent(i) + a1_a2_diff[5]->GetBinContent(i) + a1_a2_diff[6]->GetBinContent(i) + a1_a2_diff[7]->GetBinContent(i));
				  a1_a2_diff[11]->SetBinContent(i,a1_a2_diff[0]->GetBinContent(i));
				  a1_a2_diff[12]->SetBinContent(i,a1_a2_diff[4]->GetBinContent(i));
				  a1_a2_diff[13]->SetBinContent(i,a1_a2_diff[8]->GetBinContent(i));
				*/}
      // multiplicity of 1 means APM = 1, 2, 4, 8, 16, 32, 64, 128,or 256

      c1->cd(1);
      a1_a2_diff[4]->Draw();
      a1_a2_diff[0]->Draw("same");
      a1_a2_diff[8]->Draw("same");
      c1->cd(2);
      a1_a2_diff[9]->Draw();
      c1->cd(3);
      a1_a2_diff[10]->Draw();
      c1->cd(4);
      gPad->SetLogy();
      anode_pathway_multiplicity->Draw();
      c1->Modified();
      c1->Update();
    }

    if(GetDisplayWindow()->GetCurrentTabName().compare("A1-A2 Detail") == 0){       
      TCanvas* c1 = GetDisplayWindow()->GetCanvas("A1-A2 Detail");
      c1->Clear();
      c1->Divide(3,3);
      for(int i = 0; i < 9; i++){
	c1->cd(1+i);
        a1_a2_diff[i]->Draw();
      }
      c1->Modified();
      c1->Update();
    }

    if(GetDisplayWindow()->GetCurrentTabName().compare("M_Counts") == 0){       
      TCanvas* c1 = GetDisplayWindow()->GetCanvas("M_Counts");
      c1->Clear();
      c1->Divide(4,3);
      c1->cd(1);  m_counts[0]->Draw();
      c1->cd(2);  m_counts[1]->Draw();
      c1->cd(3);  m_counts[10]->Draw();
      c1->cd(4);  m_counts[11]->Draw();
      for(int i = 4; i < 12; i++){
	c1->cd(1+i);
	m_counts[i-2]->Draw();
      }
      c1->Modified();
      c1->Update();
    }

    if(GetDisplayWindow()->GetCurrentTabName().compare("XY Diffs") == 0){       
      TCanvas* c1 = GetDisplayWindow()->GetCanvas("XY Diffs");
      c1->Clear();
      c1->Divide(2,2);
      for(int i = 0; i < 4; i++){
	c1->cd(1+i);
	x_y_diff[i]->Draw();
      }
      c1->Modified();
      c1->Update();
    }

    if(GetDisplayWindow()->GetCurrentTabName().compare("XY Sums") == 0){       
      TCanvas* c1 = GetDisplayWindow()->GetCanvas("XY Sums");
      c1->Clear();
      c1->Divide(2,2);
      for(int i = 0; i < 4; i++){
	c1->cd(1+i);
	x_y_sum[i]->Draw();
      }
      c1->Modified();
      c1->Update();
    }

    if(GetDisplayWindow()->GetCurrentTabName().compare("XY Diff vs Sum (Full)") == 0){ 
      TCanvas* c1 = GetDisplayWindow()->GetCanvas("XY Diff vs Sum (Full)");
      c1->Clear();
      if(!(c1->GetShowEventStatus()))c1->ToggleEventStatus();
      if(!(c1->GetShowToolBar()))c1->ToggleToolBar();
      c1->Divide(2,2);
      for(int i = 0; i < 4; i++){
	c1->cd(1+i);
	x_y_diff_vs_sum_silly[i]->Draw("colz");
      }
      c1->Modified();
      c1->Update();
    }

    if(GetDisplayWindow()->GetCurrentTabName().compare("XY Positions") == 0){       
      TCanvas* c1 = GetDisplayWindow()->GetCanvas("XY Positions");
      c1->Clear();
      c1->Divide(2,2);
      for(int i = 0; i < 4; i++){
	c1->cd(1+i);
	x_y_diff_vs_sum_silly[i]->Draw("colz");
      }
      c1->Modified();
      c1->Update();
    }
  }
}; 

int main(int argc, char *argv[])
{
  TStyle *his_disp = new TStyle("his_disp","histogram display style");
  his_disp->SetCanvasBorderMode(0);
  his_disp->SetFrameBorderMode(0);
  his_disp->SetPadBorderMode(0);
  his_disp->SetCanvasColor(0);
  his_disp->SetPadColor(0);
  his_disp->SetStatColor(0);
  his_disp->SetTitleColor(0);
  his_disp->SetOptStat(1111111);
  his_disp->SetPalette(1,0);
  his_disp->cd();
  MyTestLoop::CreateSingleton<MyTestLoop>();  
  return MyTestLoop::Get().ExecuteLoop(argc, argv);
}
