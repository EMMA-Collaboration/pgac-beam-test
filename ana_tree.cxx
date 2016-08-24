// Default program for dealing with various standard TRIUMF VME setups:
// V792, V1190 (VME), L2249 (CAMAC), Agilent current meter
//
//

// C++ header files
#include <stdio.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <math.h>

// ROOT header files
#include "TH1D.h"
//#include "TH2D.h"
//#include "TH3D.h"
//#include "TPolyLine3D.h"
#include "TMath.h"
#include "TTree.h"
#include "TBranch.h"

// ROOTANA header files
#include "TRootanaEventLoop.hxx"

#define USE_V1190

#ifdef  USE_V1190
#include <TV1190Data.hxx>
#endif 

//global counting variables
const int Nchan = 14; //number of detector channels
const int Ncounts = 5; //set maximum number of multiple hits per channel to read out

class Analyzer: public TRootanaEventLoop {
public:
  //1D histograms
  TH1F *hdata[Nchan];
  TH1F *hmsize;

  TTree *tree1; // Raw data tree

  /* Structs for tree */
  struct PGAC_t{
    Int_t chanPGAC[7];
  };
  
  struct PGACsig_t{
    Int_t ab;
    Int_t am;
    Int_t at;
    Int_t xn;
    Int_t xf;
    Int_t yn;
    Int_t yf;
  };

  PGAC_t PGAC[2];
  
  Analyzer() {
    DisableAutoMainWindow();
    // Create histograms (if enabled)
    SetOutputFilename("output_tree/emma_ana_");
    DisableAutoMainWindow();
  };

  virtual ~Analyzer() {};

  void Initialize() {
    for(int i=0;i<Nchan;i++){
      hdata[i]=0;
    }
  }//end Initialize

  void BeginRun(int transition,int run,int time) {
    tree1 = new TTree("t1","Raw Data Tree, array");
    tree1->Branch("PGAC1",&PGAC[0],"ab/I:am:at:xn:xf:yn:yf");
    tree1->Branch("PGAC2",&PGAC[1],"ab/I:am:at:xn:xf:yn:yf");
    
    //histogram naming variables
    TString titles[Nchan]={"A1B","A1M","A1T","A2B","A2M","A2T","X1L","X1R","Y1B","Y1T","X2L","X2R","Y2B","Y2T"};
    char name[100];
    TString title;
    TString subtitle;
    TString subsubtitle;
    //hisogram ranges
    //data range, full range, un-calibrated 
    Float_t dmin=0;
    Float_t dmax=2048*10;//4096*5;
    Float_t dbin=2048;

    //hisogram definitions         
    for(int i = 0; i < Nchan; i++) { // debugging histograms, loop over channels    
      title=titles[i];
      title+=" V1190 Ch.";
      title+=i;
      
      sprintf(name,"hdata%i",i);
      subtitle=" Data";

      hdata[i] = new TH1F(name,title+subtitle,dbin,dmin,dmax);
      hdata[i]->SetXTitle("Channel");
      hdata[i]->SetYTitle("Number of entries");
    }
    Int_t mbins=30;
    hmsize = new TH1F("hmsize","Measurement Size",mbins,0,mbins);
    hmsize->SetXTitle("Measurement size");
    hmsize->SetYTitle("Number of entries");
  }

  bool ProcessMidasEvent(TDataContainer& dataContainer) {
    TV1190Data *data = dataContainer.GetEventData<TV1190Data>("EMMT");
    if(!data) return false;
    std::vector<TDCMeasurement> measurements = data->GetMeasurements();
    int chan = 0; //channel number
    int datum[Nchan][Ncounts] = {{0}}; //data word
    std::vector<int> counts(Nchan,0); //number of counts (multi-hit)

    //--- Unpack data
    for(unsigned int i = 0; i < measurements.size(); i++) { // loop over measurements
      hmsize->Fill(measurements.size());
      chan = measurements[i].GetChannel();
      if((GetCurrentRunNumber()>480)&&(GetCurrentRunNumber()<612)) {//switch inputs if second PGAC run
	if((chan>5)&&(chan<10)) {
	  chan=15-chan;
	}
      }
      counts[chan]++; //"counts" uses natural counting; i.e. 0 is 0 counts, 1 is 1 count...
      datum[chan][counts[chan]] = measurements[i].GetMeasurement();
      hdata[chan]->Fill(datum[chan][counts[chan]]);// raw data histograms
      /*      
	      switch(chan){
	      case 0 ://A1B, ab
	      PGAC[0].chanPGAC[0]=datum[chan][counts[chan]];
	      break;
	      case 1 ://A1M, am
	      PGAC[0].chanPGAC[1]=datum[chan][counts[chan]];
	      break;
	      case 2 ://A1T, at
	      PGAC[0].chanPGAC[2]=datum[chan][counts[chan]];
	      break;
	      case 3 ://A2B, ab
	      PGAC[1].chanPGAC[0]=datum[chan][counts[chan]];
	      break;
	      case 4 ://A2M, am
	      PGAC[1].chanPGAC[1]=datum[chan][counts[chan]];
	      break;
	      case 5 ://A2T, at
	      PGAC[1].chanPGAC[2]=datum[chan][counts[chan]];
	      break;
	      case 6 ://X1L, xf
	      PGAC[0].chanPGAC[4]=datum[chan][counts[chan]];
	      break;
	      case 7 ://X1R, xn
	      PGAC[0].chanPGAC[3]=datum[chan][counts[chan]];
	      break;
	      case 8 ://Y1B, yn
	      PGAC[0].chanPGAC[5]=datum[chan][counts[chan]];
	      break;
	      case 9 ://Y1T, yf
	      PGAC[0].chanPGAC[6]=datum[chan][counts[chan]];
	      break;
	      case 10 ://X2L, xf
	      PGAC[1].chanPGAC[4]=datum[chan][counts[chan]];
	      break;
	      case 11 ://X2R, xn
	      PGAC[1].chanPGAC[3]=datum[chan][counts[chan]];
	      break;
	      case 12 ://Y2B, yn
	      PGAC[1].chanPGAC[5]=datum[chan][counts[chan]];
	      break;
	      case 13 ://Y2T, yf
	      PGAC[1].chanPGAC[6]=datum[chan][counts[chan]];
	      break;
	      }
	      tree1->Fill();
      */
    }//end unpack data

    //--- Sort data
    Bool_t non_zero_event=kFALSE;
    Int_t maxcount=Ncounts;
    for(Int_t j=1;j<=maxcount;j++) {//loop over counts
	for(int i=0;i<Nchan/2;i++) {//loop over pairs of signals
	  if(i<3) {//anode data
	    PGAC[0].chanPGAC[i]=datum[i][j];
	    PGAC[1].chanPGAC[i]=datum[i+3][j];
	  }
	  else {//cathode data 
	    if(i%2) {// x-signals (the x-signals are labeled from downstream; I am "flipping" them 
	      //        to keep the relative positions defined from upstream: XR is "left" and XL 
	      //        is "right") (i=3,5)
	      PGAC[(int)((i-3)/2)].chanPGAC[3]=datum[(2*i)+1][j]; // "near" signal; labeled right, beam left
	      PGAC[(int)((i-3)/2)].chanPGAC[4]=datum[(2*i)+0][j]; // "far" signal; labeled left, beam right
	    }
	    else    {// y-signals (i=4,6)
	      PGAC[(int)((i-3)/2)].chanPGAC[5]=datum[(2*i)+0][j]; // "near" signal; bottom
	      PGAC[(int)((i-3)/2)].chanPGAC[6]=datum[(2*i)+1][j]; // "far" signal; top
	    }	
	  }//end cathode data
	}//end pairs loop
	non_zero_event=kFALSE;
	for(Int_t k=0;k<Nchan;k++) {
	  non_zero_event+=datum[k][j];
	}
	if(non_zero_event)
	  tree1->Fill();//only fill trees with non-zero data
    }//end counts loop

    return true;
  }//end ProcessMidasEvent
}; 

int main(int argc, char *argv[])
{
  Analyzer::CreateSingleton<Analyzer>();
  return Analyzer::Get().ExecuteLoop(argc, argv);
}
