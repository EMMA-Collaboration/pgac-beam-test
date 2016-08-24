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
#include "TH2D.h"
#include "TH3D.h"
#include "TPolyLine3D.h"
#include "TMath.h"

// ROOTANA header files
#include "TRootanaEventLoop.hxx"

#define USE_V1190

#ifdef  USE_V1190
#include <TV1190Data.hxx>
#endif 

//global counting variables
const int Nchan = 14; //number of detector channels
const int Ncounts = 5; //maximum number of multi-hit counts readout from TDC

class Analyzer: public TRootanaEventLoop {
public:
  //1D histograms
  TH1F *hdata[Nchan];
 
  Analyzer() {
    DisableAutoMainWindow();
    // Create histograms (if enabled)
    SetOutputFilename("output/emma_ana_bare_");
    DisableAutoMainWindow();
  };

  virtual ~Analyzer() {};

  void Initialize() {
    for(int i=0;i<Nchan;i++){
      hdata[i]=0;
    }
  }//end Initialize

  void BeginRun(int transition,int run,int time) {
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
      hdata[i]->SetYTitle("Number of Entries");
    }
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
      chan = measurements[i].GetChannel();
      if((GetCurrentRunNumber()>480)&&(GetCurrentRunNumber()<612)) {//switch inputs if second PGAC run
	if((chan>5)&&(chan<10)) {
	  chan=15-chan;
	}
      }
      counts[chan]++; //"counts" uses natural counting; i.e. 0 is 0 counts, 1 is 1 count...
      datum[chan][counts[chan]] = measurements[i].GetMeasurement();
      hdata[chan]->Fill(datum[chan][counts[chan]]);// raw data histograms
    }//end unpack data
    return true;
  }//end ProcessMidasEvent
}; 

int main(int argc, char *argv[])
{
  Analyzer::CreateSingleton<Analyzer>();
  return Analyzer::Get().ExecuteLoop(argc, argv);
}
