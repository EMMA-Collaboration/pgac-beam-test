#include "TV1190Histogram.h"

#include "TV1190Data.hxx"
#include "TDirectory.h"

const int Nchannels = 64;

/// Reset the histograms for this canvas
TV1190Histograms::TV1190Histograms(){  
  
  CreateHistograms();
}


void TV1190Histograms::CreateHistograms(){
  

  // Otherwise make histograms
  clear();
  
  std::cout << "Create Histos" << std::endl;
  for(int i = 0; i < Nchannels; i++){ // loop over channels    

    char name[100];
    char title[100];
    sprintf(name,"V1190_%i_%i",0,i);

    // Delete old histograms, if we already have them
    TH1D *old = (TH1D*)gDirectory->Get(name);
    if (old){
      delete old;
    }


    // Create new histograms
    // OK, define the channel mapping here
    if(i==0)
      sprintf(title,"V1190 A1B (ch=%i)",i);	
    else if(i==1)
      sprintf(title,"V1190 A1M (ch=%i)",i);
    else if(i==2)
      sprintf(title,"V1190 A1T (ch=%i)",i);
    else if(i==3)
      sprintf(title,"V1190 A2B (ch=%i)",i);	
    else if(i==4)
      sprintf(title,"V1190 A2M (ch=%i)",i);
    else if(i==5)
      sprintf(title,"V1190 A2T (ch=%i)",i);
    else if(i==6)
      sprintf(title,"V1190 X1L (ch=%i)",i);
    else if(i==7)
      sprintf(title,"V1190 X1R (ch=%i)",i);
    else if(i==8)
      sprintf(title,"V1190 Y1B (ch=%i)",i);
    else if(i==9)
      sprintf(title,"V1190 Y1T (ch=%i)",i);
    else if(i==10)
      sprintf(title,"V1190 X2L (ch=%i)",i);
    else if(i==11)
      sprintf(title,"V1190 X2R (ch=%i)",i);
    else if(i==12)
      sprintf(title,"V1190 Y2B (ch=%i)",i);
    else if(i==13)
      sprintf(title,"V1190 Y2T (ch=%i)",i);
    else
      sprintf(title,"V1190 unused (ch=%i)",i);
      

    TH1D *tmp = new TH1D(name,title,10000,0,4000);
    tmp->SetXTitle("Time (ns)");
    tmp->SetYTitle("Number of Entries");
    push_back(tmp);
  }

}


 #include <sys/time.h>

  
/// Update the histograms for this canvas.
void TV1190Histograms::UpdateHistograms(TDataContainer& dataContainer){


  TV1190Data *data = dataContainer.GetEventData<TV1190Data>("EMMT");
  if(!data) return;
  
  /// Get the Vector of ADC Measurements.
  std::vector<TDCMeasurement> measurements = data->GetMeasurements();
  for(unsigned int i = 0; i < measurements.size(); i++){ // loop over measurements
	
    int chan = measurements[i].GetChannel();
    if(GetHistogram(chan))
      GetHistogram(chan)->Fill(measurements[i].GetMeasurement()/5.0);
  }

}



/// Take actions at begin run
void TV1190Histograms::BeginRun(int transition,int run,int time){

  CreateHistograms();

}

/// Take actions at end run  
void TV1190Histograms::EndRun(int transition,int run,int time){

}
