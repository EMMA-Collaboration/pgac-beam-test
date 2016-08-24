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
//#include "TV1190Histogram.h"
#endif 

const int Nchan = 14; //number of detector channels
const int Ncounts = 5; //maximum number of multi-hit counts readout from TDC
const Float_t Tcal = 5.0; //calibration from number of channels to ns for the TDC
int maxcount = 0; //running count of maximum multi-hit numbrer
const int Nray = 10;
int nevent = 0; //event counter

//flags
const bool docal=true;//turn calibration on
Bool_t dogate=1;//turn gating on
const bool free_range=false;//set histogram to auto-range
const bool match_sim=true;//match histogram ranges and binning to simulation
Bool_t bdiag=kTRUE;
Bool_t doprint=false;
Bool_t set_showtest=kFALSE;
//Bool_t set_showoutput=kFALSE;

class Analyzer: public TRootanaEventLoop {

#ifdef  USE_V1190
  //  TV1190Histograms *v1190_histos;
#endif 

public:
  TH1F *hcounts[Nchan+2];  
  TH1F *hdata[Nchan];
  TH2F *hxfxn[Nchan/2];
  TH2F *hxfxnz[Nchan/2];
  TH2F *hxfxnzn[Nchan/2];
  TH2F *hxfxnznn[Nchan/2];
  
  TH1F *hsum[Nchan/2];
  TH1F *hsumz[Nchan/2];
  TH1F *hsumzg[Nchan/2];
  TH1F *hdiff[Nchan/2];
  TH1F *hdiffz[Nchan/2];
  TH1F *htime[9];
  TH1F *htimez[9];
  TH2F *haa[9];
  TH2F *haaz[9];
  TH1F *hmhit;
  TH1F *haacoinc;
  TH2F *htev[3];
  
  TH1F *hx[4];
  TH1F *hxc[4];
  TH1F *hxcg[4];
  
  //1D diagnostic plots
  TH1F *hang[4];
  //  TH1F *hangm[4];
  TH1F *hangt[4];
  TH1F *hxdiff[2];
  TH1F *hxcdiff[2];
  TH1F *hrho[2];
  TH2F *hxx[2];
  TH2F *hxxc[2];
    
  //2D position spectra
  TH2F *hhit[2];
  TH2F *hhitc[2];
  //TH2F *hhitcg[2];
  
  //2D ray tracing
  const static int planes=24;
  Float_t zp[planes];
  Float_t Z[5];
  Float_t Zp[4];
  Float_t offset_cal[2]  ;
  TH2F *hsource[planes];
  TH3F *htrace;
  TH3F *htracec;
  TH2F *hzx;
  TH2F *hzy;
  
  //2D plots with anode or time
  TH2F *haxf[12];
  TH2F *haxn[12];
  TH2F *hasum[12];
  TH2F *hasumn[12];
  TH2F *hasumnn[12];
  TH2F *hadiff[12];
  TH2F *hax[12];
  TH2F *haxc[12];
  TH2F *htx[12];
  TH2F *htxc[12];
  TH2F *htsum[12];
  
  //2D plots with sum
  //  TH2F *hsa[4];//redundant with hasum and poorly definedand and filled
  TH2F *hsdiff[Nchan/2];  
  TH2F *hsx[Nchan/2];
  TH2F *hsxc[Nchan/2];
  TH2F *hss[4];

  //Calibration constants
  Float_t hxxcal[Nchan/2][2];//linear calibration constants
  Float_t hxcal[12][2];
  Float_t hasumcal[12][3];
  Float_t hsumzcal[Nchan/2][3];
  Float_t hdiffzcal[Nchan/2][3];

  Analyzer() {
    DisableAutoMainWindow();
    // Create histograms (if enabled)
    SetOutputFilename("output/emma_ana_");
    DisableAutoMainWindow();
#ifdef  USE_V1190
    //   v1190_histos = new TV1190Histograms();
#endif 
  };

  virtual ~Analyzer() {};

  void Initialize() {
    if(bdiag){
      for(int i=0;i<Nchan;i++){
	hcounts[i]=0;
	hdata[i]=0;
      }
    }
    for(int i=0;i<Nchan/2;i++){
      if(bdiag) {
	hxfxn[i]=0;
	hxfxnzn[i]=0;
	hxfxnznn[i]=0;
      }
      hxfxnz[i]=0;
      hxxcal[i][0]=0;//offset
      if(i<3)
	hxxcal[i][1]=1;//slope (anode)
      else
	hxxcal[i][1]=-1;//slope (cathode)
      hsumzcal[i][0]=0;//mean 
      hsumzcal[i][1]=0;//sigma
      hsumzcal[i][2]=0;//width
      hdiffzcal[i][0]=0;//mean 
      hdiffzcal[i][1]=0;//sigma
      hdiffzcal[i][2]=0;//width
    }
    for(int i=0;i<2;i++){
      hhit[i]=0;
      hhitc[i]=0;
      //hhitcg[i]=0;
    }
    if(bdiag){
      for(int i=0;i<planes;i++){
	hsource[i]=0;
      }
    }
    for(int i=0;i<12;i++){
      if(i<4) {
	hx[i]=0;
	hxc[i]=0;
	hxcg[i]=0;
      }
      hxcal[i][0]=1;   //slope
      hxcal[i][1]=0;   //offset
      hasumcal[i][0]=0;//offset (not needed)
      hasumcal[i][1]=1;//slope
      hasumcal[i][2]=0;//re-centering
    }
    Z[0]=722.9; //X1 z-position
    Z[1]=716.5; //Y1 z-position
    Z[2]=686.1; //X2  z-position
    Z[3]=679.7; //Y2  z-position
    Z[4]=0;    //target z-position

    Zp[0]=(Z[0]+Z[1])/2;//anode 1 position
    Zp[1]=(Z[2]+Z[3])/2;//anode 2 position
    Zp[2]=637.0;        //mask position
    Zp[3]=Z[4];         //target position

    offset_cal[0]=118;
    offset_cal[1]=54/2;

    Float_t start = 0;
    Float_t stop = Z[3];    
    for(int i=0;i<planes;i++){
      zp[i]=start+i*((stop-start)/planes);
    }
    
  }//end Initialize

  bool testfile(const char *filename="test.cal", Bool_t showtest=1, const int numdet=24, int start_no=1, Bool_t showoutput=kTRUE)
  {
    Float_t param[numdet][50]; 
    Int_t errorline=-1;
    Int_t size=sizeof(param[0])/sizeof(param[0][0]);
    if(showtest)
      printf("  param array size is: [%d][%d].\n",(int)(sizeof(param)/sizeof(param[0])),size); 
    Bool_t fit=kFALSE;
    for(Int_t i=0;i<numdet;i++)
      for(Int_t j=0;j<size;j++)
	param[i][j]=0;//initializes all array elements to zero
 
    FILE * infile;
    
    Int_t k=1;//array length
    while(!fit&&k<=(size)) {
      infile = fopen (filename,"r");
      for(Int_t i=0;i<numdet;i++) {
	for(Int_t j=0;j<k;j++) {
	  fscanf(infile,"%f",&param[i][j]);
	}
      }
      fclose(infile);
      if(showtest) printf("  Testing array length %d:\n",k);
      if (param[0][0]==start_no) {
	fit=kTRUE;
	for(Int_t i=0;i<numdet;i++) {
	  fit=(fit&&(param[i][0]==(i+start_no)));	
	  if(fit){
	    if(showtest) printf("   %2.0f ",param[i][0]);
	    for(Int_t j=1;j<k;j++){
	      if(showtest) printf("%11.5f ",param[i][j]);
	    }
	    if(showtest) printf("\n");
	    if((i+1)>errorline)errorline=i+1;
	  }
	}
      }
      else {
	fit=kFALSE;
	k=size;
	printf("  Each line of \"%s\" is expected to start with a detector number.\n",filename);
	printf("  File \"%s\" is expected to start with %d.\n",filename,start_no);
      }
      k++;
    }//end while
    if(fit) {
      printf("  File \"%s\" has %d elements per line.\n",filename,k-1);
      if(showoutput)
	for(Int_t i=0;i<numdet;i++) {
	  printf("   %2.0f ",param[i][0]);
	  for(Int_t j=1;j<k-1;j++){
	    printf("%11.5f ",param[i][j]);
	  }
	  printf("\n");
	}
      return true;
    }   
    else{
      printf("  File \"%s\" has more than %d elements per line, or there is an error on line %d.\n",filename,size,errorline);
      return false;
    }
  }

  void readcal(Int_t run_no) {
    printf("Loading calibration data for run number %d...\n",run_no);
    Float_t detno=0;
    //1. Load gain-matching calibration file
    TString calfile = "cal/run_";
    calfile+=run_no;
    calfile+="_hxx.cal";
    ifstream dummy;
    dummy.open(calfile);
    if(!dummy.good()) {
      printf(" Input file %s not found.  Attempting default file...\n",calfile.Data());
      if((GetCurrentRunNumber()>480)&&(GetCurrentRunNumber()<612)) {
	calfile="cal/run_606_hxx.cal"; 
      }
      else
	calfile="cal/run_430_hxx.cal"; 
    }
    dummy.close();   
    dummy.clear();
    ifstream infile(calfile);
    if(infile.good()) {
      Int_t array_size=sizeof(hxxcal[0])/sizeof(hxxcal[0][0]);
      Int_t det_num=sizeof(hxxcal)/sizeof(hxxcal[0]);
      printf(" Input calibration file is: %s\n",calfile.Data());
      if(set_showtest)
	printf("  hxxcal array size is: [%d][%d].\n",det_num,array_size);
      if(testfile(calfile.Data(),set_showtest,det_num,0)) {
	for(Int_t i=0;i<det_num;i++) {
	  infile>>detno;//First number in each row is detector number
	  for(Int_t j=0;j<array_size;j++){//loop over number of fit parameters
	    infile>>hxxcal[i][j];
	  }
	}
      }
    }
    else printf(" Input file %s not found!\n",calfile.Data());
    infile.close();	
  
    //load calibration file for anode-cathode gain matching
    calfile = "cal/run_";
    calfile+=run_no;
    calfile+="_hasum.cal";
    dummy.open(calfile);
    if(!dummy.good()) {
      printf(" Input file %s not found.  Attempting default file...\n",calfile.Data());
      if((GetCurrentRunNumber()>480)&&(GetCurrentRunNumber()<612)) {
	calfile="cal/run_606_hasum.cal"; 
      }
      else
	calfile="cal/run_480_hasum.cal"; 
    }
    dummy.close();   
    dummy.clear();    
    ifstream infile3(calfile);
    if(infile3){
      printf(" Input calibration file is: %s\n",calfile.Data());
      Int_t det_num=12;
      if(testfile(calfile.Data(),set_showtest,det_num,0)) {
	for(Int_t i=0;i<det_num;i++) {//loop over number of sets
	  infile3>>detno;//First number in each row is detector number
	  for(Int_t j=0;j<3;j++){//loop over number of fit parameters
	    infile3>>hasumcal[i][j];
	  }
	}
      }
    }
    else printf(" Input file %s not found!\n",calfile.Data());
    infile3.close();

    //load position calibration file
    calfile = "cal/run_";
    calfile+=run_no;
    calfile+="_XY.cal";
    dummy.open(calfile);
    if(!dummy.good()) {
      printf(" Input file %s not found.  Attempting default file...\n",calfile.Data());
      if((GetCurrentRunNumber()>480)&&(GetCurrentRunNumber()<612)) {
	calfile="cal/run_606_XY.cal"; 
      }
      else
	calfile="cal/run_480_XY.cal"; 
    }
    dummy.close();   
    dummy.clear();
    ifstream infile2(calfile);
    if(infile2){
      printf(" Input calibration file is: %s\n",calfile.Data());
      Int_t det_num=12;      
      if(testfile(calfile.Data(),set_showtest,det_num,0)) {
	for(Int_t i=0;i<det_num;i++) {
	  infile2>>detno;//First number in each row is detector number
	  for(Int_t j=0;j<2;j++){//loop over number of fit parameters
	    infile2>>hxcal[i][j];
	  }
	}
      }
    }
    else printf("  Input file %s not found!\n",calfile.Data());
    infile2.close();
     
    if(dogate) {
      //load gating parameters
      calfile = "cal/run_";
      calfile+=run_no;
      calfile+="_hsumz.cal";
      ifstream infile4(calfile);
      if(infile4) {
	printf(" Input calibration file is: %s\n",calfile.Data());
	if(testfile(calfile.Data(),set_showtest,7,0)) {
	  for(Int_t i=0;i<7;i++) {//loop over number of sets
	    infile4>>detno;//First number in each row is detector number
	    for(Int_t j=0;j<3;j++){//loop over number of fit parameters
	      infile4>>hsumzcal[i][j];
	    }
	  }
	}
      }
      else {
	printf(" Input file %s not found!\n  Gating disabled.\n",calfile.Data());
	dogate=kFALSE;
      }
      infile4.close();
      //load gating parameters
      calfile = "cal/run_";
      calfile+=run_no;
      calfile+="_hdiffz.cal";
      ifstream infile5(calfile);
      if(infile5) {
	printf(" Input calibration file is: %s\n",calfile.Data());
	if(testfile(calfile.Data(),set_showtest,7,0)) {
	  for(Int_t i=0;i<7;i++) {//loop over number of sets
	    infile5>>detno;//First number in each row is detector number
	    for(Int_t j=0;j<3;j++){//loop over number of fit parameters
	      infile5>>hdiffzcal[i][j];
	    }
	  }
	}
      }
      else {
	printf(" Input file %s not found!\n  Gating disabled.\n",calfile.Data());
	dogate=kFALSE;
      }
      infile5.close();
    }
    else
      printf(" Gating disabled.\n");
  }//end readcal

  void BeginRun(int transition,int run,int time){
    if(docal)
      readcal(GetCurrentRunNumber());
    else
      printf("No calibration applied.\n");
    char name[100];
    char pretitle[100];
    TString title;
    TString subtitle;
    TString subsubtitle;
    //hisogram ranges
    //data range, full range, un-calibrated 
    Float_t dmin=0;
    Float_t dmax=2048*10;//4096*5;
    Float_t dbin=2048;
    //data range, zoomed, un-calibrated
    Float_t xmin=9300;
    Float_t xmax=12100;
    Float_t xbin=(Int_t)(xmax-xmin);    
    Float_t cal_bin=500*1;
    Float_t cal_max=67;
    Float_t cal_x_off=offset_cal[0];
    Float_t cal_y_off=offset_cal[1];
    Float_t xbin2D=(Int_t)(xbin/4);
    //time, zoomed and scaled data range
    Float_t tmin=(Int_t)xmin/Tcal;
    Float_t tmax=(Int_t)xmax/Tcal;
    Float_t tbin=(Int_t)xbin/Tcal;
    //range of calibrated position
    Float_t pmin=0;
    Float_t pmax=0;   
    Float_t pxmin=-8-6;
    Float_t pxmax=174-6;
    Float_t pymin=-3-6;
    Float_t pymax=69-6;
    Float_t pzmin=0;
    Float_t pzmax=730;
    Float_t pzbin=4;
    if(match_sim){
      pxmin=-cal_max+cal_x_off;
      pxmax= cal_max+cal_x_off;
      pymin=-cal_max+cal_y_off;
      pymax= cal_max+cal_y_off;
    }

    //range of uncalibrated positions
    Float_t xxmin=0.48;
    Float_t xxmax=0.55;
    Float_t xymin=0.475;
    Float_t xymax=0.525;
    //range of differences
    Float_t xdmin=0;
    Float_t xdmax=0;
    Float_t xdbin=0;
    //range of sum
    Float_t xsmin=19000;
    Float_t xsmax=22500;
    Float_t xsbin=(xsmax-xsmin);
    //range of anodes
    Float_t amin=9200;//9600;
    Float_t amax=10400;//10050;
    Int_t abin=(Int_t)(amax-amin);
    Float_t adiffmin=-100;
    Float_t adiffmax=-adiffmin;

    Float_t xsource=0;
    Float_t xbin3D=256;
    TString titles[Nchan]={"A1B","A1M","A1T","A2B","A2M","A2T","X1L","X1R","Y1B","Y1T","X2L","X2R","Y2B","Y2T"};
 
    if(bdiag){
      hmhit = new TH1F("hmhit","hmhit",10,0,10);
      haacoinc = new TH1F("haacoinc","haacoinc",9,0,9);
      
      for(int i = 0; i < Nchan; i++) { // debugging histograms, loop over channels    
	title=titles[i];
	title+=" V1190 Ch.";
	title+=i;
      
	sprintf(name,"hcounts%i",i);
	subtitle=" Counts";
	hcounts[i] = new TH1F(name,title+subtitle,Ncounts,0,Ncounts);
	hcounts[i]->SetXTitle("Hits per Event");
	hcounts[i]->SetYTitle("Number of Entries");
 
	sprintf(name,"hdata%i",i);
	subtitle=" Data";
	if(free_range)
	  hdata[i] = new TH1F(name,title+subtitle,dbin,1e5,-1e5);
	else
	  hdata[i] = new TH1F(name,title+subtitle,dbin,dmin,dmax);
	hdata[i]->SetXTitle("Channel");
	hdata[i]->SetYTitle("Number of Entries");
      }
      sprintf(name,"hcounts%i",Nchan);
      hcounts[Nchan]=new TH1F(name,"Total Counts",Nchan,0,Nchan);
      hcounts[Nchan]->SetXTitle("Channel No.");
      hcounts[Nchan]->SetYTitle("Number of Entries");    
      sprintf(name,"hcounts%i",Nchan+1);
      hcounts[Nchan+1]=new TH1F(name,"Coincidences",Nchan/2,0,Nchan/2);
      hcounts[Nchan+1]->SetXTitle("Channel Pair");
      hcounts[Nchan+1]->SetYTitle("Number of Entries");
    }

    for(int i = 0; i < Nchan/2; i++) {//correlation plots 
      if(i==0)
	sprintf(pretitle,"A1B vs A2B");
      if(i==1)
	sprintf(pretitle,"A1M vs A2M");
      if(i==2)
	sprintf(pretitle,"A1T vs A2T");
      if(i==3)
	sprintf(pretitle,"X1L vs X1R"); 
      if(i==4)
	sprintf(pretitle,"Y1T vs Y1B"); 
      if(i==5)
	sprintf(pretitle,"X2L vs X2R");
      if(i==6)
	sprintf(pretitle,"Y2T vs Y2B");
      title=pretitle;
    
      sprintf(name,"hxfxn%i",i);   
      subtitle="";
      if(bdiag)
	hxfxn[i] = new TH2F(name,title+subtitle,dbin,dmin,dmax,dbin,dmin,dmax);
     
      if(i<3) {//anode pairs
	xmin=amin;
	xmax=amax;
	xbin2D=(Int_t)(xmax-xmin);    
	if(bdiag){
	  hxfxn[i]->SetXTitle("Detector 2 Anode Channel No.");
	  hxfxn[i]->SetYTitle("Detector 1 Anode Channel No.");
	}
      }
      else{
	if(bdiag){
	  hxfxn[i]->SetXTitle("Channel No. (far side)");
	  hxfxn[i]->SetYTitle("Channel No. (near side)");
	}
      }
      // if((i==3)||(i==5)) {//X-position pairs
      if(i>2) {
	xmin= 9300;
	xmax=12100;
	xbin2D=(Int_t)(xmax-xmin);    
      }
      /* if((i==4)||(i==6)) {//Y-position pairs
	 xmin=8200;
	 xmax=11000;
	 xbin2D=(Int_t)(xmax-xmin);    
	 }*/
      while(xbin2D>512){
	xbin2D/=2;
      }
      
      sprintf(name,"hxfxnz%i",i);   
      subtitle=" Zoom";
      hxfxnz[i] = new TH2F(name,title+subtitle,xbin2D,xmin,xmax,xbin2D,xmin,xmax);
      if(i<3){
	hxfxnz[i]->SetXTitle("Channel No.");
	hxfxnz[i]->SetYTitle("Channel No.");
      }
      else {
	hxfxnz[i]->SetXTitle("Channel No. (near side)");
	hxfxnz[i]->SetYTitle("Channel No. (far side)");
      }
      if(bdiag) {
	sprintf(name,"hxfxnzn%i",i);   
	subtitle=" Zoomed, Gain-matched";
	hxfxnzn[i] = new TH2F(name,title+subtitle,xbin2D,xmin,xmax,xbin2D,xmin,xmax);
	if(i<3){
	  hxfxnzn[i]->SetXTitle("Detector 2 Anode Channel No.");
	  hxfxnzn[i]->SetYTitle("Detector 1 Anode Channel No.");
	}
	else {
	  hxfxnzn[i]->SetXTitle("Channel No. (near side)");
	  hxfxnzn[i]->SetYTitle("Channel No. (far side)");
	}
	sprintf(name,"hxfxnznn%i",i);   
	subtitle=" Zoomed, Gain-matched";
	hxfxnznn[i] = new TH2F(name,title+subtitle,xbin2D,xmin,xmax,xbin2D,xmin,xmax);
	if(i<3){
	  hxfxnznn[i]->SetXTitle("Detector 2 Anode Channel No.");
	  hxfxnznn[i]->SetYTitle("Detector 1 Anode Channel No.");
	}
	else {
	  hxfxnznn[i]->SetXTitle("Channel No. (near side)");
	  hxfxnznn[i]->SetYTitle("Channel No. (far side)");
	}
      }

      sprintf(name,"hsum%i",i);   
      subtitle=", Sum";
      if(free_range)
	hsum[i] = new TH1F(name,title+subtitle,2*xbin2D,1e5,-1e5);
      else
	hsum[i] = new TH1F(name,title+subtitle,2*xbin2D,0,2*dmax);
      hsum[i]->SetXTitle("Sum of Channels");
      hsum[i]->SetYTitle("Number of Entries");

      sprintf(name,"hdiff%i",i);   
      subtitle=", Difference";
      if(free_range)
	hdiff[i] = new TH1F(name,title+subtitle,xbin2D,1e5,-1e5);
      else
	hdiff[i] = new TH1F(name,title+subtitle,xbin2D,-xbin,xbin);
      hdiff[i]->SetXTitle("Difference of Channels");
      hdiff[i]->SetYTitle("Number of Entries");
      
      if(i<3) {//anode signals
	xdmin=adiffmin;//-300;
	xdmax=adiffmax;
	xdbin=(Int_t)(xdmax-xdmin); 

	tmin=xdmin/5.;
	tmax=xdmax/5.;
	tbin=xdbin;

	if(bdiag){
	  sprintf(name,"htev%i",i);
	  subtitle=" Time vs Event No.";
	  htev[i] = new TH2F(name,title+subtitle,1000,0,1e6,xdbin,xdmin,xdmax);
	  htev[i]->SetXTitle("Event No.");
	  htev[i]->SetYTitle("Difference in Anode Channels (Time)");
	}
      }
      else{//cathode signals
	if(i%2==0){//Y signals
	  xdmin=-1000;	 
	  xdmax=1000;	
	}
	else{//X signals
	  xdmin=-800;
	  xdmax=2100;	
	}
	xdbin=(Int_t)(xdmax-xdmin);
	//while(xdbin>1024){xdbin/=2;}
      }
      sprintf(name,"hdiffz%i",i);   
      subtitle=", Difference zoomed";
      hdiffz[i] = new TH1F(name,title+subtitle,xdbin,xdmin,xdmax);
      hdiffz[i]->SetXTitle("Difference of Channels");
      hdiffz[i]->SetYTitle("Number of Entries");

      sprintf(name,"hsumz%i",i);   
      subtitle=", Sum zoomed";
      hsumz[i] = new TH1F(name,title+subtitle,xsbin,xsmin,xsmax);
      hsumz[i]->SetXTitle("Sum of Channels");
      hsumz[i]->SetYTitle("Number of Entries");
     
      sprintf(name,"hsumzg%i",i);   
      subtitle=", Sum zoomed, gated";
      hsumzg[i] = new TH1F(name,title+subtitle,xsbin,xsmin,xsmax);
      hsumzg[i]->SetXTitle("Sum of Channels");
      hsumzg[i]->SetYTitle("Number of Entries");
     
      sprintf(name,"hsdiff%i",i);   
      subtitle=", Sum vs Difference";
      if(free_range)
	hsdiff[i] = new TH2F(name,title+subtitle,xdbin,1e5,-1e5,abin,1e5,-1e5);
      else
	//hsdiff[i] = new TH2F(name,title+subtitle,xdbin,xdmin,xdmax,xsbin,xsmin,xsmax);
	hsdiff[i] = new TH2F(name,title+subtitle,xbin2D,xdmin,xdmax,xbin2D,xsmin,xsmax);
      hsdiff[i]->SetXTitle("Difference of Cathode Channels");
      hsdiff[i]->SetYTitle("Sum of Cathode Channels");

    }//end of correlation plots

    for(int i = 0; i < 9; i++) {//anode correlation plots 
      if(i==0)
	sprintf(pretitle,"A1B vs A2B");
      if(i==1)
	sprintf(pretitle,"A1B vs A2M");
      if(i==2)
	sprintf(pretitle,"A1B vs A2T");
      if(i==3)
	sprintf(pretitle,"A1M vs A2B");
      if(i==4)
	sprintf(pretitle,"A1M vs A2M");
      if(i==5)
	sprintf(pretitle,"A1M vs A2T");
      if(i==6)
	sprintf(pretitle,"A1T vs A2B");
      if(i==7)
	sprintf(pretitle,"A1T vs A2M");
      if(i==8)
	sprintf(pretitle,"A1T vs A2T");
      title=pretitle;
     
      if(bdiag) {
	sprintf(name,"haa%i",i);   
	subtitle=" full range";
	haa[i] = new TH2F(name,title+subtitle,dbin,dmin,dmax,dbin,dmin,dmax);
	haa[i]->SetXTitle("Detector 2 Anode Channel No.");
	haa[i]->SetYTitle("Detector 1 Anode Channel No.");
            
	sprintf(name,"htime%i",i);   
	subtitle=", Difference (time)";
	if(free_range)
	  htime[i] = new TH1F(name,title+subtitle,abin,1e5,-1e5);
	else
	  htime[i] = new TH1F(name,title+subtitle,dbin/8,-dmax/10,dmax/10);
	htime[i]->SetXTitle("Difference of Channels");
	htime[i]->SetYTitle("Number of Entries");
      }

      sprintf(name,"haaz%i",i);   
      subtitle=" Zoom";
      haaz[i] = new TH2F(name,title+subtitle,abin,amin,amax,abin,amin,amax);
      haaz[i]->SetXTitle("Detector 2 Anode Channel No.");
      haaz[i]->SetYTitle("Detector 1 Anode Channel No.");
      
      xdbin=(Int_t)(adiffmax-adiffmin);
      
      sprintf(name,"htimez%i",i);   
      subtitle=", Difference zoomed (time)";
      htimez[i] = new TH1F(name,title+subtitle,xdbin,adiffmin,adiffmax);
      htimez[i]->SetXTitle("Difference of Channels");
      htimez[i]->SetYTitle("Number of Entries");
    }
   
    //re-size bins for 2D histograms
    //xsbin/=4;
    while(xsbin>512){
      xsbin/=2;
    }
    //xdbin/=4;
    for(int i = 0; i < 12; i++){//2D diagnostic histograms
      subtitle="";
      if(i<6)
	subtitle=", A1";
      else
	subtitle=", A2";
      switch(i%3)
	{
	case 0:
	  subtitle+="B";
	  break;
	case 1:
	  subtitle+="M";
	  break;
	case 2:
	  subtitle+="T";
	  break;
	default:
	  break;
	}
      subsubtitle=" vs ";
      if(i<3)
	subsubtitle+="X1";
      else
	if(i<6)
	  subsubtitle+="Y1";
	else
	  if(i<9)
	    subsubtitle+="X2";
	  else
	    if(i<12)
	      subsubtitle+="Y2";
     
      if(((int)(i/3))%2){//Y signals
	xdmin=-1000;
	xdmax=1000;
	pmin=pymin;
	pmax=pymax;
      }
      else{//X signals
	xdmin=-800;
	xdmax=2100;
	pmin=pxmin;
	pmax=pxmax;
      }
    
      xdbin=(Int_t)(xdmax-xdmin);
      while(xdbin>512){
	xdbin/=2;
      }         
      
      while(abin>512){
	abin/=2;
      }        

      sprintf(name,"hasum%i",i);  
      title="Anode vs Sum";
      if(free_range)
	hasum[i] = new TH2F(name,title+subtitle+subsubtitle,xsbin,1e5,-1e5,abin,1e5,-1e5);
      else
	hasum[i] = new TH2F(name,title+subtitle+subsubtitle,abin,amin,amax,xsbin,xsmin,xsmax);
      hasum[i]->SetXTitle("Anode Channel");
      hasum[i]->SetYTitle("Sum of Cathode Channels");
      
      if(bdiag) {
	sprintf(name,"hasumn%i",i);  
	title="Anode vs Sum, Sum Gain-matched";
	hasumn[i] = new TH2F(name,title+subtitle+subsubtitle,abin,amin,amax,xsbin,xsmin,xsmax);
	hasumn[i]->SetXTitle("Anode Channel");
	hasumn[i]->SetYTitle("Sum of Cathode Channels");


	sprintf(name,"hasumnn%i",i);  
	title="Anode vs Sum, Fullly Gain-matched";
	hasumnn[i] = new TH2F(name,title+subtitle+subsubtitle,abin*2,amin,amax,xsbin,xsmin,xsmax);
	hasumnn[i]->SetXTitle("Anode Channel");
	hasumnn[i]->SetYTitle("Sum of Cathode Channels");
      }

      sprintf(name,"haxf%i",i);  
      title="Anode vs Sum";
      haxf[i] = new TH2F(name,title+subtitle+subsubtitle+" X_{far}",xbin2D,xmin,xmax,abin,amin,amax);
      haxf[i]->SetXTitle("Cathode Channel");
      haxf[i]->SetYTitle("Anode Channel");

      sprintf(name,"haxn%i",i);  
      title="Anode vs Sum";
      haxn[i] = new TH2F(name,title+subtitle+subsubtitle+" X_{near}",xbin2D,xmin,xmax,abin,amin,amax);
      haxn[i]->SetXTitle("Cathode Channel");
      haxn[i]->SetYTitle("Anode Channel");

      sprintf(name,"hadiff%i",i);   
      title="Anode vs Difference";
      if(free_range)
	hadiff[i] = new TH2F(name,title+subtitle+subsubtitle,xdbin,1e5,-1e5,abin,1e5,-1e5);
      else
	hadiff[i] = new TH2F(name,title+subtitle+subsubtitle,xdbin,xdmin,xdmax,abin,amin,amax);
      hadiff[i]->SetXTitle("Difference of Cathode Channels");
      hadiff[i]->SetYTitle("Anode Channel");

      xbin=xbin2D;

      sprintf(name,"hax%i",i);   
      title="Anode vs Position";
      hax[i] = new TH2F(name,title+subtitle+subsubtitle,xbin*2,xxmin,xxmax,abin,amin,amax);
      hax[i]->SetXTitle("Position (relative)");
      hax[i]->SetYTitle("Anode Channel");

      sprintf(name,"haxc%i",i);   
      title="Anode vs Position, Calibrated";
      if(match_sim)
	haxc[i] = new TH2F(name,title+subtitle+subsubtitle,cal_bin*2,pmin,pmax,abin,amin,amax);
      else
	haxc[i] = new TH2F(name,title+subtitle+subsubtitle,xbin,pmin,pmax,abin,amin,amax);
      haxc[i]->SetXTitle("Position (mm)");
      haxc[i]->SetYTitle("Anode Channel");

      tmin=adiffmin;//-300.0;
      tmax=adiffmax;//-tmin;
      tbin=(int)(tmax-tmin);
      
      subtitle=", A2";
      switch(i%3)
	{
	case 0:
	  subtitle+="B-A1B";
	  break;
	case 1:
	  subtitle+="M-A1M";
	  break;
	case 2:
	  subtitle+="T-A1T";
	  break;
	default:
	  break;
	}

      sprintf(name,"htx%i",i);   
      title="Time vs Position";
      htx[i] = new TH2F(name,title+subtitle+subsubtitle,xbin,xxmin,xxmax,tbin,tmin,tmax);
      htx[i]->SetXTitle("Position (relative)");
      htx[i]->SetYTitle("Difference of Anode Channels (Time)");

      sprintf(name,"htxc%i",i);   
      title="Time vs Position, Calibrated";
      if(match_sim)
	htxc[i] = new TH2F(name,title+subtitle+subsubtitle,cal_bin,pmin,pmax,tbin,tmin,tmax);
      else
	htxc[i] = new TH2F(name,title+subtitle+subsubtitle,xbin,pmin,pmax,tbin,tmin,tmax);
      htxc[i]->SetXTitle("Position (mm)");
      htxc[i]->SetYTitle("Difference of Anode Channles (Time)");

      subsubtitle+=" sum";

      sprintf(name,"htsum%i",i);   
      title="Time vs Sum";
      htsum[i] = new TH2F(name,title+subtitle+subsubtitle,xsbin,xsmin,xsmax,tbin,tmin,tmax);
      htsum[i]->SetXTitle("Sum of Cathode Channels");
      htsum[i]->SetYTitle("Difference of Anode Channels (Time)");
    }
    
    for(int i = 0; i < 4; i++){//1D position plots
      if(i<2) {
	if(i==0) title="X";
	if(i==1) title="Y";
	subtitle=" Position Difference";
	sprintf(name,"hxdiff%i",i);
	hxdiff[i] = new TH1F(name,title+subtitle,xbin,-0.006,0.006);
	hxdiff[i]->SetXTitle("Position Difference (relative)");
	hxdiff[i]->SetYTitle("Number of Entries");
      }      
      
      if(i==0)
	sprintf(pretitle,"X1");
      if(i==1)
	sprintf(pretitle,"Y1");
      if(i==2)
	sprintf(pretitle,"X2");
      if(i==3)
	sprintf(pretitle,"Y2"); 
      title=pretitle;

      subtitle=" Position";
      sprintf(name,"hx%i",i);
      if(free_range)
	hx[i] = new TH1F(name,title+subtitle,xbin,1e5,-1e5);
      else
	hx[i] = new TH1F(name,title+subtitle,xbin,xxmin,xxmax);
      hx[i]->SetXTitle("Position (relative)");
      hx[i]->SetYTitle("Number of Entries");

      sprintf(name,"hsx%i",i);   
      subtitle=" Sum vs Position";
      hsx[i] = new TH2F(name,title+subtitle,xbin,xxmin,xxmax,xsbin,xsmin,xsmax);
      hsx[i]->SetXTitle("Position (relative)");
      hsx[i]->SetYTitle("Sum of Cathode Channels");

      /* sprintf(name,"hsa%i",i);  
	 subtitle=" Sum vs Anode";
	 if(free_range)
	 hsa[i] = new TH2F(name,title+subtitle,abin,1e5,-1e5,xsbin,1e5,-1e5);
	 else
	 hsa[i] = new TH2F(name,title+subtitle,abin,amin,amax,xsbin,xsmin,xsmax);
	 hsa[i]->SetXTitle("Anode Channel");
	 hsa[i]->SetYTitle("Sum of Cathode Channels");*/


      if(i%2){//Y signals
	pmin=pymin;
	pmax=pymax;
      }
      else{//X signals
	pmin=pxmin;
	pmax=pxmax;
      }

      subtitle=" Position, Calibrated";
      sprintf(name,"hxc%i",i);
      if(free_range)
	hxc[i] = new TH1F(name,title+subtitle,xbin,1e5,-1e5);
      else
	if(match_sim)
	  hxc[i] = new TH1F(name,title+subtitle,cal_bin,pmin,pmax);
	else
	  hxc[i] = new TH1F(name,title+subtitle,xbin,pmin,pmax);
      hxc[i]->SetXTitle("Position (mm)");
      hxc[i]->SetYTitle("Number of Entries");

      subtitle=" Position, Calibrated & Gated";
      sprintf(name,"hxcg%i",i);
      hxcg[i] = new TH1F(name,title+subtitle,cal_bin,pmin,pmax);
      hxcg[i]->SetXTitle("Position (mm)");
      hxcg[i]->SetYTitle("Number of Entries");

      if(i<2) {
	if(i==0) title="X";
	if(i==1) title="Y";
	subtitle=" Position Difference, Calibrated";
	sprintf(name,"hxcdiff%i",i);
	hxcdiff[i] = new TH1F(name,title+subtitle,xbin,-10,10);
	hxcdiff[i]->SetXTitle("Position Difference (mm)");
	hxcdiff[i]->SetYTitle("Number of Entries");

	subtitle=" ";
	sprintf(name,"hrho%i",i);
	hrho[i] = new TH1F(name,title+subtitle,xbin,-50,200);
	hrho[i]->SetXTitle("Position Difference (mm)");
	hrho[i]->SetYTitle("Number of Entries");
	
	subtitle=" Positions, Det 1 vs Det 2";
	sprintf(name,"hxx%i",i);
	if(i==0) hxx[i] = new TH2F(name,title+subtitle,cal_bin,xxmin,xxmax,
				   cal_bin,xxmin,xxmax);
	if(i==1) hxx[i] = new TH2F(name,title+subtitle,cal_bin,xymin,xymax,
				   cal_bin,xymin,xymax);
	hxx[i]->SetXTitle("Det 2 Position (mm)");
	hxx[i]->SetYTitle("Det 1 Position (mm)");

	subtitle=" Positions, Det 1 vs Det 2";
	sprintf(name,"hxxc%i",i);
	if(i==0) hxxc[i] = new TH2F(name,title+subtitle,cal_bin,pxmin,pxmax,
				    cal_bin,pxmin,pxmax);
	if(i==1) hxxc[i] = new TH2F(name,title+subtitle,cal_bin,pymin,pymax,
				    cal_bin,pymin,pymax);
	hxxc[i]->SetXTitle("Det 2 Position (mm)");
	hxxc[i]->SetYTitle("Det 1 Position (mm)");
      }

      sprintf(name,"hsxc%i",i);   
      title="Sum vs Position, Calibrated";
      if(match_sim)
	hsxc[i] = new TH2F(name,title+subtitle,cal_bin,pmin,pmax,xsbin,xsmin,xsmax);
      else
	hsxc[i] = new TH2F(name,title+subtitle,xbin,pmin,pmax,xsbin,xsmin,xsmax);
      hsxc[i]->SetXTitle("Position (mm)");
      hsxc[i]->SetYTitle("Sum of Cathode Channels");
      
      title=pretitle;

      subtitle+=" Angle to Detector 2";
      sprintf(name,"hang%i",i);   
      hang[i] = new TH1F(name,title+subtitle,xbin,-20,20);
      hang[i]->SetXTitle("Position (relative)");
      hang[i]->SetYTitle("Number of Entries");

      subtitle=" Angle to Target";
      sprintf(name,"hangt%i",i);   
      hangt[i] = new TH1F(name,title+subtitle,xbin,-7,7);
      hangt[i]->SetXTitle("Angle to Target (deg)");
      hangt[i]->SetYTitle("Number of Entries");
    }   

    for(int i=0; i<2; i++){//2D position plots 
      xbin2D=dbin/8;
      sprintf(name,"hhit%i",i);
      title = "Y vs X ";
      title+=i+1;
          
      if(free_range)
	hhit[i] = new TH2F(name,title,xbin2D,1e5,-1e5,xbin2D,1e5,-1e5);
      else
	hhit[i] = new TH2F(name,title,xbin2D,xxmin,xxmax,xbin2D,xymin,xymax);
      hhit[i]->SetXTitle("Horizontal Position (relative)");
      hhit[i]->SetYTitle("Vertical Position (relative)");
      
      sprintf(name,"hss%i",i);   
      subtitle=", Sums zoomed";
      hss[i] = new TH2F(name,title+subtitle,xsbin,xsmin,xsmax,xsbin,xsmin,xsmax);
      hss[i]->SetXTitle("Sum of X-Cathode Channels");
      hss[i]->SetYTitle("Sum of Y-Cathode Channels");

      sprintf(name,"hhitc%i",i);      
      title+=" Calibrated";
      if(free_range)
	hhitc[i] = new TH2F(name,title,xbin2D,1e5,-1e5,xbin2D,1e5,-1e5);
      else
	if(match_sim)
	  hhitc[i] = new TH2F(name,title,cal_bin*2,pxmin,pxmax,cal_bin*4,pymin,pymax);
	else
	  hhitc[i] = new TH2F(name,title,xbin2D,pxmin,pxmax,xbin2D,pymin,pymax);	  

      hhitc[i]->SetXTitle("Horizontal Position (mm)");
      hhitc[i]->SetYTitle("Vertical Position (mm)");
      //sprintf(name,"hhitcg%i",i);      
      //title+=" Gated";
      //if(match_sim)
      //hhitcg[i] = new TH2F(name,title,cal_bin,pxmin,pxmax,cal_bin,pymin,pymax);
      //else
      //hhitcg[i] = new TH2F(name,title,xbin2D,pxmin,pxmax,xbin2D,pymin,pymax);
      //hhitcg[i]->SetXTitle("Horizontal Position (mm)");
      //hhitcg[i]->SetYTitle("Vertical Position (mm)");  
    }
    if(bdiag) {
      for(int i=0; i<planes; i++) {
	sprintf(name,"hsource%i",i);      
	title="Rays traced to Z=";
	title+=(int)zp[i];
	
	if(zp[i]==0.0)
	  subtitle=" source (origin)";
	else
	  if(zp[i]==637.0)
	    subtitle=" mask";
	  else
	    if(zp[i]<0.0)
	      subtitle=" behind source (origin)";
	    else
	      subtitle="";
	xsource=160;
	hsource[i] = new TH2F(name,title+subtitle,xbin2D,-xsource,xsource,
			      xbin2D,-xsource,xsource);
      }
      //3D position plots
      pzbin=((Int_t)((pzmax-pzmin)/(Zp[0]-Zp[1])))*2;
      htrace = new TH3F("htrace","space",xbin3D,xxmin,xxmax,xbin3D,xxmin,xxmax,4,0,3);
      htracec = new TH3F("htracec","space",xbin3D,pxmin,pxmax,xbin3D,pymin,pymax,pzbin,pzmin,pzmax);
      hzx = new TH2F("hzx","Z vs X",cal_bin,pxmin,pxmax,pzbin,pzmin,pzmax);
      hzy = new TH2F("hzy","Z vs Y",cal_bin,pymin,pymax,pzbin,pzmin,pzmax);
    }
    // Begin of run calls...

#ifdef  USE_V1190
    //v1190_histos->BeginRun(transition,run,time);
#endif 

  }//end BeginRun

  bool ProcessMidasEvent(TDataContainer& dataContainer) {
    TV1190Data *data = dataContainer.GetEventData<TV1190Data>("EMMT");
    if(!data) return false;
    std::vector<TDCMeasurement> measurements = data->GetMeasurements();
    int chan = 0; //channel number
    int datum[Nchan][Ncounts] = {{0}}; //data word
    std::vector<int> counts(Nchan,0); //number of counts (multi-hit)
    
    //flags
    Bool_t bcoinc[Nchan/2]={0};//signal-pair coincidences
    //Bool_t bac[6]; //anode peak width gate
    //    Bool_t bsum[4]={0}; //cathode sum peak width gate
    Bool_t baacoinc[3][3]={{0}}; //andode-anode coincidence
    Int_t mhit=0;
    //Int_t aa_delta=0; //differene in position between anodes
    //Bool_t bsacoinc[2][3]={{0}}; //cathode sum anode coincidence
    //Int_t sa_sum=0; //sum of coincidences btween cathodes and anodes

    //anode variables
    Float_t a[6][2]={{0}};//anode signals: raw [0] and cal [1]
    Float_t atime=0;
    // Float_t arange=9956-9644;//anode range
    //Float_t awide=arange/11;//selects 6/66 of the anode range
    //Float_t amean=9891;//anode mean

    //cathode and position variables
    //...for pair of signals
    Float_t xf[4][3]={{0}};
    Float_t xn[4][3]={{0}};
    Float_t  x_sum[4][3]={{0}};
    Float_t x_diff[4][3]={{0}};
    Float_t x[4][2]={{0}};
    
    //derived variables
    Float_t ang[4]={};//relative angle
    Float_t angt[4]={};//central angle, to target
    Float_t rho=0;
    /*Float_t X0=0, Y0=0;
      Float_t xwidth=0, ywidth=0;
      Float_t xcenter=0, ycenter=0;*/

    //--- Unpack data
    for(unsigned int i = 0; i < measurements.size(); i++) { // loop over measurements
      chan = measurements[i].GetChannel();
      if((GetCurrentRunNumber()>480)&&(GetCurrentRunNumber()<612)) {//switch inputs if second PGAC run
	if((chan>5)&&(chan<10)) {
	  chan=15-chan;
	}
      }
      counts[chan]++; //"counts" uses natural counting; i.e. 0 is 0 counts, 1 is 1 count...
      if(counts[chan]>maxcount){
	maxcount=counts[chan];
	printf(" maxcount = %d\n",maxcount);
      }
      datum[chan][counts[chan]] = measurements[i].GetMeasurement();
      
      if(bdiag){
	//if(counts[chan]==1)
	hdata[chan]->Fill(datum[chan][counts[chan]]);// raw data histograms
      }
      //printf("  hdata[%2d][%d]=%5d\n",chan,counts[chan],datum[chan][counts[chan]]);
      if(counts[chan]>1) {//multi-hit diagnostic
	printf("   multihit!\n   evt = %d, hdata[%2d] = \n",nevent,chan);
	for(int j = 1; j < counts[chan]+1; j++) {
	  printf("    %d %5d\n",j,datum[chan][j]);
	  if(j>1)
	    printf("    delta = %d\n",datum[chan][j]-datum[chan][j-1]);
	  if(datum[chan][j]<datum[chan][j-1])
	    printf("     out of order!\n");
	}
	
      }
    }//end unpack data
    
    if(bdiag) {//counts histograms   
      for(int i = 0; i < Nchan; i++) {// loop over channels    
	hcounts[i]->Fill(counts[i]); //counts per multi-hit
	for(int j = 0; j < counts[i]; j++) {// loop over counts
	  hcounts[Nchan]->Fill(i); //non-zero counts per channel
	}
      }
    }
    
    //Int_t j=1;
    //maxcount=1;
    //--- Sort data and apply calibration
    for(Int_t j=1;j<=maxcount;j++)
    {//loop over counts
      for(int i=0;i<Nchan/2;i++) {//loop over pairs of signals
	if(i<3) {//anode data
	  a[i][0]=datum[i][j];
	  a[i+3][0]=datum[i+3][j];
	  if(docal) {
	    //Step 1: Gain-match time signals.
	    //the offset is included separately to easily shift the mean to zero
	    a[i][1]=a[i][0]/hxxcal[i][1];	 
	    a[i][1]-=hxxcal[i][0];
	  }
	  else 
	    a[i][1]=a[i][0];
	  a[i+3][1]=a[i+3][0];//note: only A1 is modified
	}
	else {//cathode data
	  if(i%2) {// x-signals (the x-signals are labeled from downstream; I am "flipping" them 
	    //        to keep the relative positions defined from upstream: XR is "left" and XL 
	    //        is "right")
	    xn[i-3][0]=datum[(2*i)+1][j]; // "near" signal; labeled right, beam left
	    xf[i-3][0]=datum[(2*i)+0][j]; // "far" signal; labeled left, beam right
	  }
	  else    {// y-signals
	    xn[i-3][0]=datum[(2*i)+0][j]; // "near" signal; bottom
	    xf[i-3][0]=datum[(2*i)+1][j]; // "far" signal; top
	  }	
	  if(docal) {
	    //Step 2: gain-match position signals
	    if(hxxcal[i][1]<-1.) {//this if-statement is used so that the expansion coefficient is always >1
	      xn[i-3][1]=xn[i-3][0]*(-hxxcal[i][1]);
	      xf[i-3][1]=xf[i-3][0];
	    }
	    else{
	      xn[i-3][1]=xn[i-3][0]; 
	      xf[i-3][1]=xf[i-3][0]*(-1./hxxcal[i][1]);
	    }
	    //Step 3: gain-match sum of position signals to anode signal
	    //  i.e., remove the anode-dependance of the sum of the position signals
	    //  i.e., make sum of position signals constant using anode signal
	    for(int k=0;k<3;k++) {
	      //use slope of fit to remove anode-dependance, then center
	      if(a[k+3*(int)((i-3)/2)][1]>0) {
		xn[i-3][2]=xn[i-3][1]-(a[k+3*(int)((i-3)/2)][1]*hasumcal[((i-3)*3)+k][1])/2;
		xn[i-3][2]+=hasumcal[((i-3)*3)+k][2]/2;
		//xn[i-3][2]+=100*k;
		xf[i-3][2]=xf[i-3][1]-(a[k+3*(int)((i-3)/2)][1]*hasumcal[((i-3)*3)+k][1])/2;
		xf[i-3][2]+=hasumcal[((i-3)*3)+k][2]/2;
		//xf[i-3][2]+=100*k;
	      }
	    }
	  }
	  else
	    {//no calibration applied
	      xn[i-3][1]=xn[i-3][0]; 
	      xf[i-3][1]=xf[i-3][0];
	      xn[i-3][2]=xn[i-3][0]; 
	      xf[i-3][2]=xf[i-3][0]; 
	    }
	}
      }//calibration done

      //--- Fill histograms
      for(int i=0;i<Nchan/2;i++) {//loop over pairs of signals
	if(i<3) {//anode data
	  atime=a[i+3][1]-a[i][1];
	  mhit=0;
	  for(int k=0;k<3;k++) {
	    if((a[i][1]>0)&&(a[k+3][1]>0))
	      {
		baacoinc[k][i]=kTRUE;
		haacoinc->Fill(3*i+k);
		mhit++;
		htime [(3*i)+k]->Fill(a[k+3][1]-a[i][1]);
		htimez[(3*i)+k]->Fill(a[k+3][1]-a[i][1]);
		haa   [(3*i)+k]->Fill(a[k+3][1],a[i][1]);
		haaz  [(3*i)+k]->Fill(a[k+3][1],a[i][1]);
	      }
	  }
	  hmhit->Fill(mhit);
    	  
	  if(a[i][0]&&a[i+3][0]) {//only fill histogram if both signals have a hit
	    bcoinc[i]=true;
	    //hsum[i]->Fill(a[i+3][1]+a[i][1]);//does this have any physical meaning?
	    //hsumz[i]->Fill(a[i+3][1]+a[i][1]);
	    if(bdiag){
	      hcounts[Nchan+1]->Fill(i);//coincidences
	      htev[i]->Fill(nevent,atime);
	    }
	  }
	}//end anode loop (i<3)
	else{//cathode data
	  if(xf[i-3][0]&&xn[i-3][0]) {//only fill histogram if both position signals have a hit
	    bcoinc[i]=true;
	    if(bdiag)
	      hcounts[Nchan+1]->Fill(i);
	    if(bdiag) {
	      hxfxn[i]->Fill(xn[i-3][0],xf[i-3][0]);
	      hxfxnz[i]->Fill(xn[i-3][0],xf[i-3][0]);
	      hxfxnzn[i]->Fill(xn[i-3][1],xf[i-3][1]);
	      hxfxnznn[i]->Fill(xn[i-3][2],xf[i-3][2]);
	    }
	    else
	      hxfxnz[i]->Fill(xn[i-3][2],xf[i-3][2]);
	    for(Int_t k=0;k<3;k++) {
	      x_sum [i-3][k]=xf[i-3][k]+xn[i-3][k];
	      x_diff[i-3][k]=xn[i-3][k]-xf[i-3][k];//with xf, xn, and x defined correctly, x_diff is (positively) proportional to position
	    }
   
	    hsum[i] ->Fill(x_sum[i-3][2]);
	    hsumz[i]->Fill(x_sum[i-3][2]);
	    hdiff[i] ->Fill(x_diff[i-3][2]);
	    hdiffz[i]->Fill(x_diff[i-3][2]);
	    
	    x[i-3][0]=(1/2.)*(1+((x_diff[i-3][2])/(x_sum[i-3][2]))); //Position on detector with XN@x=0 and XF@x=1.
	    hx[i-3]->Fill(x[i-3][0]);
	    
	    for(int k=0;k<3;k++) {
	      if(a[k+3*(int)((i-3)/2)][1]>0)
		x[i-3][1]=(x[i-3][0]*hxcal[((i-3)*3)+k][0])+hxcal[((i-3)*3)+k][1];//inverse of usual peakfit scaling
	    }
	    hxc[i-3]->Fill(x[i-3][1]);
	    
	    hsdiff[i]->Fill(x_diff[i-3][2],x_sum[i-3][2]);
	    hsx[i-3]->Fill(x[i-3][0],x_sum[i-3][2]);
	    hsxc[i-3]->Fill(x[i-3][1],x_sum[i-3][2]);
	 
	    angt[i-3]=TMath::RadToDeg()*TMath::ATan((x[i-3][1]-offset_cal[(i-3)%2])/Zp[(int)((i-3)/2)]);
	    hangt[i-3]->Fill(angt[i-3]);
	   
	    //--- Plots with anode data ----//
	    for(int k=0;k<3;k++) {//loop over anode segments, 0=B,1=M,2=T	   
	      //time needs to be re-defined here, otherwise the i=3 "top" time would be used
	      atime=a[k+3][1]-a[k][1];//A2-A1
	      if(i<5) {//cathodes of det 1
		if(a[k][1]>0) {
		  if(bdiag){
		    hasum  [((i-3)*3)+k]->Fill(a[k][1],x_sum[i-3][0]);
		    hasumn [((i-3)*3)+k]->Fill(a[k][1],x_sum[i-3][1]);
		    hasumnn[((i-3)*3)+k]->Fill(a[k][1],x_sum[i-3][2]);
		  }
		  else
		    hasum[((i-3)*3)+k]->Fill(a[k][1],x_sum[i-3][2]);
		  haxf   [((i-3)*3)+k]->Fill(xf[i-3][2],     a[k][1]);
		  haxn   [((i-3)*3)+k]->Fill(xn[i-3][2],     a[k][1]);
		  hadiff [((i-3)*3)+k]->Fill(x_diff[i-3][2],    a[k][1]);
		  hax    [((i-3)*3)+k]->Fill(x[i-3][0], a[k][1]);
		  haxc   [((i-3)*3)+k]->Fill(x[i-3][1], a[k][1]);
		}
	      }
	      else {//cathodes of det 2
		if(a[k+3][1]>0) {
		  if(bdiag){
		    hasum  [((i-3)*3)+k]->Fill(a[k+3][1],x_sum[i-3][0]);
		    hasumn [((i-3)*3)+k]->Fill(a[k+3][1],x_sum[i-3][1]);
		    hasumnn[((i-3)*3)+k]->Fill(a[k+3][1],x_sum[i-3][2]);
		  }
		  else
		    hasum[((i-3)*3)+k]->Fill(a[k+3][1],x_sum[i-3][2]);		   
		  haxf   [((i-3)*3)+k]->Fill(xf[i-3][2],     a[k+3][1]);
		  haxn   [((i-3)*3)+k]->Fill(xn[i-3][2],     a[k+3][1]);
		  hadiff [((i-3)*3)+k]->Fill(x_diff[i-3][2],    a[k+3][1]);
		  hax    [((i-3)*3)+k]->Fill(x[i-3][0], a[k+3][1]);
		  haxc   [((i-3)*3)+k]->Fill(x[i-3][1], a[k+3][1]);
		}
	      }
	      if(bcoinc[k]) {
		htx  [((i-3)*3)+k]->Fill(x[i-3][0],         atime);
		htxc [((i-3)*3)+k]->Fill(x[i-3][1],atime);
		htsum[((i-3)*3)+k]->Fill(x_sum[i-3][2],atime);
	      }
	      //if(bac[k]&&bac[k+3])
	      //hxcg[i-3]->Fill(x[i-3][1]);
	    }//end anode segment loop
	    
	  }//end cathode conicidence
	}//end cathode data
      }//end signal-pair loop
      if(bcoinc[3]&&bcoinc[4]) {// x-y coincidence for det 1
	hhit[0]->Fill(x[0][0],x[1][0]); 
	hhitc[0]->Fill(x[0][1],x[1][1]); 
	//if(bsum[0]&&bsum[1])
	//hhitcg[0]->Fill(x[0][1],x[1][1]); 
      }
      if(bcoinc[5]&&bcoinc[6]) {// x-y coincidence for det 2
	hhit[1]->Fill(x[2][0],x[3][0]);      
	hhitc[1]->Fill(x[2][1],x[3][1]);      
	//if(bsum[2]&&bsum[3])
	// hhitcg[1]->Fill(x[0][1],x[1][1]); 
      }
      if(bcoinc[3]&&bcoinc[4]&&bcoinc[5]&&bcoinc[6]) {// x-y coincidence for det 1 & 2
	for(int i=0;i<2;i++) {
	  hxdiff[i]->Fill(x[i][0]-x[i+2][0]);
	  hxcdiff[i]->Fill(x[i][1]-x[i+2][1]);
	  rho=TMath::Sqrt(TMath::Power(x[i][1]-offset_cal[0],2)
			  +TMath::Power(x[i+1][1]-offset_cal[1],2));
	  hrho[i]->Fill(rho);
	  hxx[i]->Fill(x[i+2][0],x[i][0]);
	  hxxc[i]->Fill(x[i+2][1],x[i][1]);
	  hss[i]->Fill(x_sum[i+2][2],x_sum[i][2]);
	  if(bdiag){
	    htrace->Fill(x[2*i][0],x[(2*i)+1][0],2*i);
	    htracec->Fill(x[2*i][1],x[(2*i)+1][1],Zp[i]);
	    hzx->Fill(x[2*i][1],Zp[i]);
	    hzy->Fill(x[(2*i)+1][1],Zp[i]);
	  }
	  ang[i]=TMath::RadToDeg()*TMath::ATan((x[i][1]-x[i+2][1])/(Zp[0]-Zp[1]));
	  hang[i]->Fill(ang[i]);
	 
	}
	/*	for(int i=0;i<planes;i++) { // loop over planes
		Z[4]=zp[i];
		X0=x[2][1]-(Z[2]-Z[4])*(x[0][1]-x[2][1])/(Z[0]-Z[2]);
		Y0=x[3][1]-(Z[3]-Z[4])*(x[1][1]-x[3][1])/(Z[1]-Z[3]);
		//X0=x[2][1]-(Z[2]-Z[4])*(x[0][1])/(Z[0]);
		//Y0=x[3][1]-(Z[3]-Z[4])*(x[1][1])/(Z[1]);
		if(bdiag)
		hsource[i]->Fill(X0,Y0);
		// hit source coincidence
		xcenter=0.262090;
		ycenter=17.4166;
		xwidth=20;
		ywidth=xwidth;
		if((fabs(X0-xcenter)<xwidth)&&(fabs(Y0-ycenter)<ywidth)){
		hhitcg[0]->Fill(x[0][1],x[1][1]); 
		hhitcg[1]->Fill(x[2][1],x[3][1]); 
		}
		}*/
      }
    }//end loop over counts
    nevent++;
    return true;
  }//end ProcessMidasEvent
}; 

int main(int argc, char *argv[])
{
  Analyzer::CreateSingleton<Analyzer>();
  return Analyzer::Get().ExecuteLoop(argc, argv);
}
