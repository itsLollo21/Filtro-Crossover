





double L = 1.01*0.001;
double RL= 150.4;
double RiL= 12.25;
double RhoL = RL+RiL;
double TauL = L/RhoL;

double C= 881*TMath::Power(10,-9);
double RC = 150.6;
double RiC = 15.53;
double RhoC = RC+RiC;
double TauC = RhoC*C;




Double_t classic_VL_Function(Double_t *x, Double_t *par){
    Float_t xx = x[0];
    Double_t val = (150.6*5)/(162.8*TMath::Sqrt(1+TMath::Power(2*M_PI*xx*par[0],2)));
    return val;
}


Double_t classic_VC_Function(Double_t *x, Double_t *par){
    Float_t xx = x[0];
    Double_t val = (150.4*5)/(165.9*TMath::Sqrt(1+TMath::Power(2*M_PI*xx*par[0],-2)));
    return val;
}



Double_t VL_Function(Double_t *x, Double_t *par){

Float_t xx = x[0];
Double_t val = par[2]*par[3]*TMath::Sqrt(1+1/(TMath::Power(2*M_PI*xx*par[1],2)))
/(par[4]*TMath::Sqrt(TMath::Power(1+par[0]/par[1]+par[6]*(par[4]+par[5])/(par[4]*par[5]),2)
+TMath::Power(2*M_PI*xx*par[0]-1/(2*M_PI*xx*par[1])+par[6]*(2*M_PI*xx*par[0]/par[5]-1/(par[4]*2*M_PI*xx*par[1])),2)));
return val;
}

Double_t VC_Function(Double_t *x, Double_t *par){

Float_t xx = x[0];
Double_t val = par[2]*par[3]*TMath::Sqrt(1+TMath::Sqrt(1+TMath::Power(2*M_PI*xx*par[0],2)))
/(par[5]*TMath::Sqrt(TMath::Power(1+par[0]/par[1]+par[6]*(par[4]+par[5])/(par[4]*par[5]),2)
+TMath::Power(2*M_PI*xx*par[0]-1/(2*M_PI*xx*par[1])+par[6]*(2*M_PI*xx*par[0]/par[5]-1/(par[4]*2*M_PI*xx*par[1])),2)));
return val;
}

Double_t Vin_Function(Double_t *x, Double_t *par){

Float_t xx = x[0];
Double_t val = par[2]*TMath::Sqrt((TMath::Power(1+par[0]/par[1],2)+TMath::Power(2*M_PI*xx*par[0]-1/(2*M_PI*xx*par[1]),2))
/(TMath::Power(1+par[0]/par[1]+par[5]*(par[3]+par[4])/(par[3]*par[4]),2)
+TMath::Power(2*M_PI*xx*par[0]-1/(2*M_PI*xx*par[1])+par[5]*(2*M_PI*xx*par[0]/par[4]-1/(par[3]*2*M_PI*xx*par[1])),2)));
return val;
}

Double_t VL_Tau_Function(Double_t *x, Double_t *par){
//tau L è par[1], tau C è par[0]
Float_t xx = x[0];
Double_t val = RL*5*TMath::Sqrt(1+TMath::Power(2*M_PI*xx*par[0],-2))
/(RhoL*TMath::Sqrt(TMath::Power(1+par[1]/par[0]+(50*(RhoL+RhoC)/(RhoL*RhoC)),2)+(2*M_PI*xx*par[1]-1
/(2*M_PI*xx*par[0])+50*(2*M_PI*xx*par[1]/RhoC-1/(2*M_PI*xx*par[0]*RhoL)))));
return val;
}

Double_t VC_Tau_Function(Double_t *x, Double_t *par){
//tau L è par[1], tau C è par[0]
Float_t xx = x[0];
Double_t val = RC*5*TMath::Sqrt(1+TMath::Power(2*M_PI*xx*par[1],2))
/(RhoC*TMath::Sqrt(TMath::Power(1+par[1]/par[0]+(50*(RhoL+RhoC)/(RhoL*RhoC)),2)+(2*M_PI*xx*par[1]-1
/(2*M_PI*xx*par[0])+50*(2*M_PI*xx*par[1]/RhoC-1/(2*M_PI*xx*par[0]*RhoL)))));
return val;
}


/*void Fit_VL(){
    TF1 *VL = new TF1("Fit_VL",classic_VL_Function);
VL->SetParameter(1,150.6);
VL->SetParameter(2,162.8);
VL->SetParameter(3,5);

VL->SetParNames("RL","V","Rho L","Tau L");
//Double_t TauL= 1.01*TMath::Power(10,-3)/(12.25+150.6);
//VL->SetParameter(4,ris);
VL->Draw();
}*/

void small_sweep2(){
    
TGraphErrors *j1 = new TGraphErrors("small_sweep(Vin).dat","%lg%lg%lg");
TGraphErrors *j2 = new TGraphErrors("small_sweep(VL).dat","%lg%lg%lg");                  //TGraphErrors * graph =new TGraphErrors("grafico.dat","%lg %lg %*lg %lg"); 
TGraphErrors *j3 = new TGraphErrors("small_sweep(VC).dat","%lg%lg%lg");
//TF1 *VL = new TF1("Fit di VL","classic_VL_Function");

TCanvas *my_c = new TCanvas("small_swee","Filtro CrossOver");
//TF1 *f = new TF1("Vin fit","Vin_Function");
//TF1 *f = new TF1("linear","x*[0]+[1]");
j1->SetTitle("Small Sweep");
j1->GetXaxis()->SetTitle("Frequency (Hz)");
j1->GetYaxis()->SetTitle("Amplitude (V)");
j1->GetYaxis()->SetRangeUser(2.65,3.15);
j1->GetXaxis()->SetRangeUser(6300,6600);
j2->GetYaxis()->SetRangeUser(2.76,2.78);
j2->GetXaxis()->SetRangeUser(6300,6600);

j1->SetLineColor(kBlue);

j2->SetLineColor(kMagenta);

j3->SetLineColor(kSpring);
//j1->Fit(f);
//j1->Draw("APE"); <---- commenta per restringere l'asse y
//TCanvas *my_can = new TCanvas("small_sw","Filtro CrossOver");
//j2->Fit(f,"L");
TF1 *VL = new TF1("Fit_VL",VL_Tau_Function,6300,6600, 2 ); //numero parametri 2 per le funzioni generali
TF1 *VC = new TF1("Fit_VC",VC_Tau_Function,6300,6600, 2 );
//TF1 *VL = (TF1 *)gROOT->GetFunction("Fit_VL");







VL->SetParameter(0,TauC); //<----commenta per fit classico
VL->SetParameter(1,TauL);   //0 per la funzione classica, 1 per quella generale------------------
//VL->Draw();
j2->Fit(VL,"R");
double TL1= VL->GetParameter(1); //idem
double TLer1 = VL->GetParError(1); //idem
double TC1= VL->GetParameter(0); //<---- commenta per f classica
double TCer1 = VL->GetParError(0); //<--- commenta per f classica
cout<< "Tau L è:" << TL1 << "+-"<<TLer1 <<endl;
j2->Draw("APE");

VC->SetParameter(0,TauC);
VC->SetParameter(1,TauL); //<---- commenta per fit classico

j3->Fit(VC,"R");
double TL2= VC->GetParameter(1); //<------ commenta per f classico
double TLer2 = VC->GetParError(1); //<---- commenta per f classico
double TC2= VC->GetParameter(0);
double TCer2 = VC->GetParError(0);
cout<< "Tau C è:" << TC2 << "+-"<<TCer2 <<endl;
j3->Draw("same");

cout<<"--------------------------------------------------------------------------------------"<<endl;
cout<< "Il Chi quadro ridotto è:" << VC->GetChisquare()/VC->GetNDF()<<endl;

double TauC_pesato = (TC1*TMath::Power(TCer1,-2)+TC2*TMath::Power(TCer2,-2))/(TMath::Power(TCer1,-2)+TMath::Power(TCer2,-2));
double TauL_pesato = (TL1*TMath::Power(TLer1,-2)+TL2*TMath::Power(TLer2,-2))/(TMath::Power(TLer1,-2)+TMath::Power(TLer2,-2));
double TauL_pesatoER = TMath::Power(TMath::Sqrt(TMath::Power(TLer1,-2)+TMath::Power(TLer2,-2)),-1);
double TauC_pesatoER = TMath::Power(TMath::Sqrt(TMath::Power(TCer1,-2)+TMath::Power(TCer2,-2)),-1);
cout<<"Tau L pesato è :" << TauL_pesato << "+-" << TauL_pesatoER << endl;
cout<<"Tau C pesato è :" << TauC_pesato <<"+-" << TauC_pesatoER << endl;


    
}