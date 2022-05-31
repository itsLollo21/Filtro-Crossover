





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


void small_sweep2(){
    
TGraphErrors *j1 = new TGraphErrors("small_sweep(Vin).dat","%lg%lg%lg");
TGraphErrors *j2 = new TGraphErrors("small_sweep(VL).dat","%lg%lg%lg");                 
TGraphErrors *j3 = new TGraphErrors("small_sweep(VC).dat","%lg%lg%lg");


TCanvas *my_c = new TCanvas("small_swee","Filtro CrossOver");

TF1 *VL = new TF1("Fit_VL",VL_Tau_Function,6300,6600, 2 ); //numero parametri 2 per le funzioni generali
TF1 *VC = new TF1("Fit_VC",VC_Tau_Function,6300,6600, 2 );


VL->SetParameter(0,TauC);
VL->SetParameter(1,TauL);   

j2->Fit(VL,"R");
double TL1= VL->GetParameter(1); 
double TLer1 = VL->GetParError(1); 
double TC1= VL->GetParameter(0);
double TCer1 = VL->GetParError(0); 
cout<< "Tau L è:" << TL1 << "+-"<<TLer1 <<endl;
j2->Draw("APE");

VC->SetParameter(0,TauC);
VC->SetParameter(1,TauL); 

j3->Fit(VC,"R");
double TL2= VC->GetParameter(1); /
double TLer2 = VC->GetParError(1); 
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
