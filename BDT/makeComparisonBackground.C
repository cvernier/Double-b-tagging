#include "CMS_lumi.C"

void makeComparisonBackground(TString var , double low =0., double max =100., TString label){

gROOT->SetStyle("Plain");
gStyle->SetPadGridX(0);
gStyle->SetPadGridY(0);
gStyle->SetOptStat(0);


int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV
int iPos =11;



TCanvas * c = new TCanvas("c","c", 700, 700);
TLegend *leg = new TLegend(0.434465625,0.7021654,0.8365625,0.8603839,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.033);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
    leg->SetTextFont(42);

TFile * f = new TFile("qcd_forTraining.root");
f->cd();
Fjets->Draw(var+">>h"," abs(int(flavour))==5 && nbHadrons<2 && massPruned>=70 && massPruned<200");//ptPruned>300 & ptPruned<500");
h->SetLineColor(kBlue+2);
h->SetLineWidth(3);
h->GetXaxis()->SetTitle(label);
h->GetXaxis()->SetTitleFont(42);
h->GetXaxis()->SetLabelFont(42);
h->GetYaxis()->SetTitleFont(42);
h->GetYaxis()->SetLabelFont(42);
h->SetTitle("");
double Ymax = h->GetMaximum();
h->GetYaxis()->SetRangeUser(0.1, Ymax*3.);
h->DrawNormalized("");
int nbin = 8;// h->GetXaxis()->GetNbins();
std::stringstream ss;
ss << nbin;
TString strNb = ss.str();
//float low =  h->GetXaxis()->GetXmin();
//if( low<0) {low=0; }
std::stringstream ss1;
ss1 << low;
TString strL=(ss1.str());
//float max = h->GetXaxis()->GetXmax();
std::stringstream ss2;
ss2 << max;
TString strM= (ss2.str());
Fjets->Draw(var+">>h("+strNb+","+strL+","+strM+")"," abs(int(flavour))==5 && nbHadrons<2 && massPruned>=70 && massPruned<200");//ptPruned>300 & ptPruned<500");
h->SetLineColor(kBlue+2);
h->SetLineWidth(3);
h->GetXaxis()->SetTitle(label);
h->GetYaxis()->SetTitle("a.u.");
h->SetTitle("");
h->GetYaxis()->SetRangeUser(10., 10.*Ymax);
h->GetXaxis()->SetTitleFont(42);
h->GetXaxis()->SetLabelFont(42);
h->GetYaxis()->SetTitleFont(42);
h->GetYaxis()->SetLabelFont(42);
h->Rebin(1);
h->DrawNormalized("");

//std::cout<<nbin<< "  "<<low<<"  "<<max<<std::endl;
leg->AddEntry("h","QCD, single b","l");
TCanvas * d = new TCanvas("d","d");
d->cd();
Fjets->Draw(var+">>h2("+strNb+","+strL+","+strM+")"," abs(int(flavour))==5 && nbHadrons>1 && massPruned>=70 && massPruned<200");//ptPruned>300 & ptPruned<500");
h2->SetLineColor(kAzure+1);
h2->SetLineWidth(3);
h2->GetXaxis()->SetTitle(var);
h2->SetTitle("");
h2->Rebin(1);
//h2->GetYaxis()->SetRangeUser(0., Ymax*2.);
c->cd();
h2->DrawNormalized("same");
leg->AddEntry("h2","QCD, gluon splitting to b#bar{b}","l");
delete h2;



d->cd();
Fjets->Draw(var+">>h3("+strNb+","+strL+","+strM+")"," abs(int(flavour))!=5 && massPruned>=70 && massPruned<200");//ptPruned>300 & ptPruned<500");
h3->SetLineColor(kBlack);
h3->SetLineStyle(2);
h3->SetLineWidth(3);
h3->GetXaxis()->SetTitle(var);
h3->SetTitle("");
h3->Rebin(1);
//h2->GetYaxis()->SetRangeUser(0., Ymax*2.);
c->cd();
h3->DrawNormalized("same");
leg->AddEntry("h3","QCD, light flavour","l");
delete h3;




d->cd();
TFile * fsignal = new TFile("bulkGrav800-4000_forTraining.root");// bulkGrav800-3000_forTraining.root");//grav_ALL_forTraining.root");
fsignal->cd();
std::cout<<var+">>h4("+strNb+","+strL+","+strM+")"<<std::endl;
Fjets->Draw(var+">>h4("+strNb+","+strL+","+strM+")","abs(int(flavour))==5 && nbHadrons>1 && massPruned>=70 && massPruned<200");//ptPruned>300 & ptPruned<500");

h4->SetLineColor(kRed+1);
h4->SetLineWidth(3);
h4->GetXaxis()->SetTitle(var);
h4->SetTitle("");
h4->Rebin(1);
c->cd();
h4->DrawNormalized("same");

leg->AddEntry("h4","H(b#bar{b})","l");
delete h4;
delete d;



leg->Draw();
if(var=="tau2IVF/tau1IVF") var = "subIVF";
if(var=="tau2/tau1") var = "sub";
 CMS_lumi( c, iPeriod, iPos );


    TLatex * l1 = new TLatex();
    l1->SetTextAlign(13);
    l1->SetTextFont(42);
    l1->SetNDC();

    l1->SetTextSize(0.027);
    l1->DrawLatex(0.13,0.23, "AK8, p_{T} > 300 GeV");
    l1->DrawLatex(0.13,0.2, "70 < m < 200 GeV");


c->Print(("plots-Background/"+var+".png"));
c->SetLogy();
c->Print(("plots-Background/"+var+"_log.pdf"));
//delete c;
//delete strL;


}
