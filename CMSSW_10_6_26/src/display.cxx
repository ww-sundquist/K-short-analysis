#include <math.h>
#include <TGraph.h>

double exponential(double *t, double *params) { //define generic exponential decay function A*exp(-t/tau)
   double tau = params[0];
   double A = params[1]/tau; //params[1] will be assigned as the total number of enteries. This statement is in effect our normalization.
   double lem = -(*t) / tau;
   double res = A*std::exp(lem);
   return res;
}

double gaussian(double *m, double *params) {
   double mean = params[0];
   double stdev = params[1];
   double num = params[2];
   double lem2 = (-1.0/2)*(pow((*m) - mean,2)/pow(stdev,2));
   double res2 = (num/pow(2*3.14159265358979262*stdev,1/2))*std::exp(lem2)/(stdev); //*pow(2*3.14159265358979262,1/2));
   return res2;
}

void display() {

   //K-short display

   TFile *ksFile = new TFile("kshortAnalyze.root");
   TH1D *h_kshort_mass = (TH1D*)ksFile->Get("kshort/h_kshort_mass");
   TCanvas *c1 = new TCanvas("c1","c1",800,800);

   Double_t ksNumEntries = h_kshort_mass->GetEntries();
   double intlKSMassParams[3] = {0.5,0.1,ksNumEntries};
   TF1 *ks_gaus = new TF1("ks_gaus",gaussian,0.49,0.51,3);
   ks_gaus->SetParameters(intlKSMassParams);
   ks_gaus->SetRange(0.485,0.515);
   h_kshort_mass->Fit(ks_gaus, "R");

   double ks_max = h_kshort_mass->GetMaximum();

   double ks_meanVal = ks_gaus->GetParameter(0);
   TLine *ks_mean = new TLine(ks_meanVal, 0, ks_meanVal, ks_max);
   ks_mean->SetLineColor(kRed);
   ks_mean->SetLineStyle(kDashed);

   double ks_sigma = ks_gaus->GetParameter(1);
   double ks_5sigma_upper = ks_meanVal + 5.0*ks_sigma;
   std::cout<<"Upper limit: "<<ks_5sigma_upper<<std::endl;
   TLine *ks_5s_u = new TLine(ks_5sigma_upper, 0, ks_5sigma_upper, ks_max);
   ks_5s_u->SetLineColor(kBlue);
   ks_5s_u->SetLineStyle(kDashed);
   double ks_5sigma_lower = ks_meanVal - 5.0*ks_sigma;
   std::cout<<"Lower limit: "<<ks_5sigma_lower<<std::endl;
   TLine *ks_5s_l = new TLine(ks_5sigma_lower, 0, ks_5sigma_lower, ks_max);
   ks_5s_l->SetLineColor(kBlue);
   ks_5s_l->SetLineStyle(kDashed);

   double ks_cut_l = ks_meanVal - 0.5*ks_sigma;
   double ks_cut_u = ks_meanVal + 0.5*ks_sigma;
   std::cout<<"Lower cut at "<<ks_cut_l<<std::endl;
   std::cout<<"Upper cut at "<<ks_cut_u<<std::endl;

   h_kshort_mass->Draw();
   ks_gaus->Draw("SAME");
   ks_mean->Draw("SAME");
   ks_5s_u->Draw("SAME");
   ks_5s_l->Draw("SAME");

   h_kshort_mass->GetXaxis()->SetTitle("Mass [GeV]");
   h_kshort_mass->GetYaxis()->SetTitle("Count");
   h_kshort_mass->SetTitle("In search of K_{S}^{0}");
   c1->Update();

   //Plot particle flow neutral hadron mass distribution
   
//   TH1D *h_pfNeutralHadron_mass = (TH1D*)ksFile->Get("kshort/h_pfNeutralHadron_mass");
//   TCanvas *c_pf_mass = new TCanvas("c_pf_mass","c_pf_mass",800,800);
//   h_pfNeutralHadron_mass->Draw();
//   c_pf_mass->Update();

   //Estimate background

   TH1D *ks_bckgnd = new TH1D("ks_bckgnd","K_{s}^{0} background estimate",h_kshort_mass->GetNbinsX(),h_kshort_mass->GetXaxis()->GetXmin(),h_kshort_mass->GetXaxis()->GetXmax());
   for (int i = 1; i <= h_kshort_mass->GetNbinsX(); i++) {
      double xVal = h_kshort_mass->GetBinCenter(i);
      double count = h_kshort_mass->GetBinContent(i);
         if (xVal < ks_5sigma_lower) {
            ks_bckgnd->Fill(xVal, count);
         } else if (ks_5sigma_upper < xVal) {
            ks_bckgnd->Fill(xVal, count);
         }
   }

   int numBins = ks_bckgnd->GetNbinsX();
       
   double sumBinContent = 0.0;
   double ks_bkgd_x_min = 0.0;
   double ks_bkgd_x_max = 0.0;
     
   for (int bin = 1; bin <= numBins; ++bin) {
       double binContent = ks_bckgnd->GetBinContent(bin);
       sumBinContent += binContent;
       double xVal = ks_bckgnd->GetBinCenter(bin);
       if (binContent > 0) {
          if (ks_bkgd_x_min == 0.0 || xVal < ks_bkgd_x_min) {
             ks_bkgd_x_min = xVal;
          }
          if (ks_bkgd_x_max == 0.0 || xVal > ks_bkgd_x_max) {
             ks_bkgd_x_max = xVal;
          }
       }
   }

   std::cout << "Background min: " << ks_bkgd_x_min << std::endl;
   std::cout << "Background max: " << ks_bkgd_x_max << std::endl;

   double rangeSum1 = ks_bckgnd->Integral(ks_bckgnd->GetXaxis()->FindBin(ks_bkgd_x_min), ks_bckgnd->GetXaxis()->FindBin(ks_5sigma_lower));
   double rangeSum2 = ks_bckgnd->Integral(ks_bckgnd->GetXaxis()->FindBin(ks_5sigma_upper),  ks_bckgnd->GetXaxis()->FindBin(ks_bkgd_x_max));
   double rangeSum = rangeSum1 + rangeSum2;
   double bkgd_avg = rangeSum / ((ks_bckgnd->GetXaxis()->FindBin(ks_5sigma_lower) - ks_bckgnd->GetXaxis()->FindBin(ks_bkgd_x_min)) + (ks_bckgnd->GetXaxis()->FindBin(ks_bkgd_x_max) - ks_bckgnd->GetXaxis()->FindBin(ks_5sigma_upper)));
   std::cout << "bkgd_avg : " << bkgd_avg << std::endl;
   double totalEst = bkgd_avg * (ks_bckgnd->GetXaxis()->FindBin(ks_bkgd_x_max) - ks_bckgnd->GetXaxis()->FindBin(ks_bkgd_x_min));

//   double shoulderFraction = rangeSum / totalEst;
   double shoulderFraction = ((ks_bkgd_x_max - ks_bkgd_x_min) - (ks_5sigma_upper - ks_5sigma_lower))/(ks_cut_u - ks_cut_l);

   std::cout << "Estimated shoulder fraction of total : " << shoulderFraction << std::endl;

   TCanvas *c_ks_mass_bckgnd = new TCanvas("c_ks_mass_bckgnd","c_ks_mass_bckgnd",800,800);
  

   TLine *l_bkgd_avg = new TLine(ks_bkgd_x_min, bkgd_avg, ks_bkgd_x_max, bkgd_avg);
   l_bkgd_avg->SetLineColor(kBlack);
   l_bkgd_avg->SetLineStyle(kDashed);

   ks_bckgnd->Draw("HIST");
   l_bkgd_avg->Draw("SAME");

   c_ks_mass_bckgnd->Update();

   //Get data from TFile Branches:
//   TTree *tree = (TTree*)ksFile->Get("kshort/tree");
//   Double_t my_kslifetime;
//   tree->SetBranchAddress("kslifetmes",&my_kslifetime);
//   Long64_t numEntries = tree->GetEntries();
//   for (Long64_t entry = 0; entry < numEntries; entry++ ) {
//      tree->GetEntry(entry);
//      std::cout << my_kslifetime << std::endl;
//   } 

//   ksFile->Close();

   TH1D *h_ks_bkgd_lifetimes = (TH1D*)ksFile->Get("kshort/h_ks_bkgd_lifetimes");
   h_ks_bkgd_lifetimes->SetTitle("Unscaled background");

   TCanvas *c_ks_bckgnd = new TCanvas("c_ks_bckgnd","c_ks_bckgnd",800,800);
   c_ks_bckgnd->SetLogy();
   h_ks_bkgd_lifetimes->Draw("HIST");
   c_ks_bckgnd->Update();

   h_ks_bkgd_lifetimes->Scale(1.0/shoulderFraction); //background scaled by the shoulder fraction to interpolate between shoulders

   //K-short distances
   TH1D *h_ksdistances = (TH1D*)ksFile->Get("kshort/h_ksdistances");
   TCanvas *cDist = new TCanvas("cDist", "K_{S}^{0} distance",800,800);
   cDist->SetLogy();
   h_ksdistances->GetXaxis()->SetTitle("Distance traveled before decay (lab) [cm]");
   h_ksdistances->GetYaxis()->SetTitle("Count");
   h_ksdistances->SetTitle("K_{S}^{0} distances");
   h_ksdistances->Draw("HIST");

   double cutoffIndex = h_ksdistances->GetMaximumBin();
   double cutoff = h_ksdistances->GetXaxis()->GetBinCenter(cutoffIndex);
   TLine *l_cutoff = new TLine(cutoff, 0, cutoff, h_ksdistances->GetMaximum());
   l_cutoff->SetLineColor(kRed);
   l_cutoff->SetLineStyle(kDotted);
   l_cutoff->Draw("SAME");

   TText* t_cutoff = new TText(cutoff, h_ksdistances->GetMaximum(), Form("%.2f", cutoff));
   t_cutoff->SetTextSize(0.03);
   t_cutoff->SetTextColor(kRed);
   t_cutoff->Draw("SAME");

   cDist->Update();


   //K-short lifetimes

   TH1D *h_kslifetimes = (TH1D*)ksFile->Get("kshort/h_kslifetimes");
   h_kslifetimes->Add(h_ks_bkgd_lifetimes, -1);//subtract background
   std::cout<<"Background: subtracted."<<std::endl;
   TCanvas *c2 = new TCanvas("c2","K_{S}^{0} lifetimes",800,800);
   c2->SetLogy();

   Double_t numEntries1 = h_kslifetimes->GetEntries();
   double initialParams1[2] = {pow(10,-11), numEntries1};
   TF1 *e_fit1 = new TF1("e_fit1",exponential,0.05*pow(10,-9),0.35*pow(10,-9),2); //higher range was 0.35*pow(10,-19)
   e_fit1->SetParameters(initialParams1);
   h_kslifetimes->Fit(e_fit1, "R");

   h_kslifetimes->GetXaxis()->SetTitle("Proper time ellapsed before decay [s]");
   h_kslifetimes->GetYaxis()->SetTitle("Count");
   h_kslifetimes->SetTitle("K_{S}^{0} lifetimes");

   h_kslifetimes->Draw("HIST");
   e_fit1->Draw("SAME");
   c2->Update();


   //Distance vs. Momentum scatterplot

   TH2D *h_ks_distVp = (TH2D*)ksFile->Get("kshort/h_ks_distVp");

   TCanvas *c_distVp = new TCanvas("c_distVp", "c_distVp", 800, 800);
   h_ks_distVp->SetTitle("Does #Deltax depend on |#vec{p}|?");
   h_ks_distVp->GetXaxis()->SetTitle("|#vec{p}| [GeV]");
   h_ks_distVp->GetYaxis()->SetTitle("#Deltax [cm]");
   h_ks_distVp->SetStats(kFALSE);
   h_ks_distVp->Draw("colorz");

   TLine *l_beampipe = new TLine(0, 2.2, 40, 2.2);
   l_beampipe->SetLineColor(kRed);
   l_beampipe->Draw("SAME");

   TText* t_beampipe = new TText(27, 2.2, Form("%.16s", "Beam pipe radius"));
   t_beampipe->SetTextSize(0.03);
   t_beampipe->SetTextColor(kRed);
   t_beampipe->Draw("SAME");

//   TLine *l_si = new TLine(0, 100, 40, 100);
//   l_si->SetLineColor(kGreen);
//   l_si->SetLineStyle(kDashed);
//   l_si->Draw("SAME");

//   TText* t_si = new TText(27, 100, Form("%.14s", "End Si tracker"));
//   t_si->SetTextSize(0.03);
//   t_si->SetTextColor(kGreen);
//   t_si->Draw("SAME");

   c_distVp->Update();


   //Lambda display

   TFile *lamFile = new TFile("kshortAnalyze.root");
   TH1D *h_lambda_mass = (TH1D*)lamFile->Get("kshort/h_lambda_mass");
   TCanvas *c3 = new TCanvas("c3","c3",800,800);

   h_lambda_mass->Draw();
   h_lambda_mass->GetXaxis()->SetTitle("Mass [GeV]");
   h_lambda_mass->GetYaxis()->SetTitle("Count");
   h_lambda_mass->SetTitle("In search of #Lambda baryons");
   c3->Update();

   TH1D *h_lamlifetimes = (TH1D*)lamFile->Get("kshort/h_lamlifetimes");
   TCanvas *c4 = new TCanvas("c4","#Lambda lifetimes",800,800);
   c4->SetLogy();

   Double_t numEntries2 = h_lamlifetimes->GetEntries();
   double initialParams2[2] = {pow(10,-10), numEntries2};
   TF1 *e_fit2 = new TF1("e_fit2",exponential,0.02*pow(10,-9),0.42*pow(10,-9),2);
   e_fit2->SetParameters(initialParams2);
   h_lamlifetimes->Fit(e_fit2, "R");

   h_lamlifetimes->GetXaxis()->SetTitle("Proper time ellapsed before decay [s]");
   h_lamlifetimes->GetYaxis()->SetTitle("Count");
   h_lamlifetimes->SetTitle("#Lambda lifetimes");

   h_lamlifetimes->Draw();
   c4->Update();

   //Dipion display

   TH1D *h_dipion_mass = (TH1D*)ksFile->Get("kshort/h_dipion_mass");
   TCanvas *dipiCan = new TCanvas("dipiCan", "dipiCan", 800, 800);
   h_dipion_mass->Draw();
   dipiCan->Update();

//   ksFile->Close();
}

int main() {
   display();
   return 0;
};
