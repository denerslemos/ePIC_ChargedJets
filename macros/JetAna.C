#include "make_plots.h"

#include <utility>  // for std::pair

std::pair<TGraphErrors*, TGraphErrors*> get_mean_stddev_customXranges(TH2D* h2d) {
    int ny = h2d->GetNbinsY();
    int nx = h2d->GetNbinsX();

    std::vector<double> yCenters, meanX, errMeanX, stddevX, errStddevX;

    for (int iy = 1; iy <= ny; ++iy) {

        TH1D* projX = h2d->ProjectionX(Form("projX_ybin%d", iy), iy, iy);

        double N = projX->GetEntries();
        if (N > 1) {
            yCenters.push_back(h2d->GetYaxis()->GetBinCenter(iy));
            meanX.push_back(projX->GetMean());
            errMeanX.push_back(projX->GetMeanError());
            stddevX.push_back(projX->GetStdDev());
            errStddevX.push_back(projX->GetStdDevError());
        }

        delete projX;
        h2d->GetXaxis()->SetRange(0, 0);
    }

    int n = yCenters.size();
    TGraphErrors* gMean = new TGraphErrors(n);
    TGraphErrors* gStd  = new TGraphErrors(n);

    for (int i = 0; i < n; ++i) {
        gMean->SetPoint(i, yCenters[i], meanX[i]);
        gMean->SetPointError(i, 0, errMeanX[i]);

        gStd->SetPoint(i, yCenters[i], stddevX[i]);
        gStd->SetPointError(i, 0, errStddevX[i]);
    }

    gMean->SetTitle("Mean of X per Y bin;Y;Mean of X");
    gStd->SetTitle("StdDev of X per Y bin;Y;StdDev of X");

    return std::make_pair(gMean, gStd);
}

TH1* SetContentToOneKeepErrors(TH1* hist, double adjustunc) {
    if (!hist) return nullptr;

    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        double oldError = hist->GetBinError(bin);
        hist->SetBinContent(bin, 1.0);
        hist->SetBinError(bin, oldError/adjustunc);
    }

    return hist;
}

void JetAna() {

	std::vector<string> histleg;
	histleg.push_back("");

	std::vector<string> histleg2;
	histleg2.push_back("eAu");
	histleg2.push_back("ep");

	double xseceAu = 0.3413; // ub
	double xsecep = 41.0; // ub

    // Input files
    TFile *ep_1to10 = TFile::Open("rootfiles/Jetep_q1to10.root");
    TFile *ep_10to100 = TFile::Open("rootfiles/Jetep_q10to100.root");
    TFile *eAu_1to10 = TFile::Open("rootfiles/JeteAu_q1to10.root");
    TFile *eAu_10to100 = TFile::Open("rootfiles/JeteAu_q10to100.root");
 
    TH1D* hEvep1to10 = (TH1D*)ep_1to10->Get("NEvents");
    TH1D* hEvep10to100 = (TH1D*)ep_10to100->Get("NEvents");
    TH1D* hEveAu1to10 = (TH1D*)eAu_1to10->Get("NEvents");
    TH1D* hEveAu10to100 = (TH1D*)eAu_10to100->Get("NEvents"); 
    
    double lumiep_1to10 = xsecep * double(hEvep1to10->GetEntries());
    double lumiep_10to100 = xsecep * double(hEvep10to100->GetEntries()); 
    double lumieAu_1to10 =  xseceAu * double(hEveAu1to10->GetEntries());
    double lumieAu_10to100 = xseceAu * double(hEveAu10to100->GetEntries());
    
    cout << lumieAu_10to100 << endl;
    
    // Let us do the JER/JES First
    THnSparseD *JESJERep1to10 = (THnSparseD*)ep_1to10->Get("mHistJESJERvsPt_noEMCP");
    THnSparseD *JESJERep10to100 = (THnSparseD*)ep_10to100->Get("mHistJESJERvsPt_noEMCP");
    THnSparseD *JESJEReAu1to10 = (THnSparseD*)eAu_1to10->Get("mHistJESJERvsPt_noEMCP");
    THnSparseD *JESJEReAu10to100 = (THnSparseD*)eAu_10to100->Get("mHistJESJERvsPt_noEMCP"); 

	JESJERep1to10->GetAxis(0)->SetRangeUser(0.9, 1.1);
	JESJERep10to100->GetAxis(0)->SetRangeUser(0.9, 1.1);
	JESJEReAu1to10->GetAxis(0)->SetRangeUser(0.9, 1.1);
	JESJEReAu10to100->GetAxis(0)->SetRangeUser(0.9, 1.1);

	JESJERep1to10->GetAxis(1)->SetRangeUser(2.0, 50.0);
	JESJERep10to100->GetAxis(1)->SetRangeUser(2.0, 50.0);
	JESJEReAu1to10->GetAxis(1)->SetRangeUser(2.0, 50.0);
	JESJEReAu10to100->GetAxis(1)->SetRangeUser(2.0, 50.0);

	JESJERep1to10->GetAxis(2)->SetRangeUser(-0.99, 0.99);
	JESJERep10to100->GetAxis(2)->SetRangeUser(-0.99, 0.99);
	JESJEReAu1to10->GetAxis(2)->SetRangeUser(-0.99, 0.99);
	JESJEReAu10to100->GetAxis(2)->SetRangeUser(-0.99, 0.99);

	TH2D* JESJERep1to10B_2D = (TH2D*) JESJERep1to10->Projection(1,0);
	TH2D* JESJERep10to100B_2D = (TH2D*) JESJERep10to100->Projection(1,0);
	TH2D* JESJEReAu1to10B_2D = (TH2D*) JESJEReAu1to10->Projection(1,0);
	TH2D* JESJEReAu10to100B_2D = (TH2D*) JESJEReAu10to100->Projection(1,0);
	
	JESJERep1to10B_2D->Draw();
	
	JESJERep1to10->GetAxis(2)->SetRangeUser(1.01, 2.99);
	JESJERep10to100->GetAxis(2)->SetRangeUser(1.01, 2.99);
	JESJEReAu1to10->GetAxis(2)->SetRangeUser(1.01, 2.99);
	JESJEReAu10to100->GetAxis(2)->SetRangeUser(1.01, 2.99);

	TH2D* JESJERep1to10E_2D = (TH2D*) JESJERep1to10->Projection(1,0);
	TH2D* JESJERep10to100E_2D = (TH2D*) JESJERep10to100->Projection(1,0);
	TH2D* JESJEReAu1to10E_2D = (TH2D*) JESJEReAu1to10->Projection(1,0);
	TH2D* JESJEReAu10to100E_2D = (TH2D*) JESJEReAu10to100->Projection(1,0);


	JESJERep1to10->GetAxis(2)->SetRangeUser(-2.99, -1.01);
	JESJERep10to100->GetAxis(2)->SetRangeUser(-2.99, -1.01);
	JESJEReAu1to10->GetAxis(2)->SetRangeUser(-2.99, -1.01);
	JESJEReAu10to100->GetAxis(2)->SetRangeUser(-2.99, -1.01);

	TH2D* JESJERep1to10E1_2D = (TH2D*) JESJERep1to10->Projection(1,0);
	TH2D* JESJERep10to100E1_2D = (TH2D*) JESJERep10to100->Projection(1,0);
	TH2D* JESJEReAu1to10E1_2D = (TH2D*) JESJEReAu1to10->Projection(1,0);
	TH2D* JESJEReAu10to100E1_2D = (TH2D*) JESJEReAu10to100->Projection(1,0);
	
	JESJERep1to10E_2D->Add(JESJERep1to10E1_2D);
	JESJERep10to100->Add(JESJERep10to100E1_2D);
	JESJEReAu1to10E_2D->Add(JESJEReAu1to10E1_2D);
	JESJEReAu10to100E_2D->Add(JESJEReAu10to100E1_2D);
	
	
	// Get mean/stddev graphs from your function
	auto [mean1, std1] = get_mean_stddev_customXranges(JESJERep1to10B_2D);     // ep 1to10
	auto [mean2, std2] = get_mean_stddev_customXranges(JESJERep10to100B_2D);      // ep 10to100
	auto [mean3, std3] = get_mean_stddev_customXranges(JESJEReAu1to10B_2D);   // eAu 1to10
	auto [mean4, std4] = get_mean_stddev_customXranges(JESJEReAu10to100B_2D); // eAu 10to100

// Set styles for mean graphs
mean1->SetMarkerStyle(25); mean1->SetMarkerColor(kBlue+3);    mean1->SetLineColor(kBlue+3);
mean2->SetMarkerStyle(24); mean2->SetMarkerColor(kRed+3);     mean2->SetLineColor(kRed+3);
mean3->SetMarkerStyle(21); mean3->SetMarkerColor(kBlue+3);    mean3->SetLineColor(kBlue+3);
mean4->SetMarkerStyle(20); mean4->SetMarkerColor(kRed+3);     mean4->SetLineColor(kRed+3);

// Set styles for stddev graphs
std1->SetMarkerStyle(25); std1->SetMarkerColor(kBlue+3);    std1->SetLineColor(kBlue+3);
std2->SetMarkerStyle(24); std2->SetMarkerColor(kRed+3);     std2->SetLineColor(kRed+3);
std3->SetMarkerStyle(21); std3->SetMarkerColor(kBlue+3);    std3->SetLineColor(kBlue+3);
std4->SetMarkerStyle(20); std4->SetMarkerColor(kRed+3);     std4->SetLineColor(kRed+3);


// Mean canvas and graph axis range
TCanvas* cMean = new TCanvas("cMean", "Mean vs Y", 1200, 600);
mean1->SetTitle("");
mean1->GetYaxis()->SetRangeUser(0.9, 1.1);   // Set Y range between 0.8 and 1.2
mean1->GetXaxis()->SetRangeUser(0.0, 18.0);  // Set Y range between 0 and 0.1
mean1->GetYaxis()->SetTitle("JES");  // Set Y range between 0 and 0.1
mean1->GetXaxis()->SetTitle("p^{Gen}_{T} [GeV]");  // Set Y range between 0 and 0.1
mean1->Draw("AP");
mean2->Draw("P SAME");
mean3->Draw("P SAME");
mean4->Draw("P SAME");
// For mean canvas
TLatex latexMean;
latexMean.SetNDC();  // Use normalized coordinates (0-1)
latexMean.SetTextSize(0.04);
latexMean.DrawLatex(0.15, 0.85, "|#eta| < 1");

// Draw horizontal line at y=1
TLine* line1 = new TLine(mean1->GetXaxis()->GetXmin(), 1.0, 18.0, 1.0);
line1->SetLineColor(kBlack);
line1->SetLineStyle(2);  // dashed line
line1->Draw("SAME");

// Legend for mean
auto legMean = new TLegend(0.7,0.7,0.9,0.9);
legMean->AddEntry(mean1, "ep, q^{2} > 1", "p");
legMean->AddEntry(mean2, "ep, q^{2} > 10", "p");
legMean->AddEntry(mean3, "eAu, q^{2} = [1-10]", "p");
legMean->AddEntry(mean4, "eAu, q^{2} = [10-100]", "p");
legMean->Draw();

// Stddev canvas and graph axis range
TCanvas* cStd = new TCanvas("cStd", "StdDev vs Y", 1200, 600);
std1->SetTitle("");
std1->GetYaxis()->SetRangeUser(0.0, 0.07);  // Set Y range between 0 and 0.1
std1->GetXaxis()->SetRangeUser(0.0, 18.0);  // Set Y range between 0 and 0.1
std1->GetYaxis()->SetTitle("JER");  // Set Y range between 0 and 0.1
std1->GetXaxis()->SetTitle("p^{Gen}_{T} [GeV]");  // Set Y range between 0 and 0.1
std1->Draw("AP");
std2->Draw("P SAME");
std3->Draw("P SAME");
std4->Draw("P SAME");
// For stddev canvas
TLatex latexStd;
latexStd.SetNDC();
latexStd.SetTextSize(0.04);
latexStd.DrawLatex(0.15, 0.85, "|#eta| < 1");

// Legend for stddev
auto legStd = new TLegend(0.7,0.7,0.9,0.9);
legStd->AddEntry(std1, "ep, q^{2} > 1", "p");
legStd->AddEntry(std2, "ep, q^{2} > 10", "p");
legStd->AddEntry(std3, "eAu, q^{2} = [1-10]", "p");
legStd->AddEntry(std4, "eAu, q^{2} = [10-100]", "p");
legStd->Draw();

cMean->Update();
cStd->Update();


// Get mean/stddev graphs from your function
auto [mean_epLow, std_epLow] = get_mean_stddev_customXranges(JESJERep1to10E_2D);     // ep 1to10
auto [mean_epHigh, std_epHigh] = get_mean_stddev_customXranges(JESJERep10to100E_2D);  // ep 10to100
auto [mean_eAuLow, std_eAuLow] = get_mean_stddev_customXranges(JESJEReAu1to10E_2D);   // eAu 1to10
auto [mean_eAuHigh, std_eAuHigh] = get_mean_stddev_customXranges(JESJEReAu10to100E_2D); // eAu 10to100

// Set styles for mean graphs
mean_epLow->SetMarkerStyle(25); mean_epLow->SetMarkerColor(kBlue+3);    mean_epLow->SetLineColor(kBlue+3);
mean_epHigh->SetMarkerStyle(24); mean_epHigh->SetMarkerColor(kRed+3);     mean_epHigh->SetLineColor(kRed+3);
mean_eAuLow->SetMarkerStyle(21); mean_eAuLow->SetMarkerColor(kBlue+3);    mean_eAuLow->SetLineColor(kBlue+3);
mean_eAuHigh->SetMarkerStyle(20); mean_eAuHigh->SetMarkerColor(kRed+3);     mean_eAuHigh->SetLineColor(kRed+3);

// Set styles for stddev graphs
std_epLow->SetMarkerStyle(25); std_epLow->SetMarkerColor(kBlue+3);    std_epLow->SetLineColor(kBlue+3);
std_epHigh->SetMarkerStyle(24); std_epHigh->SetMarkerColor(kRed+3);     std_epHigh->SetLineColor(kRed+3);
std_eAuLow->SetMarkerStyle(21); std_eAuLow->SetMarkerColor(kBlue+3);    std_eAuLow->SetLineColor(kBlue+3);
std_eAuHigh->SetMarkerStyle(20); std_eAuHigh->SetMarkerColor(kRed+3);     std_eAuHigh->SetLineColor(kRed+3);

// Mean canvas and graph axis range
TCanvas* canvasMean = new TCanvas("canvasMean", "Mean vs Y", 1200, 600);
mean_epLow->SetTitle("");
mean_epLow->GetYaxis()->SetRangeUser(0.9, 1.1);
mean_epLow->GetXaxis()->SetRangeUser(0.0, 18.0);
mean_epLow->GetYaxis()->SetTitle("JES");
mean_epLow->GetXaxis()->SetTitle("p^{Gen}_{T} [GeV]");
mean_epLow->Draw("AP");
mean_epHigh->Draw("P SAME");
mean_eAuLow->Draw("P SAME");
mean_eAuHigh->Draw("P SAME");

// Add text "1 < |η| < 3"
TLatex latexTextMean;
latexTextMean.SetNDC();
latexTextMean.SetTextSize(0.04);
latexTextMean.DrawLatex(0.15, 0.85, "1 < |#eta| < 3");

// Draw horizontal line at y=1
TLine* refLine = new TLine(mean_epLow->GetXaxis()->GetXmin(), 1.0, 18.0, 1.0);
refLine->SetLineColor(kBlack);
refLine->SetLineStyle(2);
refLine->Draw("SAME");

// Legend for mean
auto legendMean = new TLegend(0.7, 0.7, 0.9, 0.9);
legendMean->AddEntry(mean_epLow, "ep, q^{2} > 1", "p");
legendMean->AddEntry(mean_epHigh, "ep, q^{2} > 10", "p");
legendMean->AddEntry(mean_eAuLow, "eAu, q^{2} = [1-10]", "p");
legendMean->AddEntry(mean_eAuHigh, "eAu, q^{2} = [10-100]", "p");
legendMean->Draw();

// Stddev canvas and graph axis range
TCanvas* canvasStd = new TCanvas("canvasStd", "StdDev vs Y", 1200, 600);
std_epLow->SetTitle("");
std_epLow->GetYaxis()->SetRangeUser(0.0, 0.07);
std_epLow->GetXaxis()->SetRangeUser(0.0, 18.0);
std_epLow->GetYaxis()->SetTitle("JER");
std_epLow->GetXaxis()->SetTitle("p^{Gen}_{T} [GeV]");
std_epLow->Draw("AP");
std_epHigh->Draw("P SAME");
std_eAuLow->Draw("P SAME");
std_eAuHigh->Draw("P SAME");

// Add text "1 < |η| < 3"
TLatex latexTextStd;
latexTextStd.SetNDC();
latexTextStd.SetTextSize(0.04);
latexTextStd.DrawLatex(0.15, 0.85, "1 < |#eta| < 3");

// Legend for stddev
auto legendStd = new TLegend(0.7, 0.7, 0.9, 0.9);
legendStd->AddEntry(std_epLow, "ep, q^{2} > 1", "p");
legendStd->AddEntry(std_epHigh, "ep, q^{2} > 10", "p");
legendStd->AddEntry(std_eAuLow, "eAu, q^{2} = [1-10]", "p");
legendStd->AddEntry(std_eAuHigh, "eAu, q^{2} = [10-100]", "p");
legendStd->Draw();

canvasMean->Update();
canvasStd->Update();

// Jets
// -> pT for |eta| < 3
// -> eta and phi for pT > 5 GeV

// Let us do the JER/JES First
THnSparseD *JetRecoep1to10 = (THnSparseD*)ep_1to10->Get("mHistJetReco_noEMCP");
THnSparseD *JetRecoep10to100 = (THnSparseD*)ep_10to100->Get("mHistJetReco_noEMCP");
THnSparseD *JetRecoeAu1to10 = (THnSparseD*)eAu_1to10->Get("mHistJetReco_noEMCP");
THnSparseD *JetRecoeAu10to100 = (THnSparseD*)eAu_10to100->Get("mHistJetReco_noEMCP"); 

THnSparseD *JetGenep1to10 = (THnSparseD*)ep_1to10->Get("mHistJetGen_noE");
THnSparseD *JetGenep10to100 = (THnSparseD*)ep_10to100->Get("mHistJetGen_noE");
THnSparseD *JetGeneAu1to10 = (THnSparseD*)eAu_1to10->Get("mHistJetGen_noE");
THnSparseD *JetGeneAu10to100 = (THnSparseD*)eAu_10to100->Get("mHistJetGen_noE"); 

THnSparseD *TrkRecoep1to10 = (THnSparseD*)ep_1to10->Get("mHistFrgConstReco_noEMCP");
THnSparseD *TrkRecoep10to100 = (THnSparseD*)ep_10to100->Get("mHistFrgConstReco_noEMCP");
THnSparseD *TrkRecoeAu1to10 = (THnSparseD*)eAu_1to10->Get("mHistFrgConstReco_noEMCP");
THnSparseD *TrkRecoeAu10to100 = (THnSparseD*)eAu_10to100->Get("mHistFrgConstReco_noEMCP"); 

THnSparseD *TrkGenep1to10 = (THnSparseD*)ep_1to10->Get("mHistFrgConstGen_noE");
THnSparseD *TrkGenep10to100 = (THnSparseD*)ep_10to100->Get("mHistFrgConstGen_noE");
THnSparseD *TrkGeneAu1to10 = (THnSparseD*)eAu_1to10->Get("mHistFrgConstGen_noE");
THnSparseD *TrkGeneAu10to100 = (THnSparseD*)eAu_10to100->Get("mHistFrgConstGen_noE"); 


THnSparseD *Trk4VRecoep1to10 = (THnSparseD*)ep_1to10->Get("mHistTrkConstReco_noEMCP");
THnSparseD *Trk4VRecoep10to100 = (THnSparseD*)ep_10to100->Get("mHistTrkConstReco_noEMCP");
THnSparseD *Trk4VRecoeAu1to10 = (THnSparseD*)eAu_1to10->Get("mHistTrkConstReco_noEMCP");
THnSparseD *Trk4VRecoeAu10to100 = (THnSparseD*)eAu_10to100->Get("mHistTrkConstReco_noEMCP"); 

THnSparseD *Trk4VGenep1to10 = (THnSparseD*)ep_1to10->Get("mHistTrkConstGen_noE");
THnSparseD *Trk4VGenep10to100 = (THnSparseD*)ep_10to100->Get("mHistTrkConstGen_noE");
THnSparseD *Trk4VGeneAu1to10 = (THnSparseD*)eAu_1to10->Get("mHistTrkConstGen_noE");
THnSparseD *Trk4VGeneAu10to100 = (THnSparseD*)eAu_10to100->Get("mHistTrkConstGen_noE"); 

//for ReAu extrapolation
JetRecoep10to100->GetAxis(1)->SetRangeUser(-3,3);
TH1D* JetpTDenC = (TH1D*)JetRecoep10to100->Projection(0);
JetRecoeAu10to100->GetAxis(1)->SetRangeUser(-3,3);
TH1D* JetpTNumC = (TH1D*)JetRecoeAu10to100->Projection(0);
JetpTNumC->Divide(JetpTDenC);
SetContentToOneKeepErrors(JetpTNumC,1.0);

std::vector<TH1 *> hstatJetpT;
hstatJetpT.push_back(JetpTNumC);
histo_plots(hstatJetpT, "E", 0.0001, 20.0, 0.9, 1.1, false, false, "Jet p_{T} [GeV]", "R_{eAu} Stat. Uncer.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg, 0.0, 0.2, 0.0, 0.2, "JetReAupT" );

JetRecoep10to100->GetAxis(1)->SetRange(0,0);
JetRecoeAu10to100->GetAxis(1)->SetRange(0,0);
JetRecoep10to100->GetAxis(0)->SetRangeUser(5.0,50.0);
JetRecoeAu10to100->GetAxis(0)->SetRangeUser(5.0,50.0);
TH1D* JetetaDenC = (TH1D*)JetRecoep10to100->Projection(1);
TH1D* JetetaNumC = (TH1D*)JetRecoeAu10to100->Projection(1);
JetetaNumC->Divide(JetetaDenC);
SetContentToOneKeepErrors(JetetaNumC,1.0);
std::vector<TH1 *> hstatJeteta;
hstatJeteta.push_back(JetetaNumC);
histo_plots(hstatJeteta, "E", -2.99, 2.99, 0.9, 1.1, false, false, "Jet #eta ", "R_{eAu} Stat. Uncer.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg, 0.0, 0.2, 0.0, 0.2, "JetReAuEta" );

// const
TH1D* ConsDenC = (TH1D*)TrkRecoep10to100->Projection(0); ConsDenC->Rebin();
TH1D* ConsNumC = (TH1D*)TrkRecoeAu10to100->Projection(0); ConsNumC->Rebin();
ConsNumC->Divide(ConsDenC);
SetContentToOneKeepErrors(ConsNumC,1.0);
std::vector<TH1 *> hstatTrkcons;
hstatTrkcons.push_back(ConsNumC);
histo_plots(hstatTrkcons, "E", 0.0001, 0.9999, 0.9, 1.1, false, false, "z = p^{hadron}_{T} / p^{jet}_{T}", "R_{eAu} Stat. Uncer.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg, 0.0, 0.2, 0.0, 0.2, "JetReAuCons" );


JetRecoep10to100->GetAxis(1)->SetRange(0,0);
JetRecoeAu10to100->GetAxis(1)->SetRange(0,0);
JetRecoep10to100->GetAxis(0)->SetRange(0,0);
JetRecoeAu10to100->GetAxis(0)->SetRange(0,0);

JetRecoep10to100->GetAxis(1)->SetRangeUser(-2.99,2.99);
JetRecoeAu10to100->GetAxis(1)->SetRangeUser(-2.99,2.99);

TH1D* JetRecoep10to100pT1D = (TH1D*)JetRecoep10to100->Projection(0); JetRecoep10to100pT1D->Scale(1./JetRecoep10to100pT1D->Integral());
TH1D* JetRecoeAu10to100pT1D = (TH1D*)JetRecoeAu10to100->Projection(0); JetRecoeAu10to100pT1D->Scale(1./JetRecoeAu10to100pT1D->Integral());

std::vector<TH1 *> hJetpT;
hJetpT.push_back(JetRecoeAu10to100pT1D);
hJetpT.push_back(JetRecoep10to100pT1D);
histo_plots(hJetpT, "E", 0.0001, 20.0, 0.0, 1.0, false, true, "Jet p_{T} [GeV]", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.2, 0.2, 0.4, 0.4 , "JetPtReco" );


JetGenep10to100->GetAxis(1)->SetRangeUser(-2.99,2.99);
JetGeneAu10to100->GetAxis(1)->SetRangeUser(-2.99,2.99);

TH1D* JetGenep10to100pT1D = (TH1D*)JetGenep10to100->Projection(0); JetGenep10to100pT1D->Scale(1./JetGenep10to100pT1D->Integral());
TH1D* JetGeneAu10to100pT1D = (TH1D*)JetGeneAu10to100->Projection(0); JetGeneAu10to100pT1D->Scale(1./JetGeneAu10to100pT1D->Integral());

std::vector<TH1 *> hJetpTG;
hJetpTG.push_back(JetGeneAu10to100pT1D);
hJetpTG.push_back(JetGenep10to100pT1D);
histo_plots(hJetpTG, "E", 0.0001, 20.0, 0.0, 1.0, false, true, "Jet p_{T} [GeV]", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.2, 0.2, 0.4, 0.4 , "JetPtGen" );

JetRecoep10to100->GetAxis(1)->SetRange(0,0);
JetRecoeAu10to100->GetAxis(1)->SetRange(0,0);
JetGenep10to100->GetAxis(1)->SetRange(0,0);
JetGeneAu10to100->GetAxis(1)->SetRange(0,0);

JetRecoep10to100->GetAxis(0)->SetRangeUser(5.0,50.0);
JetRecoeAu10to100->GetAxis(0)->SetRangeUser(5.0,50.0);
JetRecoep10to100->GetAxis(0)->SetRangeUser(5.0,50.0);
JetRecoeAu10to100->GetAxis(0)->SetRangeUser(5.0,50.0);

TH1D* JetRecoep10to100Eta1D = (TH1D*)JetRecoep10to100->Projection(1); JetRecoep10to100Eta1D->Scale(1./JetRecoep10to100Eta1D->Integral());
TH1D* JetRecoeAu10to100Eta1D = (TH1D*)JetRecoeAu10to100->Projection(1); JetRecoeAu10to100Eta1D->Scale(1./JetRecoeAu10to100Eta1D->Integral());
std::vector<TH1 *> hJetEta;
hJetEta.push_back(JetRecoeAu10to100Eta1D);
hJetEta.push_back(JetRecoep10to100Eta1D);
histo_plots(hJetEta, "E", -4.0, 4.0, 0.0, 1.0, false, true, "Jet #eta", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.6, 0.2, 0.8, 0.4 , "JetEtaReco" );

TH1D* JetGenep10to100Eta1D = (TH1D*)JetGenep10to100->Projection(1); JetGenep10to100Eta1D->Scale(1./JetGenep10to100Eta1D->Integral());
TH1D* JetGeneAu10to100Eta1D = (TH1D*)JetGeneAu10to100->Projection(1); JetGeneAu10to100Eta1D->Scale(1./JetGeneAu10to100Eta1D->Integral());
std::vector<TH1 *> hJetEtaG;
hJetEtaG.push_back(JetGeneAu10to100Eta1D);
hJetEtaG.push_back(JetGenep10to100Eta1D);
histo_plots(hJetEtaG, "E", -4.0, 4.0, 0.0, 1.0, false, true, "Jet #eta", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.6, 0.2, 0.8, 0.4 , "JetEtaGen" );

JetRecoep10to100->GetAxis(1)->SetRangeUser(-2.99,2.99);
JetRecoeAu10to100->GetAxis(1)->SetRangeUser(-2.99,2.99);
JetGenep10to100->GetAxis(1)->SetRangeUser(-2.99,2.99);
JetGeneAu10to100->GetAxis(1)->SetRangeUser(-2.99,2.99);

TH1D* JetRecoep10to100Phi1D = (TH1D*)JetRecoep10to100->Projection(2); JetRecoep10to100Phi1D->Scale(1./JetRecoep10to100Phi1D->Integral());
TH1D* JetRecoeAu10to100Phi1D = (TH1D*)JetRecoeAu10to100->Projection(2); JetRecoeAu10to100Phi1D->Scale(1./JetRecoeAu10to100Phi1D->Integral());
std::vector<TH1 *> hJetPhi;
hJetPhi.push_back(JetRecoeAu10to100Phi1D);
hJetPhi.push_back(JetRecoep10to100Phi1D);
histo_plots(hJetPhi, "E", -3.2, 3.2, 0.0, 1.0, false, true, "Jet #phi", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.4, 0.6, 0.6, 0.8 , "JetPhiReco" );

TH1D* JetGenep10to100Phi1D = (TH1D*)JetGenep10to100->Projection(2); JetGenep10to100Phi1D->Scale(1./JetGenep10to100Phi1D->Integral());
TH1D* JetGeneAu10to100Phi1D = (TH1D*)JetGeneAu10to100->Projection(2); JetGeneAu10to100Phi1D->Scale(1./JetGeneAu10to100Phi1D->Integral());
std::vector<TH1 *> hJetPhiG;
hJetPhiG.push_back(JetGeneAu10to100Phi1D);
hJetPhiG.push_back(JetGenep10to100Phi1D);
histo_plots(hJetPhiG, "E", -3.2, 3.2, 0.0, 1.0, false, true, "Jet #phi", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.4, 0.6, 0.6, 0.8 , "JetPhiGen" );

TH1D* JetRecoep10to100Mass1D = (TH1D*)JetRecoep10to100->Projection(3); JetRecoep10to100Mass1D->Scale(1./JetRecoep10to100Mass1D->Integral());
TH1D* JetRecoeAu10to100Mass1D = (TH1D*)JetRecoeAu10to100->Projection(3); JetRecoeAu10to100Mass1D->Scale(1./JetRecoeAu10to100Mass1D->Integral());
std::vector<TH1 *> hJetMass;
hJetMass.push_back(JetRecoeAu10to100Mass1D);
hJetMass.push_back(JetRecoep10to100Mass1D);
histo_plots(hJetMass, "E", 0.001, 7.0, 0.0, 1.0, false, true, "Jet Mass [GeV]", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.2, 0.2, 0.4, 0.4  , "JetMassReco" );

TH1D* JetGenep10to100Mass1D = (TH1D*)JetGenep10to100->Projection(3); JetGenep10to100Mass1D->Scale(1./JetGenep10to100Mass1D->Integral());
TH1D* JetGeneAu10to100Mass1D = (TH1D*)JetGeneAu10to100->Projection(3); JetGeneAu10to100Mass1D->Scale(1./JetGeneAu10to100Mass1D->Integral());
std::vector<TH1 *> hJetMassG;
hJetMassG.push_back(JetGeneAu10to100Mass1D);
hJetMassG.push_back(JetGenep10to100Mass1D);
histo_plots(hJetMassG, "E", 0.001, 7.0, 0.0, 1.0, false, true, "Jet Mass [GeV]", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.2, 0.2, 0.4, 0.4 , "JetMassGen" );


TH1D* Trk4VRecoep10to100pT1D = (TH1D*)Trk4VRecoep10to100->Projection(0); Trk4VRecoep10to100pT1D->Scale(1./Trk4VRecoep10to100pT1D->Integral());
TH1D* Trk4VRecoeAu10to100pT1D = (TH1D*)Trk4VRecoeAu10to100->Projection(0); Trk4VRecoeAu10to100pT1D->Scale(1./Trk4VRecoeAu10to100pT1D->Integral());
std::vector<TH1 *> hTrk4VpT;
hTrk4VpT.push_back(Trk4VRecoeAu10to100pT1D);
hTrk4VpT.push_back(Trk4VRecoep10to100pT1D);
histo_plots(hTrk4VpT, "E", 0.0001, 20.0, 0.0, 1.0, false, true, "Trk p_{T} [GeV]", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.2, 0.2, 0.4, 0.4 , "Trk4VPtReco" );

TH1D* Trk4VGenep10to100pT1D = (TH1D*)Trk4VGenep10to100->Projection(0); Trk4VGenep10to100pT1D->Scale(1./Trk4VGenep10to100pT1D->Integral());
TH1D* Trk4VGeneAu10to100pT1D = (TH1D*)Trk4VGeneAu10to100->Projection(0); Trk4VGeneAu10to100pT1D->Scale(1./Trk4VGeneAu10to100pT1D->Integral());
std::vector<TH1 *> hTrk4VpTG;
hTrk4VpTG.push_back(Trk4VGeneAu10to100pT1D);
hTrk4VpTG.push_back(Trk4VGenep10to100pT1D);
histo_plots(hTrk4VpTG, "E", 0.0001, 20.0, 0.0, 1.0, false, true, "Trk p_{T} [GeV]", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.2, 0.2, 0.4, 0.4 , "Trk4VPtGen" );


TH1D* Trk4VRecoep10to100Eta1D = (TH1D*)Trk4VRecoep10to100->Projection(1); Trk4VRecoep10to100Eta1D->Scale(1./Trk4VRecoep10to100Eta1D->Integral());
TH1D* Trk4VRecoeAu10to100Eta1D = (TH1D*)Trk4VRecoeAu10to100->Projection(1); Trk4VRecoeAu10to100Eta1D->Scale(1./Trk4VRecoeAu10to100Eta1D->Integral());
std::vector<TH1 *> hTrk4VEta;
hTrk4VEta.push_back(Trk4VRecoeAu10to100Eta1D);
hTrk4VEta.push_back(Trk4VRecoep10to100Eta1D);
histo_plots(hTrk4VEta, "E", -4.0, 4.0, 0.0, 1.0, false, true, "Trk #eta", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.6, 0.2, 0.8, 0.4 , "Trk4VEtaReco" );

TH1D* Trk4VGenep10to100Eta1D = (TH1D*)Trk4VGenep10to100->Projection(1); Trk4VGenep10to100Eta1D->Scale(1./Trk4VGenep10to100Eta1D->Integral());
TH1D* Trk4VGeneAu10to100Eta1D = (TH1D*)Trk4VGeneAu10to100->Projection(1); Trk4VGeneAu10to100Eta1D->Scale(1./Trk4VGeneAu10to100Eta1D->Integral());
std::vector<TH1 *> hTrk4VEtaG;
hTrk4VEtaG.push_back(Trk4VGeneAu10to100Eta1D);
hTrk4VEtaG.push_back(Trk4VGenep10to100Eta1D);
histo_plots(hTrk4VEtaG, "E", -4.0, 4.0, 0.0, 1.0, false, true, "Trk #eta", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.6, 0.2, 0.8, 0.4 , "Trk4VEtaGen" );


TH1D* Trk4VRecoep10to100Phi1D = (TH1D*)Trk4VRecoep10to100->Projection(2); Trk4VRecoep10to100Phi1D->Scale(1./Trk4VRecoep10to100Phi1D->Integral());
TH1D* Trk4VRecoeAu10to100Phi1D = (TH1D*)Trk4VRecoeAu10to100->Projection(2); Trk4VRecoeAu10to100Phi1D->Scale(1./Trk4VRecoeAu10to100Phi1D->Integral());
std::vector<TH1 *> hTrk4VPhi;
hTrk4VPhi.push_back(Trk4VRecoeAu10to100Phi1D);
hTrk4VPhi.push_back(Trk4VRecoep10to100Phi1D);
histo_plots(hTrk4VPhi, "E", -3.2, 3.2, 0.0, 1.0, false, true, "Trk #phi", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.4, 0.6, 0.6, 0.8 , "Trk4VPhiReco" );

TH1D* Trk4VGenep10to100Phi1D = (TH1D*)Trk4VGenep10to100->Projection(2); Trk4VGenep10to100Phi1D->Scale(1./Trk4VGenep10to100Phi1D->Integral());
TH1D* Trk4VGeneAu10to100Phi1D = (TH1D*)Trk4VGeneAu10to100->Projection(2); Trk4VGeneAu10to100Phi1D->Scale(1./Trk4VGeneAu10to100Phi1D->Integral());
std::vector<TH1 *> hTrk4VPhiG;
hTrk4VPhiG.push_back(Trk4VGeneAu10to100Phi1D);
hTrk4VPhiG.push_back(Trk4VGenep10to100Phi1D);
histo_plots(hTrk4VPhiG, "E", -3.2, 3.2, 0.0, 1.0, false, true, "Trk #phi", "A.U.", "", false, " ", "", 0.0, 0.0, "", 0.35, 0.22, "", 0.0, 0.0, histleg2, 0.4, 0.6, 0.6, 0.8 , "Trk4VPhiGen" );
 
 
}
