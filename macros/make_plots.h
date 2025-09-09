// general calling
// ROOT stuff
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom3.h>
#include "TFrame.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "THnBase.h"
#include <TLorentzVector.h>
#include "THnSparse.h"
#include "TRandom3.h"
#include <TCanvas.h>
#include <TH1D.h>
// C++ stuff
#include <vector>
#include <cstring>
#include <ctime>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <vector>
#include <map>

// defined colors
enum DSLcolors
{
	black,
	red,
	green,
	blue,
	darkred,
	darkgreen,
	darkblue,
	yellow,
	orange,
	grey,
	purple,
	pink
};

int DSL_colors[12] = {
	TColor::GetColor(0, 0, 0),
	TColor::GetColor(230, 0, 0),
	TColor::GetColor(0, 190, 0),
	TColor::GetColor(0, 200, 255),
	TColor::GetColor(150, 0, 0),
	TColor::GetColor(0, 100, 0),
	TColor::GetColor(0, 0, 200),
	TColor::GetColor(230, 180, 30),
	TColor::GetColor(250, 100, 0),
	TColor::GetColor(80, 80, 85),
	TColor::GetColor(85, 40, 80),
	TColor::GetColor(250, 20, 250)};

void histo_plots(std::vector<TH1 *> histo, string histstyle, double h_xmin, double h_xmax, double h_ymin, double h_ymax, bool usexlogscale, bool useylogscale, string xlabel, string ylabel, string cmstext, bool insideofcanvas, string cmssystem, string randomtext, double rnd_x, double rnd_y, string ptrange, double pt_x, double pt_y, string multrange, double mult_x, double mult_y, std::vector<string> histo_legend, double leg_xmin, double leg_ymin, double leg_xmax, double leg_ymax, string outfilename)
{

	/*
	histo -> vector of histograms
	histstyle -> for plots --> Use "E" for symbols (circles, squares, ...) and "HIST" for histograms
	h_xmin, h_xmax, h_ymin, h_ymax -> ranges of x and y labels
	usexlogscale and useylogscale -> logscale in x and/or y --> true means logscale on
	cmstext -> this is the "CMS Preliminary" it depends on the status of the analysis
		--> we usualy start with Preliminary -> if needed to be shown on thesis or internal conference use Work In Progress
		--> at the end we do not need to use this (final results/plots), so just set ""
	insideofcanvas --> text inside or outside the plot. You decide what is better.
	cmssystem --> collision, energy and luminosity. I prefer like: pPb 186 nb^{-1} (8.16 TeV)
	randomtext, rnd_x and rnd_y -> some text inside of the canvas, for example MC model or so. If you do not want just keep "". x and y represents the positions in the plot.
	ptrange, pt_x and pt_y -> some text for pT range (if needed) inside of the canvas. x and y represents the positions in the plot.
	multrange mult_x and  mult_y -> similar to pT but for multiplicity
	histo_legend -> vector with name of histograms to include in the legend
	leg_xmin, leg_ymin, leg_xmax, leg_ymax -> legend position
	outfilename -> output file name
	*/

	TCanvas *cc = new TCanvas("cc", "cc", 900, 900);
	cc->cd(1);
	cc->cd(1)->SetTopMargin(0.1);
	cc->cd(1)->SetLeftMargin(0.13);
	cc->cd(1)->SetRightMargin(0.045);
	cc->cd(1)->SetBottomMargin(0.145);
	if (useylogscale)
		gPad->SetLogy();
	if (usexlogscale)
		gPad->SetLogx();
	gPad->SetTickx();
	gPad->SetTicky();
	// gPad->SetFixedAspectRatio();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *legend;
	if (histo.size() > 1)
	{
		legend = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
	}

	for (int i = 0; i < histo.size(); i++)
	{
		if (i == 0)
		{
			histo.at(i)->SetMarkerStyle(71);
		}
		else if (i == 1)
		{
			histo.at(i)->SetMarkerStyle(72);
		}
		else
		{
			histo.at(i)->SetMarkerStyle(73 + i);
		}
		histo.at(i)->SetMarkerSize(1.85);
		if (i == 0)
		{
			histo.at(i)->SetLineColor(DSL_colors[i]);
		}
		else
		{
			histo.at(i)->SetLineColor(DSL_colors[i + 3]);
		}
		if (i == 0)
		{
			histo.at(i)->SetMarkerColor(DSL_colors[i]);
		}
		else
		{
			histo.at(i)->SetMarkerColor(DSL_colors[i + 3]);
		}
		histo.at(i)->SetLineWidth(3);
		histo.at(i)->SetLineStyle(1);
		if (i == 0)
		{
			histo.at(i)->Draw(Form("%s", histstyle.c_str()));
			if (h_ymin != 0 && h_ymax != 0)
				histo.at(i)->GetYaxis()->SetRangeUser(h_ymin, h_ymax);
			if (h_xmin != 0 && h_xmax != 0)
				histo.at(i)->GetXaxis()->SetRangeUser(h_xmin, h_xmax);
			// histo.at(i)->GetYaxis()->SetMaxDigits(2);
			histo.at(i)->GetYaxis()->SetTitle(Form("#scale[1.1]{%s}", ylabel.c_str()));
			histo.at(i)->GetXaxis()->SetTitle(Form("#scale[1.1]{%s}", xlabel.c_str()));
			histo.at(i)->GetXaxis()->CenterTitle(1);
			histo.at(i)->GetYaxis()->CenterTitle(1);
			histo.at(i)->GetXaxis()->SetTickSize(0.03);
			histo.at(i)->GetYaxis()->SetTickSize(0.03);
			histo.at(i)->GetXaxis()->SetLabelSize(0.045);
			histo.at(i)->GetYaxis()->SetLabelSize(0.045);
			histo.at(i)->GetYaxis()->SetLabelOffset(0.0045);
			histo.at(i)->GetXaxis()->SetLabelOffset(0.0045);
			histo.at(i)->GetYaxis()->SetTitleOffset(1.7);
			histo.at(i)->GetXaxis()->SetTitleOffset(1.5);
		}

		if (i > 0)
		{
			histo.at(i)->Draw(Form("%ssame", histstyle.c_str()));
		}

		if (histo.size() > 1)
		{
			legend->SetBorderSize(0);
			legend->SetTextSize(0.035);
			legend->SetTextFont(42);
			if (histstyle == "E" || histstyle == "e")
				legend->AddEntry(histo.at(i), Form("%s", histo_legend[i].c_str()), "P");
			if (histstyle == "HIST" || histstyle == "hist")
				legend->AddEntry(histo.at(i), Form("%s", histo_legend[i].c_str()), "L");
		}
	}

	if (histo.size() > 1)
	{
		legend->Draw();
	}

	TLatex *texCMS = new TLatex;
	texCMS->SetNDC();
	texCMS->SetTextFont(42);
	texCMS->SetTextSize(0.05);
	texCMS->SetTextAlign(12);
	texCMS->DrawLatex(rnd_x, rnd_y, Form("#scale[0.8]{%s}", randomtext.c_str()));
	texCMS->DrawLatex(pt_x, pt_y, Form("#scale[0.8]{%s}", ptrange.c_str()));
	texCMS->DrawLatex(mult_x, mult_y, Form("#scale[0.8]{%s}", multrange.c_str()));

	if (!insideofcanvas)
		texCMS->DrawLatex(0.13, 0.93, Form("#scale[0.8]{ #it{%s}}", cmstext.c_str()));
	if (insideofcanvas)
		texCMS->DrawLatex(0.16, 0.83, Form("#scale[0.8]{#splitline{}{#it{%s}}}", cmstext.c_str()));
	if (!insideofcanvas)
		texCMS->DrawLatex(0.52, 0.93, Form("#scale[0.85]{%s}", cmssystem.c_str()));
	if (insideofcanvas)
		texCMS->DrawLatex(0.52, 0.93, Form("#scale[0.85]{%s}", cmssystem.c_str()));
	cc->Update();
	// cc->Print(Form("results/%s.pdf", outfilename.c_str()));
	cc->Print(Form("results/%s.png", outfilename.c_str()));
	cc->Update();
}

TH1D *normalization(TH1D *input)
{
	input->Scale(1. / (input->Integral("width")));
	return input;
}
void histo_plots_with_ratios(std::vector<TH1 *> histo, string histstyle, double h_xmin, double h_xmax, double h_ymin, double h_ymax, bool usexlogscale, bool useylogscale, string xlabel, string ylabel, string cmstext, bool insideofcanvas, string cmssystem, string randomtext, double rnd_x, double rnd_y, string ptrange, double pt_x, double pt_y, string multrange, double mult_x, double mult_y, std::vector<string> histo_legend, double leg_xmin, double leg_ymin, double leg_xmax, double leg_ymax, double lowratio, double highratio, string outfilename)
{

	/*
	histo -> vector of histograms
	histstyle -> for plots --> Use "E" for symbols (circles, squares, ...) and "HIST" for histograms
	h_xmin, h_xmax, h_ymin, h_ymax -> ranges of x and y labels
	usexlogscale and useylogscale -> logscale in x and/or y --> true means logscale on
	cmstext -> this is the "CMS Preliminary" it depends on the status of the analysis
		--> we usualy start with Preliminary -> if needed to be shown on thesis or internal conference use Work In Progress
		--> at the end we do not need to use this (final results/plots), so just set ""
	insideofcanvas --> text inside or outside the plot. You decide what is better.
	cmssystem --> collision, energy and luminosity. I prefer like: pPb 186 nb^{-1} (8.16 TeV)
	randomtext, rnd_x and rnd_y -> some text inside of the canvas, for example MC model or so. If you do not want just keep "". x and y represents the positions in the plot.
	ptrange, pt_x and pt_y -> some text for pT range (if needed) inside of the canvas. x and y represents the positions in the plot.
	multrange mult_x and  mult_y -> similar to pT but for multiplicity
	histo_legend -> vector with name of histograms to include in the legend
	leg_xmin, leg_ymin, leg_xmax, leg_ymax -> legend position
	lowratio, highratio -> min and max values show in the ratio plot
	outfilename -> output file name
	*/

	TCanvas *cc = new TCanvas("cc", "cc", 800, 700);
	cc->cd(1);
	cc->cd(1)->SetTopMargin(0.1);
	cc->cd(1)->SetLeftMargin(0.13);
	cc->cd(1)->SetRightMargin(0.045);
	cc->cd(1)->SetBottomMargin(0.145);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetFixedAspectRatio();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TPad *pad1 = new TPad("pad1", "top pad", 0.0, 0.25, 0.98, 1.0);
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetFixedAspectRatio();
	TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.005, 0.98, 0.25);
	pad2->SetTickx();
	pad2->SetTicky();
	pad2->SetFixedAspectRatio();
	pad1->SetTopMargin(0.1);
	pad1->SetBottomMargin(0.0);
	pad1->SetRightMargin(0.03);
	pad1->SetLeftMargin(0.15);
	pad1->Draw();
	// ad2->SetFrameFillStyle(4000);
	pad2->SetTopMargin(0.0);
	pad2->SetBottomMargin(0.4);
	pad2->SetRightMargin(0.03);
	pad2->SetLeftMargin(0.15);
	pad2->Draw();

	pad1->cd(); // plots
	if (useylogscale)
		gPad->SetLogy();
	if (usexlogscale)
		gPad->SetLogx();
	TLegend *legend;
	std::vector<TH1 *> histoclone;
	if (histo.size() > 1)
	{
		legend = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
	}

	for (int i = 0; i < histo.size(); i++)
	{
		if (i == 0)
		{
			histo.at(i)->SetMarkerStyle(71);
		}
		else if (i == 1)
		{
			histo.at(i)->SetMarkerStyle(72);
		}
		else
		{
			histo.at(i)->SetMarkerStyle(73 + i);
		}
		histo.at(i)->SetMarkerSize(1.85);
		if (i == 0)
		{
			histo.at(i)->SetLineColor(DSL_colors[i]);
		}
		else
		{
			histo.at(i)->SetLineColor(DSL_colors[i + 3]);
		}
		if (i == 0)
		{
			histo.at(i)->SetMarkerColor(DSL_colors[i]);
		}
		else
		{
			histo.at(i)->SetMarkerColor(DSL_colors[i + 3]);
		}
		histo.at(i)->SetLineWidth(3);
		histo.at(i)->SetLineStyle(1);
		if (i == 0)
		{
			histo.at(i)->Draw(Form("%s", histstyle.c_str()));
			if (h_ymin != 0 && h_ymax != 0)
				histo.at(i)->GetYaxis()->SetRangeUser(h_ymin, h_ymax);
			if (h_xmin != 0 && h_xmax != 0)
				histo.at(i)->GetXaxis()->SetRangeUser(h_xmin, h_xmax);
			histo.at(i)->GetYaxis()->SetMaxDigits(2);
			histo.at(i)->GetYaxis()->SetTitle(Form("#scale[1.3]{%s}", ylabel.c_str()));
			histo.at(i)->GetXaxis()->SetTitle(Form("#scale[1.3]{%s}", xlabel.c_str()));
			histo.at(i)->GetXaxis()->CenterTitle(1);
			histo.at(i)->GetYaxis()->CenterTitle(1);
			histo.at(i)->GetXaxis()->SetTickSize(0.035);
			histo.at(i)->GetYaxis()->SetTickSize(0.025);
			histo.at(i)->GetXaxis()->SetLabelSize(0.057);
			histo.at(i)->GetYaxis()->SetLabelSize(0.057);
			histo.at(i)->GetYaxis()->SetLabelOffset(0.004);
			histo.at(i)->GetXaxis()->SetLabelOffset(0.0045);
			histo.at(i)->GetXaxis()->SetTitleSize(0.045);
			histo.at(i)->GetYaxis()->SetTitleSize(0.043);
			histo.at(i)->GetYaxis()->SetTitleOffset(1.6);
			histo.at(i)->GetXaxis()->SetTitleOffset(1.5);
		}

		if (i > 0)
		{
			histo.at(i)->Draw(Form("%ssame", histstyle.c_str()));
		}

		if (histo.size() > 1)
		{
			legend->SetBorderSize(0);
			legend->SetTextSize(0.045);
			legend->SetTextFont(42);
			if (histstyle == "E" || histstyle == "e")
				legend->AddEntry(histo.at(i), Form("%s", histo_legend[i].c_str()), "P");
			if (histstyle == "HIST" || histstyle == "hist")
				legend->AddEntry(histo.at(i), Form("%s", histo_legend[i].c_str()), "L");
		}
		TH1 *histogram = (TH1 *)histo.at(i)->Clone();
		histoclone.push_back(histogram);
	}

	if (histo.size() > 1)
	{
		legend->Draw();
	}

	TLatex *texCMS = new TLatex;
	texCMS->SetNDC();
	texCMS->SetTextFont(42);
	texCMS->SetTextSize(0.05);
	texCMS->SetTextAlign(12);
	texCMS->DrawLatex(rnd_x, rnd_y, Form("#scale[1]{%s}", randomtext.c_str()));
	texCMS->DrawLatex(pt_x, pt_y, Form("#scale[1]{%s}", ptrange.c_str()));
	texCMS->DrawLatex(mult_x, mult_y, Form("#scale[1]{%s}", multrange.c_str()));

	if (!insideofcanvas)
		texCMS->DrawLatex(0.17, 0.94, Form("#scale[1.1]{ #it{%s}}", cmstext.c_str()));
	if (insideofcanvas)
		texCMS->DrawLatex(0.2, 0.83, Form("#scale[1.1]{#splitline{}{#it{%s}}}", cmstext.c_str()));
	if (!insideofcanvas)
		texCMS->DrawLatex(0.59, 0.94, Form("#scale[1.1]{%s}", cmssystem.c_str()));
	if (insideofcanvas)
		texCMS->DrawLatex(0.59, 0.94, Form("#scale[1.1]{%s}", cmssystem.c_str()));

	pad2->cd(); // ratio plots

	for (int j = 1; j < histoclone.size(); j++)
	{
		histoclone.at(j)->Divide(histoclone.at(j), histoclone.at(0), 1, 1, "B");
		if (j == 1)
		{
			histoclone.at(j)->GetYaxis()->SetTitle("#scale[4]{Ratio}");
			histoclone.at(j)->GetXaxis()->SetTitle(Form("#scale[4]{%s}", xlabel.c_str()));
			histoclone.at(j)->GetXaxis()->CenterTitle(1);
			histoclone.at(j)->GetYaxis()->CenterTitle(1);
			histoclone.at(j)->GetXaxis()->SetTickSize(0.11);
			histoclone.at(j)->GetYaxis()->SetTickSize(0.03);
			histoclone.at(j)->GetYaxis()->SetNdivisions(404);
			histoclone.at(j)->GetXaxis()->SetNdivisions(810);
			histoclone.at(j)->GetXaxis()->SetLabelSize(0.17);
			histoclone.at(j)->GetYaxis()->SetLabelSize(0.17);
			histoclone.at(j)->GetYaxis()->SetLabelOffset(0.005);
			histoclone.at(j)->GetXaxis()->SetLabelOffset(0.006);
			histoclone.at(j)->GetYaxis()->SetTitleOffset(1.8);
			histoclone.at(j)->GetXaxis()->SetTitleOffset(4);
			histoclone.at(j)->GetXaxis()->SetTitleSize(0.045);
			histoclone.at(j)->GetYaxis()->SetTitleSize(0.043);
			histoclone.at(j)->GetYaxis()->SetRangeUser(lowratio, highratio);
			if (h_xmin != 0 && h_xmax != 0)
				histoclone.at(j)->GetXaxis()->SetRangeUser(h_xmin, h_xmax);
		}
		if (j == 1)
		{
			histoclone.at(j)->Draw("E");
		}
		else
		{
			histoclone.at(j)->Draw("ESAME");
		}
		// if(histoclone.size()>2){legend->Draw();}

		

	}

	TLine *lineatone = new TLine(h_xmin, 1, h_xmax, 1);
	lineatone->SetLineStyle(2);
	lineatone->Draw("SAME");

	cc->Update();
	// cc->Print(Form("results/%s.pdf", outfilename.c_str()));
	cc->Print(Form("results/%s.png", outfilename.c_str()));
	cc->Update();
}

void histo_plots_with_difference(std::vector<TH1 *> histo, string histstyle, double h_xmin, double h_xmax, double h_ymin, double h_ymax, bool usexlogscale, bool useylogscale, string xlabel, string ylabel, string cmstext, bool insideofcanvas, string cmssystem, string randomtext, double rnd_x, double rnd_y, string ptrange, double pt_x, double pt_y, string multrange, double mult_x, double mult_y, std::vector<string> histo_legend, double leg_xmin, double leg_ymin, double leg_xmax, double leg_ymax, double lowratio, double highratio, string outfilename)
{

	/*
	histo -> vector of histograms
	histstyle -> for plots --> Use "E" for symbols (circles, squares, ...) and "HIST" for histograms
	h_xmin, h_xmax, h_ymin, h_ymax -> ranges of x and y labels
	usexlogscale and useylogscale -> logscale in x and/or y --> true means logscale on
	cmstext -> this is the "CMS Preliminary" it depends on the status of the analysis
		--> we usualy start with Preliminary -> if needed to be shown on thesis or internal conference use Work In Progress
		--> at the end we do not need to use this (final results/plots), so just set ""
	insideofcanvas --> text inside or outside the plot. You decide what is better.
	cmssystem --> collision, energy and luminosity. I prefer like: pPb 186 nb^{-1} (8.16 TeV)
	randomtext, rnd_x and rnd_y -> some text inside of the canvas, for example MC model or so. If you do not want just keep "". x and y represents the positions in the plot.
	ptrange, pt_x and pt_y -> some text for pT range (if needed) inside of the canvas. x and y represents the positions in the plot.
	multrange mult_x and  mult_y -> similar to pT but for multiplicity
	histo_legend -> vector with name of histograms to include in the legend
	leg_xmin, leg_ymin, leg_xmax, leg_ymax -> legend position
	lowratio, highratio -> min and max values show in the ratio plot
	outfilename -> output file name
	*/

	TCanvas *cc = new TCanvas("cc", "cc", 800, 700);
	cc->cd(1);
	cc->cd(1)->SetTopMargin(0.1);
	cc->cd(1)->SetLeftMargin(0.13);
	cc->cd(1)->SetRightMargin(0.045);
	cc->cd(1)->SetBottomMargin(0.145);
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetFixedAspectRatio();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TPad *pad1 = new TPad("pad1", "top pad", 0.0, 0.25, 0.98, 1.0);
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetFixedAspectRatio();
	TPad *pad2 = new TPad("pad2", "bottom pad", 0.0, 0.005, 0.98, 0.25);
	pad2->SetTickx();
	pad2->SetTicky();
	pad2->SetFixedAspectRatio();
	pad1->SetTopMargin(0.1);
	pad1->SetBottomMargin(0.0);
	pad1->SetRightMargin(0.03);
	pad1->SetLeftMargin(0.15);
	pad1->Draw();
	// ad2->SetFrameFillStyle(4000);
	pad2->SetTopMargin(0.0);
	pad2->SetBottomMargin(0.4);
	pad2->SetRightMargin(0.03);
	pad2->SetLeftMargin(0.15);
	pad2->Draw();

	pad1->cd(); // plots
	if (useylogscale)
		gPad->SetLogy();
	if (usexlogscale)
		gPad->SetLogx();
	TLegend *legend;
	std::vector<TH1 *> histoclone;
	if (histo.size() > 1)
	{
		legend = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
	}

	for (int i = 0; i < histo.size(); i++)
	{
		if (i == 0)
		{
			histo.at(i)->SetMarkerStyle(71);
		}
		else if (i == 1)
		{
			histo.at(i)->SetMarkerStyle(72);
		}
		else
		{
			histo.at(i)->SetMarkerStyle(73 + i);
		}
		histo.at(i)->SetMarkerSize(1.85);
		if (i == 0)
		{
			histo.at(i)->SetLineColor(DSL_colors[i]);
		}
		else
		{
			histo.at(i)->SetLineColor(DSL_colors[i + 3]);
		}
		if (i == 0)
		{
			histo.at(i)->SetMarkerColor(DSL_colors[i]);
		}
		else
		{
			histo.at(i)->SetMarkerColor(DSL_colors[i + 3]);
		}
		histo.at(i)->SetLineWidth(3);
		histo.at(i)->SetLineStyle(1);
		if (i == 0)
		{
			histo.at(i)->Draw(Form("%s", histstyle.c_str()));
			if (h_ymin != 0 && h_ymax != 0)
				histo.at(i)->GetYaxis()->SetRangeUser(h_ymin, h_ymax);
			if (h_xmin != 0 && h_xmax != 0)
				histo.at(i)->GetXaxis()->SetRangeUser(h_xmin, h_xmax);
			histo.at(i)->GetYaxis()->SetMaxDigits(2);
			histo.at(i)->GetYaxis()->SetTitle(Form("#scale[1.3]{%s}", ylabel.c_str()));
			histo.at(i)->GetXaxis()->SetTitle(Form("#scale[1.3]{%s}", xlabel.c_str()));
			histo.at(i)->GetXaxis()->CenterTitle(1);
			histo.at(i)->GetYaxis()->CenterTitle(1);
			histo.at(i)->GetXaxis()->SetTickSize(0.035);
			histo.at(i)->GetYaxis()->SetTickSize(0.025);
			histo.at(i)->GetXaxis()->SetLabelSize(0.057);
			histo.at(i)->GetYaxis()->SetLabelSize(0.057);
			histo.at(i)->GetYaxis()->SetLabelOffset(0.004);
			histo.at(i)->GetXaxis()->SetLabelOffset(0.0045);
			histo.at(i)->GetXaxis()->SetTitleSize(0.045);
			histo.at(i)->GetYaxis()->SetTitleSize(0.043);
			histo.at(i)->GetYaxis()->SetTitleOffset(1.6);
			histo.at(i)->GetXaxis()->SetTitleOffset(1.5);
		}

		if (i > 0)
		{
			histo.at(i)->Draw(Form("%ssame", histstyle.c_str()));
		}

		if (histo.size() > 1)
		{
			legend->SetBorderSize(0);
			legend->SetTextSize(0.045);
			legend->SetTextFont(42);
			if (histstyle == "E" || histstyle == "e")
				legend->AddEntry(histo.at(i), Form("%s", histo_legend[i].c_str()), "P");
			if (histstyle == "HIST" || histstyle == "hist")
				legend->AddEntry(histo.at(i), Form("%s", histo_legend[i].c_str()), "L");
		}
		TH1 *histogram = (TH1 *)histo.at(i)->Clone();
		histoclone.push_back(histogram);
	}

	if (histo.size() > 1)
	{
		legend->Draw();
	}

	TLatex *texCMS = new TLatex;
	texCMS->SetNDC();
	texCMS->SetTextFont(42);
	texCMS->SetTextSize(0.05);
	texCMS->SetTextAlign(12);
	texCMS->DrawLatex(rnd_x, rnd_y, Form("#scale[1]{%s}", randomtext.c_str()));
	texCMS->DrawLatex(pt_x, pt_y, Form("#scale[1]{%s}", ptrange.c_str()));
	texCMS->DrawLatex(mult_x, mult_y, Form("#scale[1]{%s}", multrange.c_str()));

	if (!insideofcanvas)
		texCMS->DrawLatex(0.17, 0.94, Form("#scale[1.1]{ #it{%s}}", cmstext.c_str()));
	if (insideofcanvas)
		texCMS->DrawLatex(0.2, 0.83, Form("#scale[1.1]{#splitline{}{#it{%s}}}", cmstext.c_str()));
	if (!insideofcanvas)
		texCMS->DrawLatex(0.59, 0.94, Form("#scale[1.1]{%s}", cmssystem.c_str()));
	if (insideofcanvas)
		texCMS->DrawLatex(0.59, 0.94, Form("#scale[1.1]{%s}", cmssystem.c_str()));

	pad2->cd(); // ratio plots

	for (int j = 1; j < histoclone.size(); j++)
	{
		histoclone.at(j)->Add(histoclone.at(0), -1);
		if (j == 1)
		{
			histoclone.at(j)->GetYaxis()->SetTitle("#scale[4]{Diff.}");
			histoclone.at(j)->GetXaxis()->SetTitle(Form("#scale[4]{%s}", xlabel.c_str()));
			histoclone.at(j)->GetXaxis()->CenterTitle(1);
			histoclone.at(j)->GetYaxis()->CenterTitle(1);
			histoclone.at(j)->GetXaxis()->SetTickSize(0.11);
			histoclone.at(j)->GetYaxis()->SetTickSize(0.03);
			histoclone.at(j)->GetYaxis()->SetNdivisions(404);
			histoclone.at(j)->GetXaxis()->SetNdivisions(810);
			histoclone.at(j)->GetXaxis()->SetLabelSize(0.17);
			histoclone.at(j)->GetYaxis()->SetLabelSize(0.17);
			histoclone.at(j)->GetYaxis()->SetLabelOffset(0.005);
			histoclone.at(j)->GetXaxis()->SetLabelOffset(0.006);
			histoclone.at(j)->GetYaxis()->SetTitleOffset(1.6);
			histoclone.at(j)->GetXaxis()->SetTitleOffset(4);
			histoclone.at(j)->GetXaxis()->SetTitleSize(0.045);
			histoclone.at(j)->GetYaxis()->SetTitleSize(0.043);
			histoclone.at(j)->GetYaxis()->SetRangeUser(lowratio, highratio);
			if (h_xmin != 0 && h_xmax != 0)
				histoclone.at(j)->GetXaxis()->SetRangeUser(h_xmin, h_xmax);
		}
		if (j == 1)
		{
			histoclone.at(j)->Draw("E");
		}
		else
		{
			histoclone.at(j)->Draw("ESAME");
		}
		// if(histoclone.size()>2){legend->Draw();}
	}

	TLine *lineatone = new TLine(h_xmin, 0, h_xmax, 0);
	lineatone->SetLineStyle(2);
	lineatone->Draw("SAME");

	cc->Update();
	// cc->Print(Form("results/%s.pdf", outfilename.c_str()));
	cc->Print(Form("results/%s.png", outfilename.c_str()));
	cc->Update();
}

void histo_2D_corr_plots(TH2 *histo, string histstyle, double h_xmin, double h_xmax, double h_ymin, double h_ymax, double h_zmin, double h_zmax, bool usexlogscale, bool useylogscale, bool usezlogscale, string xlabel, string ylabel, string zlabel, string cmstext, string cmssystem, string randomtext, double rnd_x, double rnd_y, string ptrange, double pt_x, double pt_y, string multrange, double mult_x, double mult_y, string outfilename)
{

	TCanvas *cc = new TCanvas("cc", "cc", 900, 700);
	cc->cd(1);
	cc->cd(1)->SetTopMargin(0.1);
	cc->cd(1)->SetLeftMargin(0.1);
	cc->cd(1)->SetRightMargin(0.17);
	cc->cd(1)->SetBottomMargin(0.145);
	if (usezlogscale)
		gPad->SetLogz();
	if (useylogscale)
		gPad->SetLogy();
	if (usexlogscale)
		gPad->SetLogx();
	gPad->SetTickx();
	gPad->SetTicky();
	// gPad->SetTickz();
	gPad->SetFixedAspectRatio();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	histo->Draw(Form("%s", histstyle.c_str()));
	if (h_ymin != 0 && h_ymax != 0)
		histo->GetYaxis()->SetRangeUser(h_ymin, h_ymax);
	if (h_zmin != 0 && h_zmax != 0)
		histo->GetZaxis()->SetRangeUser(h_zmin, h_zmax);
	if (h_xmin != 0 && h_xmax != 0)
		histo->GetXaxis()->SetRangeUser(h_xmin, h_xmax);
	histo->GetZaxis()->SetMaxDigits(2);
	histo->GetZaxis()->SetTitle(Form("#scale[1.2]{%s}", zlabel.c_str()));
	histo->GetYaxis()->SetTitle(Form("#scale[1.2]{%s}", ylabel.c_str()));
	histo->GetXaxis()->SetTitle(Form("#scale[1.2]{%s}", xlabel.c_str()));
	histo->GetXaxis()->CenterTitle(1);
	histo->GetYaxis()->CenterTitle(1);
	histo->GetZaxis()->CenterTitle(1);
	histo->GetXaxis()->SetTickSize(0.04);
	histo->GetYaxis()->SetTickSize(0.04);
	histo->GetZaxis()->SetTickSize(0.04);
	histo->GetXaxis()->SetLabelSize(0.04);
	histo->GetYaxis()->SetLabelSize(0.04);
	histo->GetZaxis()->SetLabelSize(0.04);
	histo->GetYaxis()->SetLabelOffset(0.004);
	histo->GetXaxis()->SetLabelOffset(0.004);
	histo->GetZaxis()->SetLabelOffset(0.004);
	histo->GetYaxis()->SetTitleOffset(1.3);
	histo->GetXaxis()->SetTitleOffset(1.4);
	histo->GetZaxis()->SetTitleOffset(1.7);
	histo->GetXaxis()->SetNdivisions(507);
	histo->GetYaxis()->SetNdivisions(507);
	histo->GetZaxis()->SetNdivisions(507);

	TLatex *texCMS = new TLatex;
	texCMS->SetNDC();
	texCMS->SetTextFont(42);
	texCMS->SetTextSize(0.05);
	texCMS->SetTextAlign(12);
	texCMS->DrawLatex(rnd_x, rnd_y, Form("#scale[0.8]{%s}", randomtext.c_str()));
	texCMS->DrawLatex(pt_x, pt_y, Form("#scale[0.8]{%s}", ptrange.c_str()));
	texCMS->DrawLatex(mult_x, mult_y, Form("#scale[0.8]{%s}", multrange.c_str()));
	texCMS->DrawLatex(0.13, 0.93, Form("#scale[0.8]{ #it{%s}}", cmstext.c_str()));
	texCMS->DrawLatex(0.52, 0.93, Form("#scale[0.85]{%s}", cmssystem.c_str()));

	cc->Update();
	cc->Print(Form("results/%s.png", outfilename.c_str()));
	//cc->Update();
	//cc->Print(Form("results/%s.pdf", outfilename.c_str()));
	cc->Update();
}

void histo_3D_corr_plots(TH2 *histo, string histstyle, double h_xmin, double h_xmax, double h_ymin, double h_ymax, double h_zmin, double h_zmax, bool usexlogscale, bool useylogscale, bool usezlogscale, string xlabel, string ylabel, string zlabel, string cmstext, string cmssystem, string randomtext, double rnd_x, double rnd_y, string ptrange, double pt_x, double pt_y, string multrange, double mult_x, double mult_y, string outfilename)
{

	TCanvas *cc = new TCanvas("cc", "cc", 700, 700);
	cc->cd(1);
	cc->cd(1)->SetTopMargin(0.1);
	cc->cd(1)->SetLeftMargin(0.15);
	cc->cd(1)->SetRightMargin(0.05);
	cc->cd(1)->SetBottomMargin(0.145);
	if (usezlogscale)
		gPad->SetLogz();
	if (useylogscale)
		gPad->SetLogy();
	if (usexlogscale)
		gPad->SetLogx();
	gPad->SetTickx();
	gPad->SetTicky();
	gPad->SetGrid(0);
	gPad->SetFixedAspectRatio();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	histo->Draw(Form("%s", histstyle.c_str()));
	if (h_ymin != 0 && h_ymax != 0)
		histo->GetYaxis()->SetRangeUser(h_ymin, h_ymax);
	if (h_zmin != 0 && h_zmax != 0)
		histo->GetZaxis()->SetRangeUser(h_zmin, h_zmax);
	if (h_xmin != 0 && h_xmax != 0)
		histo->GetXaxis()->SetRangeUser(h_xmin, h_xmax);
	histo->GetZaxis()->SetMaxDigits(2);
	histo->GetZaxis()->SetTitle(Form("#scale[1.3]{%s}", zlabel.c_str()));
	histo->GetYaxis()->SetTitle(Form("#scale[1.3]{%s}", ylabel.c_str()));
	histo->GetXaxis()->SetTitle(Form("#scale[1.3]{%s}", xlabel.c_str()));
	histo->GetXaxis()->CenterTitle(1);
	histo->GetYaxis()->CenterTitle(1);
	histo->GetZaxis()->CenterTitle(1);
	histo->GetXaxis()->SetTickSize(0.04);
	histo->GetYaxis()->SetTickSize(0.04);
	histo->GetZaxis()->SetTickSize(0.04);
	histo->GetXaxis()->SetLabelSize(0.04);
	histo->GetYaxis()->SetLabelSize(0.04);
	histo->GetZaxis()->SetLabelSize(0.04);
	histo->GetYaxis()->SetLabelOffset(0.004);
	histo->GetXaxis()->SetLabelOffset(0.004);
	histo->GetZaxis()->SetLabelOffset(0.004);
	histo->GetYaxis()->SetTitleOffset(2);
	histo->GetXaxis()->SetTitleOffset(2);
	histo->GetZaxis()->SetTitleOffset(2);
	histo->GetXaxis()->SetNdivisions(507);
	histo->GetYaxis()->SetNdivisions(507);
	histo->GetZaxis()->SetNdivisions(507);

	TLatex *texCMS = new TLatex;
	texCMS = new TLatex(rnd_x, rnd_y, Form("#scale[0.8]{%s}", randomtext.c_str()));
	texCMS->SetTextFont(42);
	texCMS->SetTextSize(0.05);
	texCMS->SetTextAlign(12);
	texCMS->SetTextAngle(13);
	texCMS->Draw();

	TLatex *texpt = new TLatex;
	texpt = new TLatex(pt_x, pt_y, Form("#scale[0.8]{%s}", ptrange.c_str()));
	texpt->SetTextFont(42);
	texpt->SetTextSize(0.05);
	texpt->SetTextAlign(12);
	texpt->SetTextAngle(-35);
	texpt->Draw();

	TLatex *texmult = new TLatex;
	texmult = new TLatex(mult_x, mult_y, Form("#scale[0.8]{%s}", multrange.c_str()));
	texmult->SetTextFont(42);
	texmult->SetTextSize(0.05);
	texmult->SetTextAlign(12);
	texmult->SetTextAngle(-35);
	texmult->Draw();

	TLatex *tex_sys = new TLatex;
	tex_sys = new TLatex(0.26, 0.96, Form("#scale[0.85]{%s}", cmssystem.c_str()));
	tex_sys->SetTextFont(42);
	tex_sys->SetTextSize(0.05);
	tex_sys->SetTextAlign(12);
	tex_sys->SetTextAngle(-35);
	tex_sys->Draw();

	TLatex *tex_colab = new TLatex;
	tex_colab = new TLatex(-0.5, 0.77, Form("#scale[0.8]{ #it{%s}}", cmstext.c_str()));
	tex_colab->SetTextFont(42);
	tex_colab->SetTextSize(0.05);
	tex_colab->SetTextAlign(12);
	tex_colab->SetTextAngle(13);
	tex_colab->Draw();

	// gPad->GetView()->RotateView(90,0.0);
	// gPad->Modified();
	cc->Update();
	// cc->Print(Form("results/%s.pdf", outfilename.c_str()));
	cc->Print(Form("results/%s.png", outfilename.c_str()));
	cc->Update();
}


void histo_plots_noline(std::vector<TH1 *> histo, string histstyle, double h_xmin, double h_xmax, double h_ymin, double h_ymax, bool usexlogscale, bool useylogscale, string xlabel, string ylabel, string cmstext, bool insideofcanvas, string cmssystem, string randomtext, double rnd_x, double rnd_y, string ptrange, double pt_x, double pt_y, string multrange, double mult_x, double mult_y, std::vector<string> histo_legend, double leg_xmin, double leg_ymin, double leg_xmax, double leg_ymax, string outfilename)
{

	/*
	histo -> vector of histograms
	histstyle -> for plots --> Use "E" for symbols (circles, squares, ...) and "HIST" for histograms
	h_xmin, h_xmax, h_ymin, h_ymax -> ranges of x and y labels
	usexlogscale and useylogscale -> logscale in x and/or y --> true means logscale on
	cmstext -> this is the "CMS Preliminary" it depends on the status of the analysis
		--> we usualy start with Preliminary -> if needed to be shown on thesis or internal conference use Work In Progress
		--> at the end we do not need to use this (final results/plots), so just set ""
	insideofcanvas --> text inside or outside the plot. You decide what is better.
	cmssystem --> collision, energy and luminosity. I prefer like: pPb 186 nb^{-1} (8.16 TeV)
	randomtext, rnd_x and rnd_y -> some text inside of the canvas, for example MC model or so. If you do not want just keep "". x and y represents the positions in the plot.
	ptrange, pt_x and pt_y -> some text for pT range (if needed) inside of the canvas. x and y represents the positions in the plot.
	multrange mult_x and  mult_y -> similar to pT but for multiplicity
	histo_legend -> vector with name of histograms to include in the legend
	leg_xmin, leg_ymin, leg_xmax, leg_ymax -> legend position
	outfilename -> output file name
	*/

	TCanvas *cc = new TCanvas("cc", "cc", 900, 900);
	cc->cd(1);
	cc->cd(1)->SetTopMargin(0.1);
	cc->cd(1)->SetLeftMargin(0.17);
	cc->cd(1)->SetRightMargin(0.045);
	cc->cd(1)->SetBottomMargin(0.145);
	if (useylogscale)
		gPad->SetLogy();
	if (usexlogscale)
		gPad->SetLogx();
	gPad->SetTickx();
	gPad->SetTicky();
	// gPad->SetFixedAspectRatio();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TLegend *legend;
	if (histo.size() > 1)
	{
		legend = new TLegend(leg_xmin, leg_ymin, leg_xmax, leg_ymax);
	}

	for (int i = 0; i < histo.size(); i++)
	{
		if (i == 0)
		{
			histo.at(i)->SetMarkerStyle(71);
		}
		else if (i == 1)
		{
			histo.at(i)->SetMarkerStyle(72);
		}
		else
		{
			histo.at(i)->SetMarkerStyle(73 + i);
		}
		histo.at(i)->SetMarkerSize(1.85);
		if (i == 0)
		{
			histo.at(i)->SetLineColor(DSL_colors[i]);
		}
		else
		{
			histo.at(i)->SetLineColor(DSL_colors[i + 3]);
		}
		if (i == 0)
		{
			histo.at(i)->SetMarkerColor(DSL_colors[i]);
		}
		else
		{
			histo.at(i)->SetMarkerColor(DSL_colors[i + 3]);
		}
		histo.at(i)->SetLineWidth(3);
		histo.at(i)->SetLineStyle(1);
		if (i == 0)
		{
			histo.at(i)->Draw(Form("%s", histstyle.c_str()));
			if (h_ymin != 0 && h_ymax != 0)
				histo.at(i)->GetYaxis()->SetRangeUser(h_ymin, h_ymax);
			if (h_xmin != 0 && h_xmax != 0)
				histo.at(i)->GetXaxis()->SetRangeUser(h_xmin, h_xmax);
			// histo.at(i)->GetYaxis()->SetMaxDigits(2);
			histo.at(i)->GetYaxis()->SetTitle(Form("#scale[1.]{%s}", ylabel.c_str()));
			histo.at(i)->GetXaxis()->SetTitle(Form("#scale[1.]{%s}", xlabel.c_str()));
			histo.at(i)->GetXaxis()->CenterTitle(1);
			histo.at(i)->GetYaxis()->CenterTitle(1);
			histo.at(i)->GetXaxis()->SetTickSize(0.03);
			histo.at(i)->GetYaxis()->SetTickSize(0.03);
			histo.at(i)->GetXaxis()->SetLabelSize(0.04);
			histo.at(i)->GetYaxis()->SetLabelSize(0.04);
			histo.at(i)->GetYaxis()->SetLabelOffset(0.0045);
			histo.at(i)->GetXaxis()->SetLabelOffset(0.0045);
			histo.at(i)->GetYaxis()->SetTitleOffset(1.8);
			histo.at(i)->GetXaxis()->SetTitleOffset(1.5);
		}

		if (i > 0)
		{
			histo.at(i)->Draw(Form("%ssame", histstyle.c_str()));
		}

		if (histo.size() > 1)
		{
			legend->SetBorderSize(0);
			legend->SetTextSize(0.035);
			legend->SetTextFont(42);
			if (histstyle == "E" || histstyle == "e")
				legend->AddEntry(histo.at(i), Form("%s", histo_legend[i].c_str()), "P");
			if (histstyle == "HIST" || histstyle == "hist")
				legend->AddEntry(histo.at(i), Form("%s", histo_legend[i].c_str()), "L");
		}
	}

	if (histo.size() > 1)
	{
		legend->Draw();
	}

	TLatex *texCMS = new TLatex;
	texCMS->SetNDC();
	texCMS->SetTextFont(42);
	texCMS->SetTextSize(0.05);
	texCMS->SetTextAlign(12);
	texCMS->DrawLatex(rnd_x, rnd_y, Form("#scale[0.8]{%s}", randomtext.c_str()));
	texCMS->DrawLatex(pt_x, pt_y, Form("#scale[0.8]{%s}", ptrange.c_str()));
	texCMS->DrawLatex(mult_x, mult_y, Form("#scale[0.8]{%s}", multrange.c_str()));
	

	if (!insideofcanvas)
		texCMS->DrawLatex(0.13, 0.93, Form("#scale[0.7]{      #it{%s}}", cmstext.c_str()));
	if (insideofcanvas)
		texCMS->DrawLatex(0.16, 0.83, Form("#scale[0.7]{     #splitline{}{#it{%s}}}", cmstext.c_str()));
	if (!insideofcanvas)
		texCMS->DrawLatex(0.52, 0.93, Form("#scale[0.7]{%s}", cmssystem.c_str()));
	if (insideofcanvas)
		texCMS->DrawLatex(0.52, 0.93, Form("#scale[0.7]{%s}", cmssystem.c_str()));


	cc->Update();
	// cc->Print(Form("results/%s.pdf", outfilename.c_str()));
	cc->Print(Form("results/%s.png", outfilename.c_str()));
	cc->Update();
}