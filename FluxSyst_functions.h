// Script with finctions needed for fluxsystematics script to run
#include "plot_comp_functions.h"
// ------------------------------------------------------------------------------------------------------------
class event {
	public:
		std::string type; // Gen, Sig, Bkg, Sel, dirt
		double E; // Energy
		double Theta; // Theta
		int nu_flav; // nuetrino flavour type for picking the right histogram to reweight
};
// ------------------------------------------------------------------------------------------------------------
// Class to store the ratios of the CV to new universe for each nu flavour
class HistWeights {
	public:
		HistWeights(std::string flav_) {flav = flav_;} // constructor to set type
		std::vector<std::vector<TH2D*>> HP; // Vector of histograms for each universe/cv
		std::vector<TH2D*> Beamline; // As above but for the beamline variations
		std::string flav; // nuetrino type, nue, nuebar, numu, numubar
}; 
// ------------------------------------------------------------------------------------------------------------
// Function to retun the right neutrino flavour ro weight by
std::string GetMode(int nu_flav){
	if  (nu_flav == 1 || nu_flav == 3) return "nue";
	else if  (nu_flav == 5 || nu_flav == 7) return "nuebar";
	else if  (nu_flav == 2 || nu_flav == 4) return "numu";
	else if  (nu_flav == 6 || nu_flav == 8) return "numubar";
	else return " ";
}
// ------------------------------------------------------------------------------------------------------------
void DivideHists(TH2D* hCV, TH2D* hUniv, TH2D* &htemp){

	// Loop over rows
	for (unsigned int i=1; i<hCV->GetNbinsX()+1; i++) { 
		// Loop over columns
		for (unsigned int j=1; j<hCV->GetNbinsY()+1; j++){ 
			if (hCV->GetBinContent(i, j) == 0) htemp->SetBinContent(i, j, 0);
			else htemp->SetBinContent(i, j, hUniv->GetBinContent(i, j) / hCV->GetBinContent(i, j) );
		}
	}
}
// ------------------------------------------------------------------------------------------------------------
// Function to integrate a histogram
double IntegrateHist(TH2D* h){
	double integral{0};
	// Loop over rows
	for (unsigned int i=1; i<h->GetNbinsX()+1; i++) { 
		// Loop over columns
		for (unsigned int j=1; j<h->GetNbinsY()+1; j++){ 
			integral+= h->GetBinContent(i, j);
		}
	}
	return integral;
}
// ------------------------------------------------------------------------------------------------------------
// Function to take the ratio of universe i with nominal and return a weight for an E, theta
double GetWeight(int universe, int index, event event, HistWeights &nue,  HistWeights &nuebar,  HistWeights &numu,  HistWeights &numubar, std::string param){
	
	double weight, xbin, ybin;

	// Get the neutrino flavour to reweigh the event type
	std::string mode = GetMode(event.nu_flav); 

	// if index is 0 then we are getting a weight from a HP universe
	if (param.find("PPFX") != std::string::npos){
		if (mode == "nue"){
			xbin = nue.HP.at(index).at(universe)->GetXaxis()->FindBin(event.E);
			ybin = nue.HP.at(index).at(universe)->GetYaxis()->FindBin(event.Theta);
			weight = nue.HP.at(index).at(universe)->GetBinContent(xbin, ybin);
		}
		else if (mode == "nuebar"){
			xbin = nuebar.HP.at(index).at(universe)->GetXaxis()->FindBin(event.E);
			ybin = nuebar.HP.at(index).at(universe)->GetYaxis()->FindBin(event.Theta);
			weight = nuebar.HP.at(index).at(universe)->GetBinContent(xbin, ybin);
		}
		else if (mode == "numu"){
			xbin = numu.HP.at(index).at(universe)->GetXaxis()->FindBin(event.E);
			ybin = numu.HP.at(index).at(universe)->GetYaxis()->FindBin(event.Theta);
			weight = numu.HP.at(index).at(universe)->GetBinContent(xbin, ybin);
		}
		else if (mode == "numubar"){
			xbin = numubar.HP.at(index).at(universe)->GetXaxis()->FindBin(event.E);
			ybin = numubar.HP.at(index).at(universe)->GetYaxis()->FindBin(event.Theta);
			weight = numubar.HP.at(index).at(universe)->GetBinContent(xbin, ybin);
		}
		else std::cout << "Unknown mode!"<< std::endl;
	
	}
	else {
		
		// CV 
		if (param == "CV") return weight = 1;
		
		// else we have a beamline variations
		else {
			if (mode == "nue"){
				xbin =   nue.Beamline.at(index-1)->GetXaxis()->FindBin(event.E);
				ybin =   nue.Beamline.at(index-1)->GetYaxis()->FindBin(event.Theta);
				weight = nue.Beamline.at(index-1)->GetBinContent(xbin, ybin);
			}
			else if (mode == "nuebar"){
				xbin =   nuebar.Beamline.at(index-1)->GetXaxis()->FindBin(event.E);
				ybin =   nuebar.Beamline.at(index-1)->GetYaxis()->FindBin(event.Theta);
				weight = nuebar.Beamline.at(index-1)->GetBinContent(xbin, ybin);
			}
			else if (mode == "numu"){
				xbin =   numu.Beamline.at(index-1)->GetXaxis()->FindBin(event.E);
				ybin =   numu.Beamline.at(index-1)->GetYaxis()->FindBin(event.Theta);
				weight = numu.Beamline.at(index-1)->GetBinContent(xbin, ybin);
			}
			else if (mode == "numubar"){
				xbin =   numubar.Beamline.at(index-1)->GetXaxis()->FindBin(event.E);
				ybin =   numubar.Beamline.at(index-1)->GetYaxis()->FindBin(event.Theta);
				weight = numubar.Beamline.at(index-1)->GetBinContent(xbin, ybin);
			}
			else std::cout << "Unknown mode!"<< std::endl;
		}
	}

	// std::cout << weight << std::endl; // DANGEROUS!!

	return weight;
}
// ------------------------------------------------------------------------------------------------------------
// Function to add the weights for each universe -- will want to overload this with another function to handle unisims
void AddWeights(std::vector<double> &N, int Universes , int index, event event, HistWeights &nue,  HistWeights &nuebar,  HistWeights &numu,  HistWeights &numubar, std::string param){
	double weight{0};

	// Initialise the size of the counter if it is the first event loop. 
	if (N.size() == 0 )  N.resize( Universes );  // Resize to number of Universes.

	// Loop over each universe
	for (unsigned int i = 0; i < Universes; i++){ 

		// Get the weight
		weight = GetWeight(i ,index, event, nue, nuebar, numu, numubar, param);

		N.at(i) += weight; // Add weight to vector of counters.

	} 
}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the interated flux
double IntegrateFlux(int universe, TFile* fCV, int index, double POTScale, std::string param){

	double flux{0};
	TH2D *hBeamline2d, *hBeamline2dnue, *hBeamline2dnuebar, *hHP2d, *hHP2dnue ,*hHP2dnuebar, *hCV2dnue, *hCV2dnuebar, *hCV2d;
	TFile *fBeamline;
	bool boolfile, boolhist;

	// if index is 0 then we are running over HP
	// get histogram from fCV for Masterweight PPFX universe i and divide out  by the cv
	if (param.find("PPFX") != std::string::npos) {

		boolhist = GetHist(fCV, hHP2dnue, Form("nue/Multisims/nue_%s_Uni_%i_AV_TPC_2D", param.c_str(), universe)); if (boolhist == false) gSystem->Exit(0);
		boolhist = GetHist(fCV, hHP2dnuebar, Form("nuebar/Multisims/nuebar_%s_Uni_%i_AV_TPC_2D", param.c_str(), universe)); if (boolhist == false) gSystem->Exit(0);
		
		hHP2d = (TH2D*) hHP2dnue->Clone("hHP2d");
		hHP2d->Add(hHP2dnuebar); // Combine the fluxes

		double xbin_th = hHP2d->GetXaxis()->FindBin( 0.75*0.2065); // find the x bin to integrate from (threshold)

		flux  = hHP2d->Integral( xbin_th, hHP2d->GetNbinsX()+1, 0, hHP2d->GetNbinsY()+1); // Integrate over whole phase space (not quite any more)
		
		flux*= (POTScale / (GetPOT(fCV)*1.0e4)); // Scale to cm2 and the DATA POT

	}
	// else we have a beamline variations indexes from 2 to N
	else {

		// CV
		if (param == "CV"){
			boolhist = GetHist(fCV, hCV2dnue, "nue/Detsmear/nue_CV_AV_TPC_2D"); if (boolhist == false) gSystem->Exit(0);
			boolhist = GetHist(fCV, hCV2dnuebar, "nuebar/Detsmear/nuebar_CV_AV_TPC_2D"); if (boolhist == false) gSystem->Exit(0);

			hCV2d = (TH2D*) hCV2dnue->Clone("hCV2d");
			hCV2d->Add(hCV2dnuebar); // Combine the fluxes

			double xbin_th = hCV2d->GetXaxis()->FindBin( 0.75*0.2065); // find the x bin to integrate from (threshold)
		
			flux  = hCV2d->Integral(xbin_th, hCV2d->GetNbinsX()+1, 0, hCV2d->GetNbinsY()+1); // Integrate over whole phase space (not quite any more)

			flux*= (POTScale / (GetPOT(fCV)*1.0e4)); // Scale to cm2 and the DATA POT


		}
		else {

			boolfile  = GetFile(fBeamline , Form( "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/output_uboone_run%i.root", index)); if (boolfile == false) gSystem->Exit(0);

			// Get the CV histogram in 2D
			boolhist = GetHist(fBeamline, hBeamline2dnue, "nue/Detsmear/nue_CV_AV_TPC_2D"); if (boolhist == false) gSystem->Exit(0);
			boolhist = GetHist(fBeamline, hBeamline2dnuebar, "nuebar/Detsmear/nuebar_CV_AV_TPC_2D"); if (boolhist == false) gSystem->Exit(0);
			
			hBeamline2d = (TH2D*) hBeamline2dnue->Clone("hBeamline2d");
			hBeamline2d->Add(hBeamline2dnuebar); // Combine the fluxes

			double xbin_th = hBeamline2d->GetXaxis()->FindBin( 0.75*0.2065); // find the x bin to integrate from (threshold)

			flux  = hBeamline2d->Integral(xbin_th, hBeamline2d->GetNbinsX()+1, 0, hBeamline2d->GetNbinsY()+1); // Integrate over whole phase space
			
			flux*= (POTScale / (GetPOT(fBeamline)*1.0e4)); // Scale to cm2 and the DATA POT
		}

	}

	return flux;

}
// ------------------------------------------------------------------------------------------------------------
// Function to calculate the data cross section
double CalcDataXSec(double sel, double bkg , double flux,
					double targets, double intime_cosmics_bkg, double intime_cosmic_scale_factor,
					double dirt, double dirt_scale_factor, double mc_scale_factor, double efficiency ){

	bool DEBUG{false};

	if (DEBUG) std::cout << 
	"DEBUG:\n"<<
	"sel:\t" << sel << "\n" << 
	"bkg:\t" << bkg  << "\n" << 
	"flux:\t" << flux << "\n" << 
	"targets:\t" << targets << "\n" << 
	"intime_cosmics_bkg:\t" << intime_cosmics_bkg << "\n" << 
	"intime cosmic scale factor:\t" << intime_cosmic_scale_factor << "\n" << 
	"dirt:\t" << dirt << "\n" << 
	"dirt scale factor:\t" << dirt_scale_factor << "\n" << 
	"mc scale factor:\t" << mc_scale_factor << "\n" << 
	"efficiency:\t" << efficiency << std::endl;

	if (DEBUG) std::cout << "Total Scaled background:\t" <<  (intime_cosmics_bkg * intime_cosmic_scale_factor) - (dirt * dirt_scale_factor) - (bkg * mc_scale_factor) << std::endl;	

	return (sel - (intime_cosmics_bkg * intime_cosmic_scale_factor) - (dirt * dirt_scale_factor) - (bkg * mc_scale_factor)) / (efficiency * targets * flux); 
}
// ------------------------------------------------------------------------------------------------------------
// Function to read in the event lists
void ReadEvents(const char *filename, std::vector<int> &N_evtnum ){
	
	std::ifstream fileIN; 

	fileIN.open(filename); // Open the file
	
	if (!fileIN.good()) {  // Check if the file opened correctly
		std::cerr << "Error: file:\t" << filename <<"\tcould not be opened" << std::endl;
		exit(1);
	}

	double temp{0}; // Use a temp var to get the values and push back

	if (fileIN.is_open()) { 
		
		while ( !fileIN.eof()) {   // loop over lines in file
			
			fileIN >> temp;        // Add number to temp var
			N_evtnum.push_back(temp);
		}
		
		fileIN.close();
	}
}
// ------------------------------------------------------------------------------------------------------------
// Function to convert detector co-ords to a beam coordinate theta
double GetTheta(double detx, double dety, double detz){

	// Variables
	TRotation RotDet2Beam;     		// Rotations
	TVector3  detxyz, BeamCoords; 	// Translations
	std::vector<double> rotmatrix;  // Inputs

	// input detector coordinates to translate
	detxyz = {detx, dety, detz}; 	

	// From detector to beam coords
	rotmatrix = {
		0.921038538,	4.625400126e-05,	-0.3894714486,
   		0.0227135048,	0.9982916247,		0.05383241394,
   		0.3888085752,	-0.05842798945, 	0.9194640079};

	// Return the TRotation
	TVector3 newX, newY, newZ;
	newX = TVector3(rotmatrix[0], rotmatrix[1], rotmatrix[2]);
	newY = TVector3(rotmatrix[3], rotmatrix[4], rotmatrix[5]);
	newZ = TVector3(rotmatrix[6], rotmatrix[7], rotmatrix[8]);

	RotDet2Beam.RotateAxes(newX, newY, newZ); // Return the TRotation now beam to det
	RotDet2Beam.Invert(); // Invert back to the det to beam rot matrix
	
	// Rotate to beam coords
	BeamCoords = RotDet2Beam * detxyz;

	TVector3 beam_dir = {0 , 0 , 1};
	double theta = BeamCoords.Angle(beam_dir) * 180 / 3.1415926;

	// std::cout << theta << std::endl;

	return theta;
}
// ------------------------------------------------------------------------------------------------------------
// Updated function for flux lists
void ReadEventList(const char *filename, std::vector<int> &N_evtnum,
					std::vector<std::string> &class_type, std::vector<int> &mc_nu_id,
					std::vector<double> &mc_nu_dir_x, std::vector<double> &mc_nu_dir_y, std::vector<double> &mc_nu_dir_z, std::vector<double> &mc_nu_energy    ){
	
	// event number, classifier type, mc_nu_id, mc_nu_dir_x, mc_nu_dir_y, mc_nu_dir_z, mc_nu_energy

	std::ifstream fileIN; 

	fileIN.open(filename); // Open the file
	
	if (!fileIN.good()) {  // Check if the file opened correctly
		std::cerr << "Error: file:\t" << filename <<"\tcould not be opened" << std::endl;
		exit(1);
	}

	int temp_evtnum,  temp_mc_nu_id;
	double temp_mc_nu_dir_x, temp_mc_nu_dir_y, temp_mc_nu_dir_z, temp_mc_nu_energy;
	std::string temp_class_type;

	if (fileIN.is_open()) { 
		
		// loop over lines in file
		while ( fileIN >> temp_evtnum >> temp_class_type >> temp_mc_nu_id >> temp_mc_nu_dir_x >> temp_mc_nu_dir_y >> temp_mc_nu_dir_z >> temp_mc_nu_energy) {  
			
			N_evtnum.push_back(temp_evtnum);
			class_type.push_back(temp_class_type);
			mc_nu_id.push_back(temp_mc_nu_id);
			mc_nu_dir_x.push_back(temp_mc_nu_dir_x);
			mc_nu_dir_y.push_back(temp_mc_nu_dir_y);
			mc_nu_dir_z.push_back(temp_mc_nu_dir_z);
			mc_nu_energy.push_back(temp_mc_nu_energy);
		}
		
		fileIN.close();
	}
}
// ------------------------------------------------------------------------------------------------------------
//  Function to calculate histogram ratios and return them as a vector of Th2D
void PrecalcHistRatio(HistWeights &flav, const char* mode, std::vector<std::string> params){
	TH2D *hBeamline2d , *hHP2d, *hCV2d, *hCV2d_Beam;
	TFile *fBeamline, *fCV, *fCV_Beam;
	bool boolhist, boolfile;
	int b=1, h=0;

	std::cout << "Calculating Hist ratios for flavour:\t" << mode << std::endl;

	// File with CV
	boolfile  = GetFile(fCV , "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/output_uboone_run0.root"); if (boolfile == false) gSystem->Exit(0); // Most up to date version of CV
	boolhist = GetHist(fCV, hCV2d, Form("%s/Detsmear/%s_CV_AV_TPC_2D", mode, mode)); if (boolhist == false) gSystem->Exit(0); // Get the CV

	// Loop over the parameters
	for (unsigned int i=1; i < params.size(); i++){
		// std::cout << params.at(i) << std::endl;

		// --- ---- ----- HP ---- --- --- ----- //
		if (params.at(i).find("PPFX") != std::string::npos) {
			// std::cout << params.at(i) << std::endl;
			
			for (unsigned int k=0; k<100; k++){
				boolhist = GetHist(fCV, hHP2d, Form("%s/Multisims/%s_%s_Uni_%i_AV_TPC_2D",mode, mode , params.at(i).c_str() , k)); if (boolhist == false) gSystem->Exit(0);
				hHP2d->Divide(hCV2d); // Divide hists
				flav.HP.at(h).push_back(hHP2d); // push back to vector
			}
			h++;
		}
		// --- ---- ----- Beamline ---- --- --- ----- //
		else {
			// std::cout << params.at(i) << std::endl;
			boolfile  = GetFile(fBeamline , Form("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/output_uboone_run%i.root", b)); if (boolfile == false) gSystem->Exit(0);
			boolhist  = GetHist(fBeamline, hBeamline2d, Form("%s/Detsmear/%s_CV_AV_TPC_2D",mode, mode)); if (boolhist == false) gSystem->Exit(0);
			hBeamline2d->Divide(hCV2d); // Divide hists
			flav.Beamline.push_back(hBeamline2d); // push back to vector
			b++;
		}
	}

}
// ------------------------------------------------------------------------------------------------------------
// Function to caluclate the standard deviation 
double STD_Calculator(std::vector<double> vec, double CV){
	double Err{0};
	
	for (unsigned int i = 0; i < vec.size(); i++ ){

        Err +=  (vec[i] - CV) * (vec[i]- CV);  
    }

    return (std::sqrt( Err / vec.size() ) );

}
// ------------------------------------------------------------------------------------------------------------
// Make a plot of the beamline variation uncertainties
void Make_Beamline_Plot(){

	std::vector<std::string> param_max = {
		"Horn 2kA",
		"Horn 1 x 3mm",
		"Horn 1 y 3mm",
		"Beam Spot 2mm",
		"Horn 2 x 3mm",
		"Horn 2 y 3mm",
		"Horns water",
		"Beam Shift x 1mm",
		"Beam Shift y 1mm",
		"Target z 7mm",
		"Horn1 Refined Descr.",
		"Decay Pipe B Field",
		"Old Horn",
		" ",
		"Total Error"
	};

	std::vector<double> param_max_val = {
		0.3,
		1.1,
		0.98,
		2.8,
		0.4,
		0.53,
		0.86,
		3.2,
		3.2,
		1.3,
		0.65,
		1.4,
		2.3,
		0
	};

	double tot_beamline_err{0};

	for (unsigned int i = 0; i < param_max_val.size(); i++){

		tot_beamline_err+=param_max_val[i]*param_max_val[i];
	
	}
	tot_beamline_err = std::sqrt(tot_beamline_err);
	param_max_val.push_back(tot_beamline_err);

	TH1D *hBeamline = new TH1D("Beamline","", param_max_val.size()+1, 0, param_max_val.size()+1);

	for (unsigned int i = 0; i < param_max_val.size(); i++){
		hBeamline->Fill(param_max[i].c_str(), param_max_val[i]);

	}
	gStyle->SetOptStat(0); // say no to stats box
	TCanvas* c = new TCanvas();
	hBeamline->SetLineColor(kViolet-6);
	hBeamline->SetLineWidth(3);
	hBeamline->GetYaxis()->SetTitle("Percentage Uncertainty %");
	hBeamline->LabelsOption("v");
	gPad->SetBottomMargin(0.33);

	hBeamline->GetXaxis()->SetLabelSize(0.05);
	hBeamline->GetXaxis()->SetTitleSize(0.05);
	hBeamline->GetYaxis()->SetLabelSize(0.05);
	hBeamline->GetYaxis()->SetTitleSize(0.05);
	hBeamline->SetMarkerSize(1.8);
	gPad->SetLeftMargin(0.15);

	hBeamline->Draw("hist, text00");

	c->Print("plots/Beamline_Uncertainties.pdf");

	return;
}
// ------------------------------------------------------------------------------------------------------------
// Make a plot of the hadron production variation uncertainties
void Make_HP_Plot(std::vector<double> HP_uncertainties){

	std::vector<std::string> HP_names = {
		"Other",
		"Targ. Atten",
		"Thin Kaon",
		"Thin Meson",
		"Thin Neutron",
		"Thin Nuc. A",
		"Thin Nuc",
		"Thin Pion",
		"Tot. Absorp"
	};

	double HP_Quad{0};
	// Get the Quadrature Sum, start from 1 to skip master
	for (unsigned int i = 1; i < HP_uncertainties.size(); i++){

		HP_Quad+=  HP_uncertainties.at(i) * HP_uncertainties.at(i);
	}

	HP_Quad = std::sqrt(HP_Quad);

	TH1D *hHP = new TH1D("HP","", 14, 0, 14);

	for (unsigned int i = 1; i < HP_uncertainties.size(); i++){
		hHP->Fill( HP_names.at(i-1).c_str(), HP_uncertainties.at(i) );

	}
	hHP->Fill( " ", 0 );
	hHP->Fill( "Master", HP_uncertainties.at(0) );
	hHP->Fill( "  ", 0 );
	hHP->Fill( "Quadrature Sum", HP_Quad );
	hHP->Fill( "  ", 0 );
	
	gStyle->SetPaintTextFormat("4.2f");
	gStyle->SetOptStat(0); // say no to stats box
	TCanvas* c_HP = new TCanvas();
	hHP->SetLineColor(kViolet-6);
	hHP->SetLineWidth(3);
	hHP->GetYaxis()->SetTitle("Percentage Uncertainty %");
	hHP->LabelsOption("v");
	gPad->SetBottomMargin(0.25);

	hHP->GetXaxis()->SetLabelSize(0.05);
	hHP->GetXaxis()->SetTitleSize(0.05);
	hHP->GetYaxis()->SetLabelSize(0.05);
	hHP->GetYaxis()->SetTitleSize(0.05);
	hHP->SetMarkerSize(1.8);
	gPad->SetLeftMargin(0.15);

	hHP->Draw("hist, text00");

	c_HP->Print("plots/HP_Uncertainties.pdf");



}
