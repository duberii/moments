#include <TBranch.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TH1D.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TLorentzRotation.h>
#include "FSAnalysis/FSHistogramAnalysis.h"
#include "FSBasic/FSHistogram.h"
#include "FSBasic/FSString.h"
#include "FSBasic/FSTree.h"
#include "FSBasic/FSMath.h"
#include "FSMode/FSModeInfo.h"


double d_lm0(int L, int M, double costheta)
{
    return pow(-1., M) * sqrt(TMath::Factorial(L - M) / TMath::Factorial(L + M)) * std::assoc_legendre(L, M, costheta);
}

double get_polarization(double labPolAngle, double photon_energy)
{
    static TH1D* pol_0_hist = nullptr;
    static TH1D* pol_45_hist = nullptr;
    static TH1D* pol_90_hist = nullptr;
    static TH1D* pol_135_hist = nullptr;

    // Thread-safe initialization - guaranteed to run exactly once
    if (!pol_0_hist) {
        TFile* pols = TFile::Open("/N/u/rdube/Quartz/data/makePolVals2020V2/outFiles/sp20TPol.root");
        if (!pols || pols->IsZombie()) {
            std::cerr << "Error: Failed to open polarization file!" << std::endl;
            return -9999999.;
        }

        // Get and clone histograms (they will persist after file closes)
        TH1D* h0 = dynamic_cast<TH1D*>(pols->Get("hPol0"));
        TH1D* h45 = dynamic_cast<TH1D*>(pols->Get("hPol45"));
        TH1D* h90 = dynamic_cast<TH1D*>(pols->Get("hPol90"));
        TH1D* h135 = dynamic_cast<TH1D*>(pols->Get("hPol135"));

        if (!h0 || !h45 || !h90 || !h135) {
            std::cerr << "Error: Failed to find one or more polarization histograms!" << std::endl;
            pols->Close();
            delete pols;
            return -9999999.;
        }

        // Clone histograms so they survive file closure
        pol_0_hist = static_cast<TH1D*>(h0->Clone("pol_0_hist_clone"));
        pol_45_hist = static_cast<TH1D*>(h45->Clone("pol_45_hist_clone"));
        pol_90_hist = static_cast<TH1D*>(h90->Clone("pol_90_hist_clone"));
        pol_135_hist = static_cast<TH1D*>(h135->Clone("pol_135_hist_clone"));

        // Set directory to nullptr so ROOT doesn't try to manage them
        pol_0_hist->SetDirectory(nullptr);
        pol_45_hist->SetDirectory(nullptr);
        pol_90_hist->SetDirectory(nullptr);
        pol_135_hist->SetDirectory(nullptr);

        pols->Close();
        delete pols;
    }

    // Check if initialization failed
    if (!pol_0_hist || !pol_45_hist || !pol_90_hist || !pol_135_hist) {
        return 999999999.; // Error value
    }

    // Select appropriate histogram
    TH1D* selected_hist = nullptr;
    if (labPolAngle == 0.) {
        selected_hist = pol_0_hist;
    } else if (labPolAngle == 45.) {
        selected_hist = pol_45_hist;
    } else if (labPolAngle == 90.) {
        selected_hist = pol_90_hist;
    } else if (labPolAngle == 135.) {
        selected_hist = pol_135_hist;
    } else {
        return 999999999.; // Invalid angle
    }

    // Get bin content (thread-safe read operation)
    int bin = selected_hist->FindBin(photon_energy);
    return selected_hist->GetBinContent(bin);
}

double polphi(const double& PxPA, const double& PyPA, const double& PzPA, const double& EnPA,
    const double& PxPB, const double& PyPB, const double& PzPB, const double& EnPB,
    const double& PxPC, const double& PyPC, const double& PzPC, const double& EnPC,
    const double& PxPD, const double& PyPD, const double& PzPD, const double& EnPD, const double& labPolAngle)
{
    // https://github.com/JeffersonLab/halld_sim/blob/e84dafa3acdcdc126a300b5fc1f02ba89e29a490/src/libraries/AMPTOOLS_AMPS/TwoPiAngles.cc
    TLorentzVector beam(PxPA, PyPA, PzPA, EnPA);
    TLorentzVector recoil(PxPB, PyPB, PzPB, EnPB);
    TLorentzVector p1(PxPC, PyPC, PzPC, EnPC);
    TLorentzVector p2(PxPD, PyPD, PzPD, EnPD);

    TLorentzVector resonance = p1 + p2;
    TLorentzRotation resonanceBoost(-resonance.BoostVector());
    TLorentzVector recoil_res = resonanceBoost * recoil;

    // normal to the production plane
    TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();

    TVector3 eps(1.0, 0.0, 0.0); // reference beam polarization vector at 0 degrees
    double correctPolAngle = atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y))) + labPolAngle * 0.017453293;
    return correctPolAngle;
}

void calculate_columns(
    const char* input_filename,
    const char* output_filename,
    const char* final_state,
    bool isGenMC)
{
    int maxL=8;
    std::cout << "========================================" << std::endl;
    std::cout << "Calculating angles from FSRoot tree" << std::endl;
    std::cout << "Input:  " << input_filename << std::endl;
    std::cout << "Output: " << output_filename << std::endl;
    std::cout << "Final State:   " << final_state << std::endl;
    std::cout << "is gen MC?      " << isGenMC << std::endl;
    std::cout << "========================================" << std::endl;
    
    bool ks_sideband = false;
    bool use_chi2diff = false;
    std::string tree_name_str;
    if(std::string(final_state) == "kskl") {
        tree_name_str = "ntFSGlueX_100_1000";
        ks_sideband = true;
    } else if(std::string(final_state) == "kpkm") {
        tree_name_str = "ntFSGlueX_100_110000";
        use_chi2diff=true;
    }
    const char* tree_name = tree_name_str.c_str();

    // Open input file
    TFile* input_file = TFile::Open(input_filename, "READ");
    if (!input_file || input_file->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << input_filename << std::endl;
        return;
    }

    // Get input tree
    TTree* input_tree = (TTree*)input_file->Get(tree_name);
    if (!input_tree) {
        std::cerr << "Error: Cannot find tree " << tree_name << " in input file" << std::endl;
        input_file->Close();
        return;
    }
    TFile* input_file_friend;
    TTree* friend_tree;
    //get friend tree
    if ( !isGenMC ) {
        input_file_friend= TFile::Open(TString::Format("%s.Chi2DOFRank",input_filename).Data(), "READ");
        if (!input_file_friend || input_file_friend->IsZombie()) {
            std::cerr << "Error: Cannot open input file " << input_filename << ".Chi2DOFRank" <<std::endl;
            return;
        }
        friend_tree = (TTree*)input_file_friend->Get(TString::Format("%s_Chi2DOFRank",tree_name).Data());
        if (!friend_tree) {
            std::cerr << "Error: Cannot find tree " << tree_name << "_Chi2DOFRank in input file" << std::endl;
            input_file->Close();
            input_file_friend->Close();
            return;
        }
    }
    // Create output file
    TFile* output_file = TFile::Open(output_filename, "RECREATE");
    if (!output_file || output_file->IsZombie()) {
        std::cerr << "Error: Cannot create output file " << output_filename << std::endl;
        input_file->Close();
        return;
    }

    // Create output tree
    TTree* output_tree = new TTree("decayAngles", "calculated decay angles for KKbar final states");

    // ============================================
    // Set up variables for input branches
    // ============================================

    std::cout<< "Setting up variables for input branches" << std::endl;

    //For both genMC and accMC/data:

    double EnP1 = -999.0, EnP2 = -999.0, EnP3 = -999.0, EnPB = -999.0;
    double PxP1 = -999.0, PxP2 = -999.0, PxP3 = -999.0, PxPB = -999.0;
    double PyP1 = -999.0, PyP2 = -999.0, PyP3 = -999.0, PyPB = -999.0;
    double PzP1 = -999.0, PzP2 = -999.0, PzP3 = -999.0, PzPB = -999.0;
    double Weight = 1.0;

    //Only for accMC/data:

    double REnP1 = -999.0, REnP2 = -999.0, REnP3 = -999.0, REnPB = -999.0;
    double RPxP1 = -999.0, RPxP2 = -999.0, RPxP3 = -999.0, RPxPB = -999.0;
    double RPyP1 = -999.0, RPyP2 = -999.0, RPyP3 = -999.0, RPyPB = -999.0;
    double RPzP1 = -999.0, RPzP2 = -999.0, RPzP3 = -999.0, RPzPB = -999.0;
    double Chi2NDF = -999.0;
    int Chi2NDFpipi = -999;
    double RFDeltaT = -999.0;
    double AccidentalScale = -999.0;
    double KSFlightSignificance = -999.0;
    double ProtonVertexZ = -999.0;

    // ============================================
    // Set up variables for output branches
    // ============================================

    std::cout<< "Setting up variables for output branches" << std::endl;

    // For both genMC and accMC/data:

    double PolarizationDegree = -999.0;
    double MandelstamNegT = 999.0;
    double MesonResonanceMass = -999.0;
    double BaryonResonanceMass2 = -999.0;
    double BaryonResonanceMass3 = -999.0;

    double PolarizationAngleReac = -999.0;
    double helCosTheta = -999.0;
    double helPhi = -999.0;
    double gjCosTheta = -999.0;
    double gjPhi = -999.0;

    double PolarizationAngleLab = -999.0;
    double labCosTheta = -999.0;
    double labPhi=-999.0;

    std::vector<double> d;
    d.resize((maxL+1)*(maxL+2)/2);
    std::vector<double> cosMphis;
    cosMphis.resize(maxL+1);

    //Only for accMC/data:

    double MissingMass = -999.0;
    double Chi2pipiMinusChi2KK = -999.0;
    double NumUnusedShowers = -999.0;
    double NumUnusedTracks = -999.0;
    
    int Chi2DOFRankGlobal = -999;
    int Chi2DOFRank = -999;

    // ============================================
    // Setting Input Branch Addresses
    // ============================================

    std::cout<< "Setting input branch addresses" << std::endl;

    input_tree->SetBranchAddress("PolarizationAngle", &PolarizationAngleLab);
    if (!isGenMC) {
        input_tree->SetBranchAddress("RVzP1",&ProtonVertexZ);
        input_tree->SetBranchAddress("EnPB", &EnPB);
        input_tree->SetBranchAddress("PxPB", &PxPB);
        input_tree->SetBranchAddress("PyPB", &PyPB);
        input_tree->SetBranchAddress("PzPB", &PzPB);
        input_tree->SetBranchAddress("EnP1", &EnP1);
        input_tree->SetBranchAddress("PxP1", &PxP1);
        input_tree->SetBranchAddress("PyP1", &PyP1);
        input_tree->SetBranchAddress("PzP1", &PzP1);
        input_tree->SetBranchAddress("EnP2", &EnP2);
        input_tree->SetBranchAddress("PxP2", &PxP2);
        input_tree->SetBranchAddress("PyP2", &PyP2);
        input_tree->SetBranchAddress("PzP2", &PzP2);
        input_tree->SetBranchAddress("REnPB", &REnPB);
        input_tree->SetBranchAddress("RPxPB", &RPxPB);
        input_tree->SetBranchAddress("RPyPB", &RPyPB);
        input_tree->SetBranchAddress("RPzPB", &RPzPB);
        input_tree->SetBranchAddress("REnP1", &REnP1);
        input_tree->SetBranchAddress("RPxP1", &RPxP1);
        input_tree->SetBranchAddress("RPyP1", &RPyP1);
        input_tree->SetBranchAddress("RPzP1", &RPzP1);
        input_tree->SetBranchAddress("REnP2", &REnP2);
        input_tree->SetBranchAddress("RPxP2", &RPxP2);
        input_tree->SetBranchAddress("RPyP2", &RPyP2);
        input_tree->SetBranchAddress("RPzP2", &RPzP2);
        input_tree->SetBranchAddress("AccidentalScale",&AccidentalScale);
        input_tree->SetBranchAddress("RFDeltaT",&RFDeltaT);
        input_tree->SetBranchAddress("NumUnusedTracks",&NumUnusedTracks);
        input_tree->SetBranchAddress("NumNeutralHypos",&NumUnusedShowers);
        input_tree->SetBranchAddress("Chi2DOF",&Chi2NDF);
        friend_tree->SetBranchAddress("Chi2DOFRank",&Chi2DOFRank);
        friend_tree->SetBranchAddress("Chi2DOFRankGlobal",&Chi2DOFRankGlobal);
        if (use_chi2diff) {
            friend_tree->SetBranchAddress("Chi2DOFRankVarBestOther", &Chi2NDFpipi);
        }
        if (ks_sideband) {
            input_tree->SetBranchAddress("VeeLSigmaP2",&KSFlightSignificance);
        }
        if (std::string(final_state) == "kpkm") {
            input_tree->SetBranchAddress("EnP3", &EnP3);
            input_tree->SetBranchAddress("PxP3", &PxP3);
            input_tree->SetBranchAddress("PyP3", &PyP3);
            input_tree->SetBranchAddress("PzP3", &PzP3);
            input_tree->SetBranchAddress("REnP3", &REnP3);
            input_tree->SetBranchAddress("RPxP3", &RPxP3);
            input_tree->SetBranchAddress("RPyP3", &RPyP3);
            input_tree->SetBranchAddress("RPzP3", &RPzP3);
        }
        input_tree->SetBranchAddress("AccidentalScale", &AccidentalScale);
    } else {
        input_tree->SetBranchAddress("MCEnPB", &EnPB);
        input_tree->SetBranchAddress("MCPxPB", &PxPB);
        input_tree->SetBranchAddress("MCPyPB", &PyPB);
        input_tree->SetBranchAddress("MCPzPB", &PzPB);
        input_tree->SetBranchAddress("MCEnP1", &EnP1);
        input_tree->SetBranchAddress("MCPxP1", &PxP1);
        input_tree->SetBranchAddress("MCPyP1", &PyP1);
        input_tree->SetBranchAddress("MCPzP1", &PzP1);
        input_tree->SetBranchAddress("MCEnP2", &EnP2);
        input_tree->SetBranchAddress("MCPxP2", &PxP2);
        input_tree->SetBranchAddress("MCPyP2", &PyP2);
        input_tree->SetBranchAddress("MCPzP2", &PzP2);
        if (std::string(final_state) == "kpkm") {
            input_tree->SetBranchAddress("MCEnP3", &EnP3);
            input_tree->SetBranchAddress("MCPxP3", &PxP3);
            input_tree->SetBranchAddress("MCPyP3", &PyP3);
            input_tree->SetBranchAddress("MCPzP3", &PzP3);
        }
    }

    // ============================================
    // Set up output branches
    // ============================================

    std::cout<< "Setting up output branches" << std::endl;

    // For both genMC and accMC/data:

    output_tree->Branch("PolarizationAngleReac", &PolarizationAngleReac, "PolarizationAngleReac/D");
    output_tree->Branch("MandelstamNegT", &MandelstamNegT, "MandelstamNegT/D");
    output_tree->Branch("PolarizationDegree", &PolarizationDegree, "PolarizationDegree/D");

    output_tree->Branch("helCosTheta", &helCosTheta, "helCosTheta/D");
    output_tree->Branch("helPhi", &helPhi, "helPhi/D");
    output_tree->Branch("gjCosTheta", &gjCosTheta, "gjCosTheta/D");
    output_tree->Branch("gjPhi", &gjPhi, "gjPhi/D");

    output_tree->Branch("BaryonResonanceMass2", &BaryonResonanceMass2, "BaryonResonanceMass2/D");
    output_tree->Branch("BaryonResonanceMass3", &BaryonResonanceMass3, "BaryonResonanceMass3/D");
    output_tree->Branch("MesonResonanceMass", &MesonResonanceMass, "MesonResonanceMass/D");

    output_tree->Branch("BeamEnergy", &EnPB, "BeamEnergy/D");
    
    output_tree->Branch("PolarizationAngleLab", &PolarizationAngleLab, "PolarizationAngleLab/D");
    output_tree->Branch("labCosTheta",&labCosTheta, "labCosTheta/D");
    output_tree->Branch("labPhi",&labPhi, "labPhi/D");
    for (int L = 0; L<= maxL; L++) {
        output_tree->Branch(TString::Format("cos_%iphi", L),&cosMphis[L],TString::Format("cos_%iphi/D", L));
        for (int M = 0; M <= L; M++) {
            int index = M+L*(L+1)/2;
            output_tree->Branch(TString::Format("d_%i_%i_0", L,M),&d[index], TString::Format("d_%i_%i_0/D", L,M));
        }
    }
    output_tree->Branch("Weight", &Weight,"Weight/D");

    // Only for accMC:
    
    if (!isGenMC) {
        output_tree->Branch("ProtonVertexZ",&ProtonVertexZ,"ProtonVertexZ/D");
        output_tree->Branch("NumUnusedShowers",&NumUnusedShowers,"NumUnusedShowers/D");
        output_tree->Branch("NumUnusedTracks",&NumUnusedTracks,"NumUnusedTracks/D");
        output_tree->Branch("Chi2NDF",&Chi2NDF, "Chi2NDF/D");
        if (ks_sideband) {
            output_tree->Branch("KSFlightSignificance",&KSFlightSignificance,"KSFlightSignificance/D");
            output_tree->Branch("MissingMass",&MissingMass,"MissingMass/D");
        }
        if (use_chi2diff) {
            output_tree->Branch("Chi2pipiMinusChi2KK", &Chi2pipiMinusChi2KK, "Chi2pipiMinusChi2KK/D");
        }
    }
    
    

    // ============================================
    // Loop over events and fill output tree
    // ============================================

    std::cout<< "Beginning event loop" << std::endl;

    Long64_t nEntries = input_tree->GetEntries();
    if (!isGenMC) {
        Long64_t nEntries_friend = friend_tree->GetEntries();
        if (nEntries != nEntries_friend) {
            std::cout << "Mismatch between friend entries and tree entries. Main tree has " << nEntries << " entries. Friend Tree has " << nEntries_friend << " entries." << std::endl;
            return;
        }
    }
    std::cout << "\nProcessing " << nEntries << " entries..." << std::endl;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        // Load entry
        input_tree->GetEntry(entry);

        // Skip events early on (if accMC or data) if their chi2 rank is not 1.
        
        if (!isGenMC) {
            friend_tree->GetEntry(entry);
            if (Chi2DOFRank != 1) {
                continue;
            }
            if (use_chi2diff) {
                Chi2pipiMinusChi2KK = (double(Chi2NDFpipi)/1000.)-Chi2NDF;
                if (Chi2DOFRankGlobal != 1) {
                    continue;
                }
            }
        }

        // Set up KL 4-vector if final state is kskl.

        if (std::string(final_state) == "kskl") {
            PxP3 = PxPB - PxP1 - PxP2;
            PyP3 = PyPB - PyP1 - PyP2;
            PzP3 = PzPB - PzP1 - PzP2;
            EnP3 = 0.938 + EnPB - EnP1 - EnP2;
        }

        MesonResonanceMass = sqrt((EnP2+EnP3)*(EnP2+EnP3)-(PxP2+PxP3)*(PxP2+PxP3)-(PyP2+PyP3)*(PyP2+PyP3)-(PzP2+PzP3)*(PzP2+PzP3));

	    if (std::string(final_state) == "kskl") {
            MandelstamNegT=2*0.938*(0.938-EnP1)-(pow(MesonResonanceMass,4)/(4*0.938)-pow(sqrt(pow(FSMath::boostEnergy(PxP1,PyP1,PxP1,EnP1,PxPB,PyPB,PzPB,EnPB+0.938),2)-0.938*0.938)-sqrt(pow(FSMath::boostEnergy(0,0,0,0.938,PxPB,PyPB,PzPB,EnPB+0.938),2)-0.938*0.938),2)); 
            if (!isGenMC) {
                MissingMass = sqrt(pow(((REnPB+0.938272)-(REnP1+REnP2)),2)-pow(((RPxPB+0.0)-(RPxP1+RPxP2)),2)-pow(((RPyPB+0.0)-(RPyP1+RPyP2)),2)-pow(((RPzPB+0.0)-(RPzP1+RPzP2)),2));
            }
        } else {
            MandelstamNegT = (PxPB-PxP2-PxP3)*(PxPB-PxP2-PxP3)+(PyPB-PyP2-PyP3)*(PyPB-PyP2-PyP3)+(PzPB-PzP2-PzP3)*(PzPB-PzP2-PzP3)-(EnPB-EnP2-EnP3)*(EnPB-EnP2-EnP3);
        }
        BaryonResonanceMass2 = sqrt((EnP1+EnP2)*(EnP1+EnP2)-(PxP1+PxP2)*(PxP1+PxP2)-(PyP1+PyP2)*(PyP1+PyP2)-(PzP2+PzP1)*(PzP2+PzP1));
        BaryonResonanceMass3 = sqrt((EnP1+EnP3)*(EnP1+EnP3)-(PxP1+PxP3)*(PxP1+PxP3)-(PyP1+PyP3)*(PyP1+PyP3)-(PzP3+PzP1)*(PzP3+PzP1));
        PolarizationDegree = get_polarization(PolarizationAngleLab, EnPB);
        
        //Getting the helicity angles
        TLorentzVector beam(PxPB, PyPB, PzPB, EnPB); 
        TLorentzVector recoil(PxP1, PyP1, PzP1, EnP1);
        TLorentzVector p1(PxP2, PyP2, PzP2, EnP2);
        TLorentzVector p2(PxP3, PyP3, PzP3, EnP3);
        TLorentzVector resonance = p1 + p2;
        TLorentzRotation resonanceBoost( -resonance.BoostVector() );
        TLorentzVector beam_res = resonanceBoost * beam;
        TLorentzVector recoil_res = resonanceBoost * recoil;
        TLorentzVector p1_res = resonanceBoost * p1;
        // normal to the production plane
        TVector3 y = (beam.Vect().Unit().Cross(-recoil.Vect().Unit())).Unit();
        // choose helicity frame: z-axis opposite recoil proton in rho rest frame
        TVector3 z = -1. * recoil_res.Vect().Unit();
        TVector3 x = y.Cross(z).Unit();
        TVector3 angles( (p1_res.Vect()).Dot(x),(p1_res.Vect()).Dot(y),(p1_res.Vect()).Dot(z) );
        TVector3 eps(1.0, 0.0, 0.0); // reference beam polarization vector at 0 degrees
        PolarizationAngleReac = std::fmod(atan2(y.Dot(eps), beam.Vect().Unit().Dot(eps.Cross(y))) + PolarizationAngleLab * 0.017453293+6.2832,6.2832);
        helCosTheta = angles.CosTheta();
        helPhi = angles.Phi();
        gjCosTheta = FSMath::gjcostheta(PxP2, PyP2, PzP2, EnP2,PxP3, PyP3, PzP3, EnP3,PxPB, PyPB, PzPB, EnPB);
        gjPhi = FSMath::gjphi(PxP2, PyP2, PzP2, EnP2, PxP3, PyP3, PzP3, EnP3,PxP1, PyP1, PzP1, EnP1,PxPB, PyPB, PzPB, EnPB);
        labPhi = atan2(PyP2, PxP2);
        labCosTheta = PzP2/sqrt(PxP2*PxP2+PyP2*PyP2+PzP2*PzP2);
        for (int L = 0; L<= maxL; L++) {
            cosMphis[L] = cos(L*helPhi);
            for (int M = 0; M <= L; M++) {
                int index = M+L*(L+1)/2;
                d[index] = d_lm0(L,M,helCosTheta);
            }
        }
        if (!isGenMC) {
            Weight = 1.0;
            if (abs(RFDeltaT) > 6 && abs(RFDeltaT) < 14) {
                Weight *= -0.25;
            } else if (abs(RFDeltaT) > 2 && abs(RFDeltaT) < 6) {
                continue;
            }
            if (ks_sideband) {
                double KS_mass = sqrt(EnP2*EnP2 - PxP2*PxP2-PyP2*PyP2 - PzP2*PzP2);
                if (abs(KS_mass-0.5)>0.04 && abs(KS_mass-0.5) < 0.06) {
                    Weight *= -1.0;
                } else if (abs(KS_mass-0.5) > 0.02 && abs(KS_mass-0.5)<0.04) {
                    continue;
                } else if (abs(KS_mass-0.5)>0.06) {
                    continue;
                }
            }
        }
        // Fill output tree
        output_tree->Fill();
    }

    // ============================================
    // Save and cleanup
    // ============================================

    std::cout << "\nWriting output file..." << std::endl;
    output_file->cd();
    output_tree->Write();

    std::cout << "\nOutput tree statistics:" << std::endl;
    std::cout << "  Entries: " << output_tree->GetEntries() << std::endl;
    std::cout << "  Branches: " << output_tree->GetListOfBranches()->GetEntries() << std::endl;

    output_file->Close();
    input_file->Close();

    std::cout << "\nConversion complete!" << std::endl;
    std::cout << "Output file: " << output_filename << std::endl;
    std::cout << "Output tree: decayAngles" << std::endl;
}
