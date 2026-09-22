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
        TFile* pols = TFile::Open("/N/u/rdube/Quartz/software/makePolValsV9/outFiles/sp18TPol75.root");
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

void calculate_columns_kskl_GlueX_I(
    const char* input_filename,
    const char* output_filename,
    bool isBkg)
{
    int maxL=8;
    std::cout << "========================================" << std::endl;
    std::cout << "Calculating angles from FSRoot tree" << std::endl;
    std::cout << "Input:  " << input_filename << std::endl;
    std::cout << "Output: " << output_filename << std::endl;
    std::cout << "is bkg?      " << isBkg << std::endl;
    std::cout << "========================================" << std::endl;
    
    const char* tree_name = "kin";

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

    float Px_arr[3], Py_arr[3], Pz_arr[3], En_arr[3];
    float EnPB_f, PxPB_f, PyPB_f, PzPB_f;

    double EnP1 = -999.0, EnP2 = -999.0, EnP3 = -999.0, EnPB = -999.0;
    double PxP1 = -999.0, PxP2 = -999.0, PxP3 = -999.0, PxPB = -999.0;
    double PyP1 = -999.0, PyP2 = -999.0, PyP3 = -999.0, PyPB = -999.0;
    double PzP1 = -999.0, PzP2 = -999.0, PzP3 = -999.0, PzPB = -999.0;
    float Weight_f = -999.0;
    int pol_angle_i = -1.0;

    // ============================================
    // Set up variables for output branches
    // ============================================

    std::cout<< "Setting up variables for output branches" << std::endl;

    // For both genMC and accMC/data:
    double Weight = 1.0;
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

    double VanHoveX = -999.0;
    double VanHoveY = -999.0;
    double VanHoveOmega = -999.0;

    double PolarizationAngleLab = -999.0;
    double labCosTheta = -999.0;
    double labPhi=-999.0;
    double KSMass = -999.0;

    std::vector<double> d;
    d.resize((maxL+1)*(maxL+2)/2);
    std::vector<double> cosMphis;
    cosMphis.resize(maxL+1);

    // ============================================
    // Setting Input Branch Addresses
    // ============================================

    std::cout<< "Setting input branch addresses" << std::endl;

    input_tree->SetBranchAddress("pol_angle", &pol_angle_i);
    input_tree->SetBranchAddress("Px_FinalState",&Px_arr);
    input_tree->SetBranchAddress("Py_FinalState",&Py_arr);
    input_tree->SetBranchAddress("Pz_FinalState",&Pz_arr);
    input_tree->SetBranchAddress("E_FinalState",&En_arr);

    input_tree->SetBranchAddress("Px_FinalState",&Px_arr);
    input_tree->SetBranchAddress("Py_FinalState",&Py_arr);
    input_tree->SetBranchAddress("Pz_FinalState",&Pz_arr);
    input_tree->SetBranchAddress("E_Beam",&EnPB_f);
    input_tree->SetBranchAddress("Px_Beam",&PxPB_f);
    input_tree->SetBranchAddress("Py_Beam",&PyPB_f);
    input_tree->SetBranchAddress("Pz_Beam",&PzPB_f);

    input_tree->SetBranchAddress("Weight",&Weight_f);

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

    output_tree->Branch("VanHoveX",&VanHoveX, "VanHoveX/D");
    output_tree->Branch("VanHoveY",&VanHoveY, "VanHoveY/D");
    output_tree->Branch("VanHoveOmega",&VanHoveOmega, "VanHoveOmega/D");

    output_tree->Branch("BeamEnergy", &EnPB, "BeamEnergy/D");
    output_tree->Branch("KSMass", &KSMass, "KSMass/D");
    
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

    // ============================================
    // Loop over events and fill output tree
    // ============================================

    std::cout<< "Beginning event loop" << std::endl;

    Long64_t nEntries = input_tree->GetEntries();
    std::cout << "\nProcessing " << nEntries << " entries..." << std::endl;

    for (Long64_t entry = 0; entry < nEntries; entry++) {
        // Load entry
        input_tree->GetEntry(entry);

        // Set up KL 4-vector if final state is kskl.
        EnP1 = double((En_arr)[0]);
        EnP2 = double((En_arr)[1]);
        EnP3 = double((En_arr)[2]);
        PxP1 = double((Px_arr)[0]);
        PxP2 = double((Px_arr)[1]);
        PxP3 = double((Px_arr)[2]);
        PyP1 = double((Py_arr)[0]);
        PyP2 = double((Py_arr)[1]);
        PyP3 = double((Py_arr)[2]);
        PzP1 = double((Pz_arr)[0]);
        PzP2 = double((Pz_arr)[1]);
        PzP3 = double((Pz_arr)[2]);

        EnPB = EnPB_f;
        PxPB = PxPB_f;
        PyPB = PyPB_f;
        PzPB = PzPB_f;

        Weight = Weight_f;

        KSMass = sqrt(EnP2*EnP2 - PxP2*PxP2 - PyP2*PyP2 - PzP2*PzP2);

        PolarizationAngleLab = double(pol_angle_i);

        MesonResonanceMass = sqrt((EnP2+EnP3)*(EnP2+EnP3)-(PxP2+PxP3)*(PxP2+PxP3)-(PyP2+PyP3)*(PyP2+PyP3)-(PzP2+PzP3)*(PzP2+PzP3));
        MandelstamNegT=-1*(pow(((0.938272)-(EnP1)),2)-pow(((0.0)-(PxP1)),2)-pow(((0.0)-(PyP1)),2)-pow(((0.0)-(PzP1)),2))+pow((pow(((EnPB+0.938272)-(EnP1)),2)-pow(((PxPB+0.0)-(PxP1)),2)-pow(((PyPB+0.0)-(PyP1)),2)-pow(((PzPB+0.0)-(PzP1)),2))/(2*(sqrt(pow(((EnPB+0.938272)),2)-pow(((PxPB+0.0)),2)-pow(((PyPB+0.0)),2)-pow(((PzPB+0.0)),2)))),2)-pow(sqrt(pow(FSMath::boostEnergy(0.0,0.0,0.0,0.938272,PxPB,PyPB,PzPB,EnPB+0.938272),2)-0.938*0.938)-sqrt(pow(FSMath::boostEnergy(PxP1,PyP1,PzP1,EnP1,PxPB,PyPB,PzPB,EnPB+0.938272),2)-0.938*0.938),2);
        BaryonResonanceMass2 = sqrt((EnP1+EnP2)*(EnP1+EnP2)-(PxP1+PxP2)*(PxP1+PxP2)-(PyP1+PyP2)*(PyP1+PyP2)-(PzP2+PzP1)*(PzP2+PzP1));
        BaryonResonanceMass3 = sqrt((EnP1+EnP3)*(EnP1+EnP3)-(PxP1+PxP3)*(PxP1+PxP3)-(PyP1+PyP3)*(PyP1+PyP3)-(PzP3+PzP1)*(PzP3+PzP1));
        PolarizationDegree = get_polarization(PolarizationAngleLab, EnPB);

        //Calculating van hove coordinates
        VanHoveX = FSMath::vanHoveX(PxP1,PyP1,PzP1,EnP1,PxP2,PyP2,PzP2,EnP2,PxP3,PyP3,PzP3,EnP3);
        VanHoveY = FSMath::vanHoveY(PxP1,PyP1,PzP1,EnP1,PxP2,PyP2,PzP2,EnP2,PxP3,PyP3,PzP3,EnP3);
        VanHoveOmega = FSMath::vanHoveomega(PxP1,PyP1,PzP1,EnP1,PxP2,PyP2,PzP2,EnP2,PxP3,PyP3,PzP3,EnP3);
        
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
        if (isBkg) {
            Weight *= -1;
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
