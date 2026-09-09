#ifndef SHARED_FUNCS
#define SHARED_FUNCS
#include "FSBasic/FSMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TString.h"
#include "TVector3.h"
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "TMatrixD.h"
#include "TSystem.h"


double d_lm0(int L, int M, double costheta)
{
    return pow(-1., M) * sqrt(TMath::Factorial(L - M) / TMath::Factorial(L + M)) * std::assoc_legendre(L, M, costheta);
}

class ReactionSpecs{
    public:
    std::vector<TString> initial_state_particles;
    std::vector<TString> final_state_particles;
    ReactionSpecs()= default;
    ReactionSpecs(std::vector<TString> isp, std::vector<TString> fsp) : initial_state_particles(isp), final_state_particles(fsp) {}
    void addInitialStateParticle(TString isp) {
        initial_state_particles.push_back(isp);
    }
    void addFinalStateParticle(TString fsp) {
        final_state_particles.push_back(fsp);
    }
    TString getFinalStateParticle(int i) {
        return final_state_particles[i-1];
    }
    TString getInitialStateParticle(int i) {
        return initial_state_particles[i-1];
    }
    TString massString(std::vector<int> fsps) {
        TString massString = "";
        for (int i : fsps) {
            massString += getFinalStateParticle(i);
            massString += " ";
        }
        massString += "Mass";
        return massString;
    }
    TString getReaction() {
        TString reactionString = "";
        bool firstParticle = true;
        for (TString isp : initial_state_particles) {
            if (!firstParticle) {
                reactionString += " ";
            } else {
                firstParticle = false;
            }
            reactionString += isp;
        }
        reactionString += " -> ";
        firstParticle = true;
        for (TString fsp : final_state_particles) {
            if (!firstParticle) {
                reactionString += " ";
            } else {
                firstParticle = false;
            }
            reactionString += fsp;
        }
        return reactionString;
    }
};

class AnalysisInfo {
    public:
    AnalysisInfo() = default;
    TString genMC_path;
    TString accMC_path;
    TString data_path;
    TString output_path;
    TString tree_name;

    double meson_mass_min;
    double meson_mass_max;
    TString global_cuts;
    int n_bins;
    TString non_gen_MC_cuts;
    ReactionSpecs reaction;
    bool background_subtract;
};

AnalysisInfo getAnalysisInfo(std::string dataset_name) {
    AnalysisInfo analysis_info;
    analysis_info.tree_name = "decayAngles";
    analysis_info.output_path = "/N/u/rdube/Quartz/work/moments_workflow/results/";
    if (dataset_name == "kpkm_phi") {
            analysis_info.genMC_path = "/home/rdube/scratch/moments_workflow/data/genMC_kpkm_phi_pol*.root";
            analysis_info.accMC_path = "/home/rdube/scratch/moments_workflow/data/accMC_kpkm_phi_pol*.root";
            analysis_info.data_path = "/home/rdube/scratch/moments_workflow/data/data_kpkm_phi_pol*.root";
            analysis_info.background_subtract = false;
            analysis_info.meson_mass_min = 1.005;
            analysis_info.meson_mass_max = 1.035;
            analysis_info.global_cuts= TString::Format("MesonResonanceMass>%.3f&&MesonResonanceMass<%.3f&&abs(BeamEnergy - 8.3)< 0.3&&MandelstamNegT<0.4", analysis_info.meson_mass_min, analysis_info.meson_mass_max);
            analysis_info.non_gen_MC_cuts = "Chi2pipiMinusChi2KK > 10 && Chi2NDF < 5";
            analysis_info.n_bins = 100;
            analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{+}","K^{-}"});
    } else if (dataset_name == "kpkm_highM") {
            analysis_info.genMC_path = "/home/rdube/scratch/moments_workflow/data/genMC_kpkm_highM_pol*.root";
            analysis_info.accMC_path = "/home/rdube/scratch/moments_workflow/data/accMC_kpkm_highM_pol*.root";
            analysis_info.data_path = "/home/rdube/scratch/moments_workflow/data/data_kpkm_highM_pol*.root";
            analysis_info.background_subtract = false;
            analysis_info.meson_mass_min = 1.2;
            analysis_info.meson_mass_max = 2.6;
            analysis_info.global_cuts= TString::Format("MesonResonanceMass>%.3f&&MesonResonanceMass<%.3f&&abs(BeamEnergy - 8.3)< 0.3&&MandelstamNegT<1", analysis_info.meson_mass_min, analysis_info.meson_mass_max);
            analysis_info.non_gen_MC_cuts = "Chi2pipiMinusChi2KK > 10 && Chi2NDF < 5";
            analysis_info.n_bins = 70;
            analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{+}","K^{-}"});
        } else if (dataset_name == "kskl") {
            analysis_info.genMC_path = "/home/rdube/scratch/moments_workflow/data/genMC_kskl_pol*.root";
            analysis_info.accMC_path = "/home/rdube/scratch/moments_workflow/data/accMC_kskl_pol*.root";
            analysis_info.data_path = "/home/rdube/scratch/moments_workflow/data/data_kskl_pol*.root";
            analysis_info.background_subtract = true;
            analysis_info.meson_mass_min = 1.2;
            analysis_info.meson_mass_max = 2.6;
            analysis_info.n_bins = 70;
            analysis_info.global_cuts= "";
            analysis_info.non_gen_MC_cuts = "";
            analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{0}_{S}","K^{0}_{L}"});
        } else if (dataset_name =="kskl_spring_2020") {
            analysis_info.genMC_path = "/N/slate/rdube/moments/kskl/Spring_2020/genMC_kskl_pol*.root";
            analysis_info.accMC_path = "/N/slate/rdube/moments/kskl/Spring_2020/accMC_kskl_pol*.root";
            analysis_info.data_path = "/N/slate/rdube/moments/kskl/Spring_2020/data_kskl_pol*.root";
            analysis_info.background_subtract = true;
            analysis_info.meson_mass_min = 1.2;
            analysis_info.meson_mass_max = 2.56;
            analysis_info.n_bins = 17;
            analysis_info.global_cuts= "";
            analysis_info.non_gen_MC_cuts = "KSFlightSignificance>6&&abs(BeamEnergy-8.3)<0.35&&abs(ProtonVertexZ-65)<23&&NumUnusedTracks==1&&NumUnusedShowers<3&&Chi2NDF<2&&abs(MissingMass-0.5)<0.3";
            analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{0}_{S}","K^{0}_{L}"});
        } else if (dataset_name =="XX") {
            analysis_info.genMC_path = "/home/rdube/scratch/moments_workflow/data/genMC_XX_pol*.root";
            analysis_info.accMC_path = "/home/rdube/scratch/moments_workflow/data/accMC_XX_pol*.root";
            analysis_info.data_path = "/home/rdube/scratch/moments_workflow/data/data_XX_pol*.root";
            analysis_info.meson_mass_min = 1.005;
            analysis_info.meson_mass_max = 1.035;
            analysis_info.global_cuts= "";
            analysis_info.non_gen_MC_cuts = "";
            analysis_info.background_subtract = false;
            analysis_info.n_bins = 100;
            analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{+}","K^{-}"});
        } else if (dataset_name =="kpkm_phi_old") {
            analysis_info.genMC_path = "/home/rdube/scratch/KK_analysis/data/precomputed/genMC_kpkm_phi_pol*.root";
            analysis_info.accMC_path = "/home/rdube/scratch/KK_analysis/data/precomputed/accMC_kpkm_phi_pol*.root";
            analysis_info.data_path = "/home/rdube/scratch/KK_analysis/data/precomputed/data_kpkm_phi_pol*.root";
            analysis_info.meson_mass_min = 1.005;
            analysis_info.meson_mass_max = 1.035;
            analysis_info.n_bins = 100;
            analysis_info.global_cuts= TString::Format("MesonResonanceMass>%.3f&&MesonResonanceMass<%.3f&&abs(BeamEnergy - 8.3)< 0.3&&MandelstamNegT<0.4", analysis_info.meson_mass_min, analysis_info.meson_mass_max);
            analysis_info.non_gen_MC_cuts = "Chi2pipiMinusChi2KK > 10 && Chi2NDF < 5";
            analysis_info.background_subtract = false;
            analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{+}","K^{-}"});
        } else if (dataset_name =="phi_genr8") {
            analysis_info.genMC_path = "/N/u/rdube/Quartz/data/flattened_MC/phi_genr8/final/genMC_kpkm_phi_genr8.root";
            analysis_info.accMC_path = "/N/u/rdube/Quartz/data/flattened_MC/phi_genr8/final/accMC_kpkm_phi_genr8.root";
            analysis_info.data_path = "/N/u/rdube/Quartz/data/data/data_kpkm_phi_pol*.root";
            analysis_info.meson_mass_min = 1.005;
            analysis_info.meson_mass_max = 1.035;
            analysis_info.global_cuts= "abs(BeamEnergy-8.3)<0.3&&MandelstamNegT<0.4";
            analysis_info.non_gen_MC_cuts = "Chi2NDF<5&&Chi2pipiMinusChi2KK>10";
            analysis_info.background_subtract = false;
            analysis_info.n_bins = 30;
            analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{+}","K^{-}"});
        } else if (dataset_name =="phi_gen_amp_v2"){
            analysis_info.genMC_path = "/N/u/rdube/Quartz/data/flattened_MC/phi_gen_amp_v2/final/genMC_kpkm_phi_gen_amp_v2.root";
            analysis_info.accMC_path = "/N/u/rdube/Quartz/data/flattened_MC/phi_gen_amp_v2/final/accMC_kpkm_phi_gen_amp_v2.root";
            analysis_info.data_path = "/N/u/rdube/Quartz/data/data/data_kpkm_phi_pol*.root";
            analysis_info.meson_mass_min = 1.005;
            analysis_info.meson_mass_max = 1.035;
            analysis_info.global_cuts= "abs(BeamEnergy-8.3)<0.3&&MandelstamNegT<0.4";
            analysis_info.non_gen_MC_cuts = "Chi2NDF<5&&Chi2pipiMinusChi2KK>10";
            analysis_info.background_subtract = false;
            analysis_info.n_bins = 30;
            analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{+}","K^{-}"});
        } else {
            std::cout << "Analysis name not found." << std::endl;
    }
    return analysis_info;

}
std::vector<std::pair<TString, TString>> get_columns_to_define(int maxL, std::string frame, bool is_bootstrap)
{
    std::vector<std::pair<TString, TString>> columns_to_define;
    TString frame_to_use;
    if (frame=="helicity" || frame == "hel") {
        frame_to_use = "hel";
    } else {
        frame_to_use = "gj";
    }
    TString bootstrap_weight = "";
    if (is_bootstrap) {
        columns_to_define.push_back({TString("bootstrap_weight"),TString("int(gRandom->Poisson(1))")});
        bootstrap_weight = "bootstrap_weight*";
    }
    for (int l = 0; l <= maxL; l++) {
        for (int m = 0; m <= l; m++) {
            columns_to_define.push_back({ TString::Format("H_0_%i_%i", l, m), TString::Format("%sd_%i_%i_0*cos_%iphi/6.28318530718", bootstrap_weight.Data(),l, m, m) });
            columns_to_define.push_back({ TString::Format("H_1_%i_%i", l, m), TString::Format("%sd_%i_%i_0*cos_%iphi*cos(2*PolarizationAngleReac)/(3.1415926535*PolarizationDegree)", bootstrap_weight.Data(), l, m, m) });
            if (m != 0) {
                columns_to_define.push_back({ TString::Format("H_2_%i_%i", l, m), TString::Format("%sd_%i_%i_0*sin(%i*helPhi)*sin(2*PolarizationAngleReac)/(-3.1415926535*PolarizationDegree)", bootstrap_weight.Data(), l, m, m) });
            }
            
        }
    }
    return columns_to_define;
}

class MomentsArray {
    public:
    int maxL;
    double helCosTheta;
    double helPhi;
    double PolarizationAngleReac;
    double PolarizationDegree;
    std::vector<std::vector<std::vector<int>>> moments_indices;
    std::vector<std::vector<int>> inverse_map;
    std::vector<double> H;
    int dimension;
    MomentsArray(int max_l, double hct, double hp, double par, double pd) 
    : maxL(max_l)
    , helCosTheta(hct)
    , helPhi(hp)
    , PolarizationAngleReac(par)
    , PolarizationDegree(pd) 
    {
        int counter = 0;
        moments_indices.resize(3);
        for (int alpha = 0; alpha < 3; alpha ++) {
            moments_indices[alpha].resize(maxL+1);
            for (int L = 0; L <= maxL; L++) {
                for (int M = 0; M <= L; M++) {
                    if (alpha == 2 && M == 0) {
                        moments_indices[alpha][L].push_back(-999);
                        continue;
                    }
                    moments_indices[alpha][L].push_back(counter);
                    
                    inverse_map.push_back({alpha, L, M});
                    counter+=1;
                }
            }
        }
        dimension = counter;
        H.resize(dimension);
        for (int L = 0; L <= maxL; L++) {
            for (int M=0; M<= L; M++) {
                double d_l_m_0 = d_lm0(L, M, helCosTheta);
                double cosMphi = cos(M*helPhi);
                H[getMomentIndex(0,L,M)] = d_l_m_0*cosMphi/(2*TMath::Pi());
                if (PolarizationDegree <= 1e-6) {
                    H[getMomentIndex(1,L,M)]=0;
                } else {
                    H[getMomentIndex(1,L,M)] = d_l_m_0*cosMphi*cos(2.*PolarizationAngleReac)/(PolarizationDegree*TMath::Pi());
                }
                if (M != 0) {
                    if (PolarizationDegree <= 1e-6) {
                        H[getMomentIndex(2, L, M)]=0;
                    } else {
                        H[getMomentIndex(2, L, M)] = -1.0*d_l_m_0*sin(M*helPhi)*sin(2.*PolarizationAngleReac)/(PolarizationDegree*TMath::Pi());
                    }
                }
            }
        }
    }
    MomentsArray(int max_l) 
    : maxL(max_l)
    , helCosTheta(0)
    , helPhi(0)
    , PolarizationAngleReac(0)
    , PolarizationDegree(0) 
    {
        int counter = 0;
        moments_indices.resize(3);
        for (int alpha = 0; alpha < 3; alpha ++) {
            moments_indices[alpha].resize(maxL+1);
            for (int L = 0; L <= maxL; L++) {
                for (int M = 0; M <= L; M++) {
                    if (alpha == 2 && M == 0) {
                        moments_indices[alpha][L].push_back(-999);
                        continue;
                    }
                    moments_indices[alpha][L].push_back(counter);
                    
                    inverse_map.push_back({alpha, L, M});
                    counter+=1;
                }
            }
        }
        dimension = counter;
    }
    int getMomentIndex(int alpha, int L, int M) {
        if (L > maxL) {
            throw std::invalid_argument(std::string("L is bigger than Lmax"));
        } else if (L < 0) {
            throw std::invalid_argument(std::string("L is negative"));
        } else if (M>L) {
            throw std::invalid_argument(std::string("M is bigger than L"));
        }
        if (alpha == 2 && M==0 ) {
            throw std::out_of_range("Tried to access alpha==2 M==0 index");
        }
        return moments_indices.at(alpha).at(L).at(M);
    }
    double getMoment(int alpha, int L, int M) {
        int moment_to_get = getMomentIndex(alpha,L,M);
        return H.at(moment_to_get);
    }
    double getMoment(int moment_to_get) {
        return H.at(moment_to_get);
    }
    std::vector<int> getAlphaLM(int index) {
        return inverse_map.at(index);
    }
};

void draw_on_same_canvas(std::vector<TH1D*> hists, std::vector<TString> hist_names, TString fig_name, TFile* fout, TCanvas* canvas, TString output_directory_base) {
    canvas->cd();
    std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.2,0.15);
    std::unique_ptr<TLine> zero_line = std::make_unique<TLine>(hists[0]->GetXaxis()->GetXmin(), 0, hists[0]->GetXaxis()->GetXmax(), 0);
    std::vector<Color_t> colors = {kRed+1, kBlue+1, kGreen+3, kMagenta+2, kOrange+7};
    double yMin = 0.;
    double yMax = 0.;
    for (unsigned int i = 0; i < hists.size(); i++) {
        TH1D* hist = hists[i];
        TString hist_name = hist_names[i];
        hist->SetLineColor(colors[i]);
        hist->SetMarkerColor(colors[i]);
        hist->SetStats(0);
        yMin = min(yMin, hist->GetMinimum());
        yMax = max(yMax, hist->GetMaximum());
        legend->AddEntry(hist, hist_name, "lp");
    }

    for (unsigned int i = 0; i < hists.size(); i++) {
        TH1D* hist = hists[i];
        if (i==0) {
            hist->SetMinimum(yMin);
            hist->SetMaximum(yMax);
            hist->Draw("E1 X0");
        } else {
            hist->Draw("SAME E1 X0");
        }
    }
    zero_line->Draw();
    legend->Draw();
    fout->cd();
    canvas->Write(fig_name);
    canvas->SaveAs(TString::Format("%s%s.png", output_directory_base.Data(), fig_name.Data()));
}

void makeDir(TString dirString)
{
    TString dir = gSystem->ExpandPathName(dirString.Data());
    std::cout << "Checking/creating directory: " << dir << std::endl;

    if (gSystem->AccessPathName(dir, kFileExists)) {
        std::cout << "Directory doesn't exist. Creating..." << std::endl;
        int result = gSystem->mkdir(dir, kTRUE);
        if (result != 0) {
            std::cerr << "mkdir failed with code: " << result << std::endl;
        }
    } else {
        std::cout << "Directory already exists." << std::endl;
    }
}
void draw_on_canvas_and_save(TCanvas* canvas, TH1D* histogram, TString fout) {
    canvas->cd();
    std::unique_ptr<TLine> zero_line = std::make_unique<TLine>(histogram->GetXaxis()->GetXmin(),0,histogram->GetXaxis()->GetXmax(),0);
    double y_max = 0.0;
    double y_min = 0.0;
    for (int i=1; i <= histogram->GetNbinsX(); i++) {
        double value = histogram->GetBinContent(i);
        double error = histogram->GetBinError(i);
        y_max = max(y_max, value+error);
        y_min = min(y_min, value-error);
    }
    histogram->SetMinimum(y_min);
    histogram->SetMaximum(y_max);
    histogram->Draw("E1");
    zero_line->Draw("SAME");
    canvas->SaveAs(fout);
    canvas->Clear();
}

class AcceptanceMatrix {
    public:
    int maxL;
    int nThreads;
    TH1D* genMC_hist;
    int nBins;
    double xMin;
    double xMax;
    double bin_size;
    int dimension;
    bool is_bootstrap;
    std::vector<std::vector<TMatrixD>> thread_matrices;
    std::vector<TMatrixD> acceptance_matrices;
    
    AcceptanceMatrix(int n_bins, double x_min, double x_max, int max_l, int n_threads, TH1D* genMC_mass_figure, bool is_bootstrap) 
    : maxL(max_l)
    , nThreads(n_threads)
    , genMC_hist(genMC_mass_figure)
    , nBins(n_bins)
    , xMin(x_min)
    , xMax(x_max)
    , dimension(3*(max_l+1)*(max_l+2)/2 - max_l-1)
    , is_bootstrap(is_bootstrap)
    , thread_matrices(std::vector<std::vector<TMatrixD>>(nThreads, std::vector<TMatrixD>(n_bins,TMatrixD(3*(max_l+1)*(max_l+2)/2- max_l-1,3*(max_l+1)*(max_l+2)/2- max_l-1))))
    , acceptance_matrices(std::vector<TMatrixD>(n_bins, TMatrixD(3*(max_l+1)*(max_l+2)/2 - max_l-1,3*(max_l+1)*(max_l+2)/2 - max_l-1))) {
        bin_size = (xMax-xMin)/nBins;
    }
    void make_matrix_by_event(unsigned int thread_number, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, int bootstrap_weight, double weight) {
        if (bootstrap_weight == 0) {
            return;
        }
        MomentsArray moments(maxL, helCosTheta, helPhi, PolarizationAngleReac, PolarizationDegree);
        int binNumber = floor((MesonMass-xMin)/bin_size)+1;
        if (binNumber < 1 || binNumber > nBins) {
            return;
        }
        int dimension = moments.dimension;
        for (int rowNo = 0; rowNo < dimension; rowNo++) {
            std::vector<int> indices_1 = moments.getAlphaLM(rowNo);
            int alpha = indices_1[0];
            int L = indices_1[1];
            int M = indices_1[2];
            double moment_1 = moments.getMoment(rowNo);
            double transpose_prefactor = matrix_element_prefactor(alpha, L, M, PolarizationDegree);
            for (int colNo = 0; colNo <= rowNo; colNo++) {
                std::vector<int> indices_2 = moments.getAlphaLM(colNo);
                int alphaPrime = indices_2[0];
                int LPrime = indices_2[1];
                int MPrime = indices_2[2];
                double moment_2 = moments.getMoment(colNo);

                double nominal_prefactor = matrix_element_prefactor(alphaPrime, LPrime, MPrime, PolarizationDegree);
                thread_matrices[thread_number][binNumber-1](rowNo,colNo) += weight*bootstrap_weight*nominal_prefactor*moment_1*moment_2;
                if (rowNo != colNo) {
                    thread_matrices[thread_number][binNumber-1](colNo,rowNo) += weight*bootstrap_weight*transpose_prefactor*moment_1*moment_2;
                }
            }
        }
    }
    void makeMatrix() {
        
        std::cout << "finished foreachslot" << std::endl;
        for (int bin=0; bin< nBins; bin++) {
            for (int thread_num=0; thread_num < nThreads; thread_num++) {
                acceptance_matrices[bin]+= thread_matrices[thread_num][bin];
            }
            double n_gen_events_in_bin = genMC_hist->GetBinContent(bin+1);
            acceptance_matrices[bin] *= 1./n_gen_events_in_bin;
        }
    }
    void saveMatrices(TFile* fout) {
        fout->cd();
        for (int bin = 0; bin < nBins; bin++) {
            acceptance_matrices[bin].Write(TString::Format("acceptanceMatrix_%i",bin));
            TH2D matrixHist(acceptance_matrices[bin]);
            matrixHist.SetDirectory(nullptr);
            matrixHist.Write(TString::Format("acceptanceMatrixHistogram_%i",bin));
            TMatrixD inverse_matrix = acceptance_matrices[bin];
            inverse_matrix.Invert();
            inverse_matrix.Write(TString::Format("acceptanceMatrixInverse_%i",bin));
            TH2D matrixHistInverse(inverse_matrix);
            matrixHistInverse.SetDirectory(nullptr);
            matrixHistInverse.Write(TString::Format("acceptanceMatrixHistogramInverse_%i",bin));
        }
    }
    double matrix_element_prefactor(int alpha, int L, int M, double PolarizationDegree) {
        double prefactor = 2.0*TMath::Pi()*TMath::Pi()*(2.0*L+1);
        if (M != 0) {
            prefactor *= 2.0;
        }
        if (alpha == 0) {
            prefactor *= 2.0;
        } else {
            prefactor *= PolarizationDegree*PolarizationDegree;
        }
        return prefactor;
    }
};

#endif
