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
#include "THStack.h"


double d_lm0(int L, int M, double costheta)
{
    return ((M & 1) ? -1.0 : 1.0) * sqrt(TMath::Factorial(L - M) / TMath::Factorial(L + M)) * std::assoc_legendre(L, M, costheta);
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
    TString input_path;
    std::vector<TString> independent_files;
    TString output_path;
    TString tree_name;

    std::map<TString,double> luminosities;
    bool use_luminosities;

    double meson_mass_min;
    double meson_mass_max;
    TString global_cuts;
    int n_bins;
    TString non_gen_MC_cuts;
    ReactionSpecs reaction;
    bool background_subtract;

    double branching_fraction;

};

AnalysisInfo getAnalysisInfo(std::string dataset_name) {
    AnalysisInfo analysis_info;
    analysis_info.tree_name = "decayAngles";
    analysis_info.output_path = "/N/u/rdube/Quartz/work/moments_workflow/results/";
    analysis_info.branching_fraction=1.0;
    if (dataset_name == "kpkm_highM") {
        analysis_info.independent_files  = {"kpkm_pol0","kpkm_pol45","kpkm_pol90","kpkm_pol135"};
        analysis_info.input_path = "/N/slate/rdube/moments/kpkm/highM";
        analysis_info.background_subtract = true;
        analysis_info.meson_mass_min = 1.1;
        analysis_info.meson_mass_max = 2.6;
        analysis_info.global_cuts= TString::Format("abs(BeamEnergy - 8.3)< 0.3&&abs(MandelstamNegT-0.45)<0.25&&MesonResonanceMass>%.3f&&MesonResonanceMass<%.3f", analysis_info.meson_mass_min, analysis_info.meson_mass_max);
        analysis_info.non_gen_MC_cuts = "Chi2pipiMinusChi2KK>10&&Chi2NDF<5&&NumUnusedTracks==0&&NumUnusedShowers==0";
        analysis_info.n_bins = 30;
        analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{+}","K^{-}"});
        analysis_info.use_luminosities = false;
    } else if (dataset_name =="kskl_Spring_2020") {
        analysis_info.independent_files  = {"kskl_pol0","kskl_pol45","kskl_pol90","kskl_pol135"};
        analysis_info.luminosities = std::map<TString,double>({{"kskl_pol0",33.2179},{"kskl_pol45",32.9354},{"kskl_pol90",31.1341},{"kskl_pol135",31.9134}});
        analysis_info.input_path = "/N/slate/rdube/moments/kskl/Spring_2020";
        analysis_info.background_subtract = true;
        analysis_info.meson_mass_min = 1.2;
        analysis_info.meson_mass_max = 2.6;
        analysis_info.n_bins = 28;
        analysis_info.global_cuts= TString::Format("abs(BeamEnergy-8.3)<0.3&&abs(MandelstamNegT-0.45)<0.25&&MesonResonanceMass>%.3f&&MesonResonanceMass<%.3f", analysis_info.meson_mass_min, analysis_info.meson_mass_max);
        analysis_info.non_gen_MC_cuts = "abs(ProtonVertexZ-65)<13&&NumUnusedTracks==0&&NumUnusedShowers<3&&KSFlightSignificance>6&&Chi2NDF<2&&abs(MissingMass-0.5)<0.2";
        analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{0}_{S}","K^{0}_{L}"});
        analysis_info.use_luminosities = true;
        analysis_info.branching_fraction = 0.6920;
    } else if (dataset_name =="kskl_GlueX_I") {
        analysis_info.independent_files  = {"kskl_fall_2018_pol0","kskl_fall_2018_pol45","kskl_fall_2018_pol90","kskl_fall_2018_pol135","kskl_spring_2017_pol0","kskl_spring_2017_pol45","kskl_spring_2017_pol90","kskl_spring_2017_pol135","kskl_spring_2018_pol0","kskl_spring_2018_pol45","kskl_spring_2018_pol90","kskl_spring_2018_pol135"};
        analysis_info.input_path = "/N/slate/rdube/moments/kskl/GlueX_I";
        analysis_info.background_subtract = true;
        analysis_info.meson_mass_min = 1.1;
        analysis_info.meson_mass_max = 2.6;
        analysis_info.n_bins = 30;
        analysis_info.global_cuts= TString::Format("abs(BeamEnergy-8.5)<0.3&&abs(MandelstamNegT-0.45)<0.25&&MesonResonanceMass>%.3f&&MesonResonanceMass<%.3f", analysis_info.meson_mass_min, analysis_info.meson_mass_max);
        analysis_info.non_gen_MC_cuts = "abs(ProtonVertexZ-65)<13&&NumUnusedTracks==0&&NumUnusedShowers<3&&KSFlightSignificance>6&&Chi2NDF<2&&abs(MissingMass-0.5)<0.2";
        analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{0}_{S}","K^{0}_{L}"});
        analysis_info.luminosities = std::map<TString,double>({{"kskl_spring_2017_pol0",5.1977},{"kskl_spring_2017_pol45",4.95599},{"kskl_spring_2017_pol90",5.42077},{"kskl_spring_2017_pol135",5.09782},{"kskl_spring_2018_pol0",15.0207},{"kskl_spring_2018_pol45",14.6655},{"kskl_spring_2018_pol90",14.472},{"kskl_spring_2018_pol135",15.2797},{"kskl_fall_2018_pol0",10.1941},{"kskl_fall_2018_pol45",9.28027},{"kskl_fall_2018_pol90",9.36105},{"kskl_fall_2018_pol135",9.01873}});
        analysis_info.use_luminosities = true;
        analysis_info.branching_fraction = 0.6920;
    } else if (dataset_name =="kskl_GlueX_I_combined_runperiods") {
        analysis_info.independent_files  = {"kskl_*_pol0","kskl_*_pol45","kskl_*_pol90","kskl_*_pol135"};
        analysis_info.luminosities = std::map<TString,double>({{"kskl_*_pol0",5.1977+15.0207+10.1941},{"kskl_*_pol45",4.95599+14.6655+9.28027},{"kskl_*_pol90",5.42077+14.472+9.36105},{"kskl_*_pol135",5.09782+15.2797+9.01873}});
        analysis_info.input_path = "/N/slate/rdube/moments/kskl/GlueX_I";
        analysis_info.background_subtract = true;
        analysis_info.meson_mass_min = 1.1;
        analysis_info.meson_mass_max = 2.6;
        analysis_info.n_bins = 30;
        analysis_info.global_cuts= TString::Format("abs(BeamEnergy-8.5)<0.3&&abs(MandelstamNegT-0.45)<0.25&&MesonResonanceMass>%.3f&&MesonResonanceMass<%.3f", analysis_info.meson_mass_min, analysis_info.meson_mass_max);
        analysis_info.non_gen_MC_cuts = "abs(ProtonVertexZ-65)<13&&NumUnusedTracks==0&&NumUnusedShowers<3&&KSFlightSignificance>6&&Chi2NDF<2&&abs(MissingMass-0.5)<0.2";
        analysis_info.reaction = ReactionSpecs({"#gamma","p"},{"p","K^{0}_{S}","K^{0}_{L}"});
        analysis_info.use_luminosities = true;
        analysis_info.branching_fraction = 0.6920;
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
    std::vector<std::array<int,3>> inverse_map;
    std::vector<double> H;
    int dimension;
    MomentsArray(int max_l, double hct, double hp, double par, double pd) 
    : maxL(max_l)
    , helCosTheta(hct)
    , helPhi(hp)
    , PolarizationAngleReac(par)
    , PolarizationDegree(pd) 
    {
        dimension = maxL*(3*maxL+7)/2 + 2;
        inverse_map.reserve(dimension);
        int counter = 0;
        moments_indices.resize(3);
        for (int alpha = 0; alpha < 3; alpha ++) {
            moments_indices[alpha].resize(maxL+1);
            for (int L = 0; L <= maxL; L++) {
                moments_indices[alpha][L].reserve(L+1);
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
        H.resize(dimension);
        for (int L = 0; L <= maxL; L++) {
            for (int M=0; M<= L; M++) {
                double d_l_m_0 = d_lm0(L, M, helCosTheta);
                double cosMphi = cos(M*helPhi);
                H[getMomentIndex(0,L,M)] = d_l_m_0*cosMphi/(2*TMath::Pi());
                H[getMomentIndex(1,L,M)] = d_l_m_0*cosMphi*cos(2.*PolarizationAngleReac)/(PolarizationDegree*TMath::Pi());
                if (M != 0) {
                    H[getMomentIndex(2, L, M)] = -1.0*d_l_m_0*sin(M*helPhi)*sin(2.*PolarizationAngleReac)/(PolarizationDegree*TMath::Pi());
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
        dimension = maxL*(3*maxL+7)/2 + 2;
        inverse_map.reserve(dimension);
        moments_indices.resize(3);
        for (int alpha = 0; alpha < 3; alpha ++) {
            moments_indices[alpha].resize(maxL+1);
            for (int L = 0; L <= maxL; L++) {
                moments_indices[alpha][L].reserve(L+1);
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
    }
    int getMomentIndex(int alpha, int L, int M) {
        return moments_indices[alpha][L][M];
    }
    double getMoment(int alpha, int L, int M) {
        int moment_to_get = getMomentIndex(alpha,L,M);
        return H[moment_to_get];
    }
    double getMoment(int moment_to_get) {
        return H[moment_to_get];
    }
    const std::array<int,3>& getAlphaLM(int index) const {
        return inverse_map[index];
    }
    void Reset(double hct, double hp, double par, double pd) 
    {
        helCosTheta = hct;
        helPhi = hp;
        PolarizationAngleReac = par;
        PolarizationDegree = pd;
        H.resize(dimension);
        for (int L = 0; L <= maxL; L++) {
            for (int M=0; M<= L; M++) {
                double d_l_m_0 = d_lm0(L, M, helCosTheta);
                double cosMphi = cos(M*helPhi);
                H[getMomentIndex(0,L,M)] = d_l_m_0*cosMphi/(2*TMath::Pi());
                H[getMomentIndex(1,L,M)] = d_l_m_0*cosMphi*cos(2.*PolarizationAngleReac)/(PolarizationDegree*TMath::Pi());
                if (M != 0) {
                    H[getMomentIndex(2, L, M)] = -1.0*d_l_m_0*sin(M*helPhi)*sin(2.*PolarizationAngleReac)/(PolarizationDegree*TMath::Pi());
                }
            }
        }
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

void draw_on_same_canvas_with_statistical_errors(std::vector<TH1D*> hists, std::vector<TH1D*> hists_statistical,std::vector<TString> hist_names, TString fig_name, TFile* fout, TCanvas* canvas, TString output_directory_base) {
    canvas->cd();
    std::unique_ptr<TLegend> legend = std::make_unique<TLegend>(0.2,0.15);
    std::unique_ptr<TLine> zero_line = std::make_unique<TLine>(hists[0]->GetXaxis()->GetXmin(), 0, hists[0]->GetXaxis()->GetXmax(), 0);
    std::vector<Color_t> colors = {kRed+1, kBlue+1, kGreen+3, kMagenta+2, kOrange+7};
    double yMin = 0.;
    double yMax = 0.;
    for (unsigned int i = 0; i < hists.size(); i++) {
        TH1D* hist = hists[i];
        TH1D* hist_statistical = hists_statistical[i];
        TString hist_name = hist_names[i];
        hist->SetLineColor(colors[i]);
        hist->SetMarkerColor(colors[i]);
        hist->SetStats(0);
        hist_statistical->SetLineColor(colors[i]);
        hist_statistical->SetMarkerColor(colors[i]);
        hist_statistical->SetStats(0);
        yMin = min(yMin, hist->GetMinimum());
        yMax = max(yMax, hist->GetMaximum());
        yMin = min(yMin, hist_statistical->GetMinimum());
        yMax = max(yMax, hist_statistical->GetMaximum());
        legend->AddEntry(hist, hist_name, "lp");
    }

    for (unsigned int i = 0; i < hists.size(); i++) {
        TH1D* hist = hists[i];
        TH1D* hist_statistical = hists_statistical[i];
        if (i==0) {
            hist->SetMinimum(yMin);
            hist->SetMaximum(yMax);
            hist->Draw("E1 X0");
        } else {
            hist->Draw("SAME E1 X0");
        }
        hist_statistical->Draw("SAME E1 X0");
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
    std::vector<std::vector<std::vector<double>>> thread_matrices;
    std::vector<TMatrixD> acceptance_matrices;
    std::vector<MomentsArray> moments_arrays;
    
    AcceptanceMatrix(int n_bins, double x_min, double x_max, int max_l, int n_threads, TH1D* genMC_mass_figure, bool is_bootstrap) 
    : maxL(max_l)
    , nThreads(n_threads)
    , genMC_hist(genMC_mass_figure)
    , nBins(n_bins)
    , xMin(x_min)
    , xMax(x_max)
    , dimension(3*(max_l+1)*(max_l+2)/2 - max_l-1)
    , is_bootstrap(is_bootstrap)
    , thread_matrices(std::vector<std::vector<std::vector<double>>>(nThreads, std::vector<std::vector<double>>(n_bins,std::vector<double>(dimension*dimension))))
    , acceptance_matrices(std::vector<TMatrixD>(n_bins, TMatrixD(3*(max_l+1)*(max_l+2)/2 - max_l-1,3*(max_l+1)*(max_l+2)/2 - max_l-1)))
    , moments_arrays(std::vector<MomentsArray>(nThreads,MomentsArray(max_l))) {
        bin_size = (xMax-xMin)/nBins;
    }
    void make_matrix_by_event(unsigned int thread_number, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double full_weight) {
        if (full_weight == 0) {
            return;
        }
        moments_arrays[thread_number].Reset(helCosTheta, helPhi, PolarizationAngleReac, PolarizationDegree);
        int binNumber = floor((MesonMass-xMin)/bin_size)+1;
        if (binNumber < 1 || binNumber > nBins) {
            return;
        }
        int dimension = moments_arrays[thread_number].dimension;
        for (int rowNo = 0; rowNo < dimension; rowNo++) {
            const std::array<int,3>& indices_1 = moments_arrays[thread_number].getAlphaLM(rowNo);
            int alpha = indices_1[0];
            int L = indices_1[1];
            int M = indices_1[2];
            double moment_1 = moments_arrays[thread_number].getMoment(rowNo);
            double transpose_prefactor = matrix_element_prefactor(alpha, L, M, PolarizationDegree);
            for (int colNo = 0; colNo <= rowNo; colNo++) {
                const std::array<int,3>& indices_2 = moments_arrays[thread_number].getAlphaLM(colNo);
                int alphaPrime = indices_2[0];
                int LPrime = indices_2[1];
                int MPrime = indices_2[2];
                double moment_2 = moments_arrays[thread_number].getMoment(colNo);

                double nominal_prefactor = matrix_element_prefactor(alphaPrime, LPrime, MPrime, PolarizationDegree);
                thread_matrices[thread_number][binNumber-1][rowNo*dimension+colNo] += full_weight*nominal_prefactor*moment_1*moment_2;
                if (rowNo != colNo) {
                    thread_matrices[thread_number][binNumber-1][colNo*dimension+rowNo] += full_weight*transpose_prefactor*moment_1*moment_2;
                }
            }
        }
    }
    void makeMatrix() {
        std::cout << "finished foreachslot" << std::endl;
        for (int bin=0; bin< nBins; bin++) {
            for (int thread_num=0; thread_num < nThreads; thread_num++) {
                for (int rowNo = 0; rowNo < dimension; rowNo++) {
                    for (int colNo = 0; colNo < dimension; colNo++) {
                        acceptance_matrices[bin](rowNo,colNo)+=thread_matrices[thread_num][bin][rowNo*dimension+colNo];
                    }
                }
            }
            double n_gen_events_in_bin = genMC_hist->GetBinContent(bin+1);
            acceptance_matrices[bin] *= 1./n_gen_events_in_bin;
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

class MomentsVector {
    public:
    int maxL;
    int nThreads;
    int nBins;
    double xMin;
    double xMax;
    double bin_size;
    int dimension;
    bool isBootstrap;
    std::vector<std::vector<std::vector<double>>> thread_moment_values;
    std::vector<std::vector<std::vector<double>>> thread_sumw2;
    std::vector<std::vector<double>> compiled_moments;
    std::vector<std::vector<double>> compiled_sumw2;
    ReactionSpecs* reaction;
    std::vector<std::unique_ptr<TH1D>> moments_histograms;
    std::vector<std::unique_ptr<TH1D>> moments_histograms_raw;
    std::vector<MomentsArray> moments_arrays;
    
    MomentsVector(double n_bins, double x_min, double x_max, int max_l, int n_threads, ReactionSpecs* reac, bool is_bootstrap) 
    : maxL(max_l)
    , nThreads(n_threads)
    , nBins(n_bins)
    , xMin(x_min)
    , xMax(x_max)
    , dimension(3*(max_l+1)*(max_l+2)/2 - max_l-1)
    , isBootstrap(is_bootstrap)
    , thread_moment_values(std::vector<std::vector<std::vector<double>>>(n_threads,std::vector<std::vector<double>>(n_bins, std::vector<double>(3*(max_l+1)*(max_l+2)/2 - max_l-1))))
    , thread_sumw2(std::vector<std::vector<std::vector<double>>>(n_threads,std::vector<std::vector<double>>(n_bins, std::vector<double>(3*(max_l+1)*(max_l+2)/2 - max_l-1))))
    , compiled_moments(std::vector<std::vector<double>>(n_bins, std::vector<double>(3*(max_l+1)*(max_l+2)/2 - max_l-1)))
    , compiled_sumw2(std::vector<std::vector<double>>(n_bins, std::vector<double>(3*(max_l+1)*(max_l+2)/2 - max_l-1)))
    , reaction(reac)
    , moments_arrays(std::vector<MomentsArray>(n_threads,MomentsArray(max_l)))
    {
        bin_size = (xMax-xMin)/nBins;
    }
    void make_moments_by_event(unsigned int thread_number, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double full_weight) {
        if (full_weight == 0) {
            return;
        }
        moments_arrays[thread_number].Reset(helCosTheta, helPhi, PolarizationAngleReac, PolarizationDegree);
        int binNumber = floor((MesonMass-xMin)/bin_size);
        if (binNumber < 0 || binNumber >= nBins) {
            return;
        }
        for (int rowNo = 0; rowNo < dimension; rowNo++) {
            double moment = moments_arrays[thread_number].getMoment(rowNo);
            thread_moment_values[thread_number][binNumber][rowNo] += full_weight*moment;
            thread_sumw2[thread_number][binNumber][rowNo] += full_weight*moment*moment;
        }
    }
    void compile_thread_results() {
        for (int bin=0; bin< nBins; bin++) {
            for (int thread_num=0; thread_num < nThreads; thread_num++) {
                for (int index = 0; index < dimension; index++) {
                    compiled_moments[bin][index]+= thread_moment_values[thread_num][bin][index];
                    compiled_sumw2[bin][index] += thread_sumw2[thread_num][bin][index];
                }
            }
        }
        for (int index = 0; index < dimension; index++) {
            const std::array<int,3>& indices = moments_arrays[0].getAlphaLM(index);
            int alpha = indices[0];
            int L = indices[1];
            int M = indices[2];
            TString fig_name_raw = TString::Format("H_%i_%i_%i_raw", alpha, L, M);
            TString fig_name = TString::Format("H_%i_%i_%i", alpha, L, M);
            TString title_raw = TString::Format("H^{%i}(%i%i) for %s; %s (GeV/c^{2}) ; H^{%i}(%i%i)",alpha, L, M, reaction->getReaction().Data(), reaction->massString({2,3}).Data(),alpha, L, M);
            TString title = TString::Format("#LT_{}H^{%i}(%i%i)#GT for %s; %s (GeV/c^{2}) ; #LT_{}H^{%i}(%i%i)#GT",alpha, L, M, reaction->getReaction().Data(), reaction->massString({2,3}).Data(),alpha, L, M);
            moments_histograms.push_back(std::make_unique<TH1D>(fig_name,title,nBins, xMin, xMax));
            moments_histograms_raw.push_back(std::make_unique<TH1D>(fig_name_raw,title_raw,nBins, xMin, xMax));
            moments_histograms.back()->SetStats(0);
            moments_histograms_raw.back()->SetStats(0);
            moments_histograms.back()->SetDirectory(nullptr);
            moments_histograms_raw.back()->SetDirectory(nullptr);
            for (int bin=0; bin < nBins; bin++) {
                double moment_val = compiled_moments[bin][index];
                double moment_sumw2 = compiled_sumw2[bin][index];
                moments_histograms_raw.back()->SetBinContent(bin+1,moment_val);
                moments_histograms_raw.back()->SetBinError(bin+1, sqrt(moment_sumw2));
                if (index != 0) {
                    double normalized_moment_val = (compiled_moments[bin][0] == 0) ? -999. : moment_val/compiled_moments[bin][0];
                    double N_events = 2*TMath::Pi() * compiled_moments[bin][0];
                    double uncertainty = 2*TMath::Pi()*sqrt((moment_sumw2)/(N_events*N_events) - (moment_val*moment_val)/(N_events*N_events*N_events));
                    moments_histograms.back()->SetBinContent(bin+1, normalized_moment_val);
                    moments_histograms.back()->SetBinError(bin+1, uncertainty);
                } else {
                    moments_histograms.back()->SetBinContent(bin+1, 1.00);
                    moments_histograms.back()->SetBinError(bin+1, 0.00);
                }
            }
        }
    }
};

void compare_results_files(std::vector<TString> file_paths, std::vector<TString> names, TString output_dir, int L_max=6) {
    std::vector<Color_t> colors = {kRed+1,kBlue+1, kGreen+3, kMagenta+2, kYellow+2};
    std::vector<std::unique_ptr<TFile>> files_to_compare;
    for (TString file_name : file_paths) {
        files_to_compare.push_back(std::unique_ptr<TFile>(TFile::Open(file_name,"READ")));
    }
    TCanvas canvas("canvas","canvas",700,500);
    TCanvas canvas_raw("canvas_raw","canvas_raw",700,500);
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L = 0; L <= L_max; L++) {
            for (int M = 0; M <= L; M++) {
                if (alpha == 2 && M == 0) {
                    continue;
                }
                THStack moments("moments_stack",TString::Format("#LT_{}H^{%i}(%i%i)#GT for #gamma p -> p K_{S} K_{L}",alpha, L, M).Data());
                
                THStack moments_raw("moments_stack_raw",TString::Format("H^{%i}(%i%i) for #gamma p -> p K_{S} K_{L}",alpha, L, M).Data());
                

                TLegend legend(0.2,0.15);
                TLegend legend_raw(0.2,0.15);
                bool any_present = false;
                for (unsigned int i = 0; i < files_to_compare.size(); i++) {
                    TH1D* moment_histogram = (TH1D*)files_to_compare[i]->Get(TString::Format("H_%i_%i_%i_acceptance_corrected",alpha,L,M));
                    if (!moment_histogram) {
                        continue;
                    }
                    moment_histogram->SetLineColor(colors[i]);
                    moment_histogram->SetMarkerColor(colors[i]);
                    legend.AddEntry(moment_histogram, names[i], "lp");
                    moments.Add(moment_histogram);
                    TH1D* moment_histogram_raw = (TH1D*)files_to_compare[i]->Get(TString::Format("H_%i_%i_%i_raw_acceptance_corrected",alpha,L,M));
                    moment_histogram_raw->SetLineColor(colors[i]);
                    moment_histogram_raw->SetMarkerColor(colors[i]);
                    legend_raw.AddEntry(moment_histogram_raw, names[i], "lp");
                    moments_raw.Add(moment_histogram_raw);
                    any_present = true;
                }
                if (!any_present) {
                    continue;
                }
                canvas.cd();
                moments.Draw("nostack E1 X0");
                TLine line(moments.GetXaxis()->GetXmin(),0,moments.GetXaxis()->GetXmax(),0);
                moments.GetXaxis()->SetTitle("K_{S} K_{L} Mass (GeV/c^{2})");
                moments.GetYaxis()->SetTitle(TString::Format("#LT_{}H^{%i}(%i%i)#GT",alpha,L,M).Data());
                moments.Draw("nostack noclear E1 X0");
                line.Draw();
                legend.Draw();
                canvas.SaveAs(TString::Format("%s/H_%i_%i_%i.png",output_dir.Data(), alpha, L, M));

                canvas_raw.cd();
                moments_raw.Draw("nostack E1 X0");
                moments_raw.GetXaxis()->SetTitle("K_{S} K_{L} Mass (GeV/c^{2})");
                moments_raw.GetYaxis()->SetTitle(TString::Format("H^{%i}(%i%i)",alpha,L,M).Data());
                moments_raw.Draw("nostack noclear E1 X0");
                line.Draw();
                legend_raw.Draw();
                canvas_raw.SaveAs(TString::Format("%s/H_%i_%i_%i_raw.png",output_dir.Data(), alpha, L, M));
            }
        }
    }
    for (unsigned int i = 0; i < files_to_compare.size(); i++) {
        files_to_compare[i]->Close();
    }
}

#endif
