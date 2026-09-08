#include <stdexcept>
#include <TString.h>
#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

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
    {
        bin_size = (xMax-xMin)/nBins;
    }
    void make_moments_by_event(unsigned int thread_number, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, int bootstrap_weight, double weight) {
        if (bootstrap_weight == 0) {
            return;
        }
        MomentsArray moments(maxL, helCosTheta, helPhi, PolarizationAngleReac, PolarizationDegree);
        int binNumber = floor((MesonMass-xMin)/bin_size);
        if (binNumber < 0 || binNumber >= nBins) {
            return;
        }
        for (int rowNo = 0; rowNo < dimension; rowNo++) {
            double moment = moments.getMoment(rowNo);
            thread_moment_values[thread_number][binNumber][rowNo] += weight*bootstrap_weight*moment;
            thread_sumw2[thread_number][binNumber][rowNo] += weight*bootstrap_weight*moment*moment;
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
    }
    void saveMoments(TFile* fout) {
        MomentsArray moments(maxL, 0, 0, 0, 0);
        fout->cd();
        TH1D moment_hist_raw("moment_hist_raw", "---", nBins, xMin, xMax);
        TH1D moment_hist("moment_hist", "---", nBins, xMin, xMax);
        moment_hist.SetDirectory(nullptr);
        moment_hist_raw.SetDirectory(nullptr);
        moment_hist.SetStats(0);
        moment_hist_raw.SetStats(0);
        for (int index = 0; index < dimension; index++) {
            std::vector<int> indices = moments.getAlphaLM(index);
            int alpha = indices[0];
            int L = indices[1];
            int M = indices[2];
            TString fig_name_raw = TString::Format("H_%i_%i_%i_raw", alpha, L, M);
            TString fig_name = TString::Format("H_%i_%i_%i", alpha, L, M);
            TString title_raw = TString::Format("H^{%i}(%i%i) for %s; %s (GeV/c^{2}) ; H^{%i}(%i%i)",alpha, L, M, reaction->getReaction().Data(), reaction->massString({2,3}).Data(),alpha, L, M);
            TString title = TString::Format("#LT_{}H^{%i}(%i%i)#GT for %s; %s (GeV/c^{2}) ; #LT_{}H^{%i}(%i%i)#GT",alpha, L, M, reaction->getReaction().Data(), reaction->massString({2,3}).Data(),alpha, L, M);
            for (int bin=0; bin < nBins; bin++) {
                double moment_val = compiled_moments[bin][index];
                double moment_sumw2 = compiled_sumw2[bin][index];
                moment_hist_raw.SetBinContent(bin+1,moment_val);
                moment_hist_raw.SetBinError(bin+1, sqrt(moment_sumw2));
                if (index != 0) {
                    double normalized_moment_val = (compiled_moments[bin][0] == 0) ? -999. : moment_val/compiled_moments[bin][0];
                    double N_events = 2*TMath::Pi() * compiled_moments[bin][0];
                    double uncertainty = 2*TMath::Pi()*sqrt((moment_sumw2)/(N_events*N_events) - (moment_val*moment_val)/(N_events*N_events*N_events));
                    moment_hist.SetBinContent(bin+1, normalized_moment_val);
                    moment_hist.SetBinError(bin+1, uncertainty);
                } else {
                    moment_hist.SetBinContent(bin+1, 1.00);
                    moment_hist.SetBinError(bin+1, 0.00);
                }
                
            }
            moment_hist.SetTitle(title);
            moment_hist_raw.SetTitle(title_raw);
            moment_hist.Write(fig_name);
            moment_hist_raw.Write(fig_name_raw);
            moment_hist.Reset();
            moment_hist_raw.Reset();
        }
    }
};


void calculate_measured_moments(const char* dataset_name, const char* data_type, int maxL, int bootstrap_number) {
    bool is_bootstrap = (bootstrap_number>-1);

    // Sets the number of threads based on the SLURM environment variable
    int nThreads = std::stoi(std::getenv("SLURM_CPUS_PER_TASK"));
    ROOT::EnableImplicitMT(nThreads);

    //Get analysis-specific information

    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    if (is_bootstrap) {
        analysis_info.output_path=TString::Format("%s%s/%i/bootstrap_%i/%s/",analysis_info.output_path.Data(), dataset_name,maxL, bootstrap_number, data_type);
    } else {
        analysis_info.output_path=TString::Format("%s%s/%i/full/%s/",analysis_info.output_path.Data(), dataset_name,maxL, data_type);
    }
    TString fin_path;
    TString additional_cuts = "";
    if (std::string(data_type) == "data") {
        fin_path = analysis_info.data_path;
        additional_cuts = analysis_info.non_gen_MC_cuts;
    } else if (std::string(data_type) == "accMC") {
        fin_path = analysis_info.accMC_path;
        additional_cuts = analysis_info.non_gen_MC_cuts;
    } else if (std::string(data_type) == "genMC") {
        fin_path = analysis_info.genMC_path;
    } else {
        throw std::invalid_argument("Invalid data type. Must be data, genMC, or accMC");
    }

    //Create dataframe, filter, and define
    ROOT::RDataFrame dataframe("decayAngles",fin_path);
    ROOT::RDF::Experimental::AddProgressBar(dataframe);
    ROOT::RDF::RNode moments_dataframe = dataframe;
    if (analysis_info.global_cuts.Length() > 0) {
        moments_dataframe = dataframe.Filter(analysis_info.global_cuts.Data());
    }
    if (additional_cuts.Length() > 0) {
        moments_dataframe = moments_dataframe.Filter(additional_cuts.Data());
    }
    if (is_bootstrap) {
        moments_dataframe = moments_dataframe.Define("bootstrap_weight","int(gRandom->Poisson(1))");
        if (analysis_info.background_subtract) {
            moments_dataframe = moments_dataframe.Define("compositeWeight","Weight*bootstrap_weight");   
        }
    }

    //open output file
    std::unique_ptr<TFile> fout(TFile::Open(analysis_info.output_path + "moments.root","RECREATE"));

    //make histograms
    ROOT::RDF::RResultPtr<TH1D> massHistPtr;
    ROOT::RDF::RResultPtr<TH1D> mass2HistPtr;
    ROOT::RDF::RResultPtr<TH1D> mass3HistPtr;
    ROOT::RDF::RResultPtr<TH2D> angularDistPtr;
    if (is_bootstrap) {
        if (analysis_info.background_subtract) {
            massHistPtr = moments_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s Mass (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass", "compositeWeight");
        } else {
            massHistPtr = moments_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass", "bootstrap_weight");
        }
    } else {
        if (analysis_info.background_subtract) {
            massHistPtr = moments_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass", "Weight");
            mass2HistPtr = moments_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon2Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,2}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass2", "Weight");
            mass3HistPtr = moments_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon3Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,3}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass3", "Weight");
            angularDistPtr = moments_dataframe.Histo2D<double>(ROOT::RDF::TH2DModel("AngularDistribution", TString::Format("Angular Distribution for  %s; cos #theta; #phi", analysis_info.reaction.getReaction().Data()).Data(), analysis_info.n_bins, -1, 1, analysis_info.n_bins, -3.1416, 3.1416), "helCosTheta", "helPhi", "Weight");
        } else {
            massHistPtr = moments_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass");
            mass2HistPtr = moments_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon2Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,2}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass2");
            mass3HistPtr = moments_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon3Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,3}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass3");
            angularDistPtr = moments_dataframe.Histo2D<double>(ROOT::RDF::TH2DModel("AngularDistribution", TString::Format("Angular Distribution for  %s; cos #theta; #phi", analysis_info.reaction.getReaction().Data()).Data(), analysis_info.n_bins, -1, 1, analysis_info.n_bins, -3.1416, 3.1416), "helCosTheta", "helPhi");
        }
    }
    MomentsVector mv(analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max, maxL, nThreads, &analysis_info.reaction, is_bootstrap);
    MomentsVector* mv_ptr = &mv;
    if (is_bootstrap) {
        if (analysis_info.background_subtract) {
            moments_dataframe.ForeachSlot([mv_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, int bootstrap_weight, double weight){mv_ptr->make_moments_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, bootstrap_weight, weight);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","bootstrap_weight","Weight"});
        } else {
            moments_dataframe.ForeachSlot([mv_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, int bootstrap_weight){mv_ptr->make_moments_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, bootstrap_weight, 1);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","bootstrap_weight"});
        }
        
    } else {
        if (analysis_info.background_subtract) {
            moments_dataframe.ForeachSlot([mv_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double weight){mv_ptr->make_moments_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, 1, weight);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","Weight"});
        } else {
            moments_dataframe.ForeachSlot([mv_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass){mv_ptr->make_moments_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, 1, 1);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass"});
        }
    }
    mv.compile_thread_results();
    mv.saveMoments(fout.get());
    TH1D* massHist = massHistPtr.GetPtr();
    massHist->Write("MesonMass");
    if (!is_bootstrap) {
        TH1D* baryonHist2 = mass2HistPtr.GetPtr();
        TH1D* baryonHist3 = mass3HistPtr.GetPtr();
        TH2D* angularDist = angularDistPtr.GetPtr();
        baryonHist2->Write("Baryon2Mass");
        baryonHist3->Write("Baryon3Mass");
        angularDist->Write("AngularDist");
    }
    fout->Close();
}