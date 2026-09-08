#include <stdexcept>
#include <TString.h>
#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

void acceptance_matrix(const char* dataset_name, int maxL, int bootstrap_number) {
    bool is_bootstrap = (bootstrap_number>-1);

    // Sets the number of threads based on the SLURM environment variable
    int nThreads = std::stoi(std::getenv("SLURM_CPUS_PER_TASK"));
    ROOT::EnableImplicitMT(nThreads);

    //Get analysis-specific information

    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    if (is_bootstrap) {
        analysis_info.output_path=TString::Format("%s%s/%i/bootstrap_%i/accMC/",analysis_info.output_path.Data(), dataset_name,maxL, bootstrap_number);
    } else {
        analysis_info.output_path=TString::Format("%s%s/%i/full/accMC/",analysis_info.output_path.Data(), dataset_name,maxL);
    }

    //Create dataframe, filter, and define
    ROOT::RDataFrame generated_dataframe("decayAngles",analysis_info.genMC_path);
    ROOT::RDataFrame accepted_dataframe("decayAngles",analysis_info.accMC_path);
    ROOT::RDF::Experimental::AddProgressBar(generated_dataframe);
    ROOT::RDF::Experimental::AddProgressBar(accepted_dataframe);
    ROOT::RDF::RNode acceptance_matrix_dataframe = accepted_dataframe;
    ROOT::RDF::RNode mass_hist_dataframe = generated_dataframe;
    if (analysis_info.global_cuts.Length() > 0) {
        acceptance_matrix_dataframe = acceptance_matrix_dataframe.Filter(analysis_info.global_cuts.Data());
        mass_hist_dataframe = mass_hist_dataframe.Filter(analysis_info.global_cuts.Data());
    }
    if (analysis_info.non_gen_MC_cuts.Length() > 0) {
        acceptance_matrix_dataframe = acceptance_matrix_dataframe.Filter(analysis_info.non_gen_MC_cuts.Data());
    }
    if (is_bootstrap) {
        acceptance_matrix_dataframe = acceptance_matrix_dataframe.Define("bootstrap_weight","int(gRandom->Poisson(1))");
        mass_hist_dataframe = mass_hist_dataframe.Define("bootstrap_weight","int(gRandom->Poisson(1))");
    }

    //open output file
    std::unique_ptr<TFile> fout(TFile::Open(analysis_info.output_path + "acceptance_matrix.root","RECREATE"));

    //make histograms
    ROOT::RDF::RResultPtr<TH1D> massHistPtr;
    if (is_bootstrap) {
        massHistPtr = mass_hist_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s Mass (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass", "bootstrap_weight");
    } else {
        massHistPtr = mass_hist_dataframe.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s Mass (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass");
    }
    TH1D* massHist = massHistPtr.GetPtr();

    AcceptanceMatrix acceptance_matrix(analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max, maxL, nThreads, massHist, is_bootstrap);
    AcceptanceMatrix* acceptance_matrix_ptr = &acceptance_matrix;
    if (is_bootstrap) {
        if (analysis_info.background_subtract) {
            acceptance_matrix_dataframe.ForeachSlot([acceptance_matrix_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, int bootstrap_weight, double weight){acceptance_matrix_ptr->make_matrix_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, bootstrap_weight, weight);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","bootstrap_weight", "Weight"});
        } else {
            acceptance_matrix_dataframe.ForeachSlot([acceptance_matrix_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, int bootstrap_weight){acceptance_matrix_ptr->make_matrix_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, bootstrap_weight, 1);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","bootstrap_weight"});
        }
        
    } else {
        if (analysis_info.background_subtract) {
            acceptance_matrix_dataframe.ForeachSlot([acceptance_matrix_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double weight){acceptance_matrix_ptr->make_matrix_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, 1, weight);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","Weight"});
        } else {
            acceptance_matrix_dataframe.ForeachSlot([acceptance_matrix_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass){acceptance_matrix_ptr->make_matrix_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, 1, 1);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass"});
        }
    }
    acceptance_matrix.makeMatrix();
    acceptance_matrix.saveMatrices(fout.get());
    fout->Close();
}