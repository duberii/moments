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
    std::vector<std::unique_ptr<TH1D>> moments_histograms;
    std::vector<std::unique_ptr<TH1D>> moments_histograms_raw;
    MomentsArray moments_array;
    
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
    , moments_array(max_l)
    {
        bin_size = (xMax-xMin)/nBins;
    }
    void make_moments_by_event(unsigned int thread_number, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double full_weight) {
        if (full_weight == 0) {
            return;
        }
        MomentsArray moments_array_local(maxL, helCosTheta, helPhi, PolarizationAngleReac, PolarizationDegree);
        int binNumber = floor((MesonMass-xMin)/bin_size);
        if (binNumber < 0 || binNumber >= nBins) {
            return;
        }
        for (int rowNo = 0; rowNo < dimension; rowNo++) {
            double moment = moments_array_local.getMoment(rowNo);
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
            std::vector<int> indices = moments_array.getAlphaLM(index);
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

void single_step_calculate_moments(const char* dataset_name, int maxL, int bootstrap_number) {
    bool is_bootstrap = (bootstrap_number>-1);

    // Sets the number of threads based on the SLURM environment variable
    int nThreads = std::stoi(std::getenv("SLURM_CPUS_PER_TASK"));
    ROOT::EnableImplicitMT(nThreads);

    //Get analysis-specific information
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    if (is_bootstrap) {
        analysis_info.output_path=TString::Format("%s%s/%i/bootstrap_%i",analysis_info.output_path.Data(), dataset_name,maxL, bootstrap_number);
    } else {
        analysis_info.output_path=TString::Format("%s%s/%i/full",analysis_info.output_path.Data(), dataset_name,maxL);
    }

    //Create dataframe, filter, and define
    ROOT::RDataFrame dataframe_data("decayAngles",analysis_info.data_path);
    ROOT::RDataFrame dataframe_accMC("decayAngles",analysis_info.accMC_path);
    ROOT::RDataFrame dataframe_genMC("decayAngles",analysis_info.genMC_path);
    ROOT::RDF::Experimental::AddProgressBar(dataframe_data);
    ROOT::RDF::Experimental::AddProgressBar(dataframe_accMC);
    ROOT::RDF::Experimental::AddProgressBar(dataframe_genMC);
    ROOT::RDF::RNode moments_df_data = dataframe_data;
    ROOT::RDF::RNode moments_df_accMC = dataframe_accMC;
    ROOT::RDF::RNode moments_df_genMC = dataframe_genMC;
    if (analysis_info.global_cuts.Length() > 0) {
        moments_df_data = moments_df_data.Filter(analysis_info.global_cuts.Data());
        moments_df_accMC = moments_df_accMC.Filter(analysis_info.global_cuts.Data());
        moments_df_genMC = moments_df_genMC.Filter(analysis_info.global_cuts.Data());
    }
    if (analysis_info.non_gen_MC_cuts.Length() > 0) {
        moments_df_data = moments_df_data.Filter(analysis_info.non_gen_MC_cuts.Data());
        moments_df_accMC = moments_df_accMC.Filter(analysis_info.non_gen_MC_cuts.Data());
    }
    TString full_weight_string = "1";
    if (is_bootstrap) {
        full_weight_string += "*int(gRandom->Poisson(1))";
    }
    if (analysis_info.background_subtract) {
        full_weight_string += "*Weight";
    }
    moments_df_data = moments_df_data.Define("FullWeight",full_weight_string.Data());
    moments_df_accMC = moments_df_accMC.Define("FullWeight",full_weight_string.Data());   
    moments_df_genMC = moments_df_genMC.Define("FullWeight",full_weight_string.Data());

    //make histograms
    ROOT::RDF::RResultPtr<TH1D> massHistPtr_data = moments_df_data.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s Mass (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass", "FullWeight");
    ROOT::RDF::RResultPtr<TH1D> mass2HistPtr_data = moments_df_data.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon2Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,2}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass2", "FullWeight");
    ROOT::RDF::RResultPtr<TH1D> mass3HistPtr_data = moments_df_data.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon3Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,3}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass3", "FullWeight");
    ROOT::RDF::RResultPtr<TH2D> angularDistPtr_data = moments_df_data.Histo2D<double>(ROOT::RDF::TH2DModel("AngularDistribution", TString::Format("Angular Distribution for  %s; cos #theta; #phi", analysis_info.reaction.getReaction().Data()).Data(), analysis_info.n_bins, -1, 1, analysis_info.n_bins, -3.1416, 3.1416), "helCosTheta", "helPhi", "FullWeight");

    ROOT::RDF::RResultPtr<TH1D> massHistPtr_accMC = moments_df_accMC.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s Mass (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass", "FullWeight");
    ROOT::RDF::RResultPtr<TH1D> mass2HistPtr_accMC = moments_df_accMC.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon2Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,2}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass2", "FullWeight");
    ROOT::RDF::RResultPtr<TH1D> mass3HistPtr_accMC = moments_df_accMC.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon3Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,3}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass3", "FullWeight");
    ROOT::RDF::RResultPtr<TH2D> angularDistPtr_accMC = moments_df_accMC.Histo2D<double>(ROOT::RDF::TH2DModel("AngularDistribution", TString::Format("Angular Distribution for  %s; cos #theta; #phi", analysis_info.reaction.getReaction().Data()).Data(), analysis_info.n_bins, -1, 1, analysis_info.n_bins, -3.1416, 3.1416), "helCosTheta", "helPhi", "FullWeight");

    ROOT::RDF::RResultPtr<TH1D> massHistPtr_genMC = moments_df_genMC.Histo1D<double>(ROOT::RDF::TH1DModel("MesonMass", TString::Format("Meson Mass for %s; %s Mass (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({2,3}).Data()).Data(), analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max), "MesonResonanceMass", "FullWeight");
    ROOT::RDF::RResultPtr<TH1D> mass2HistPtr_genMC = moments_df_genMC.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon2Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,2}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass2", "FullWeight");
    ROOT::RDF::RResultPtr<TH1D> mass3HistPtr_genMC = moments_df_genMC.Histo1D<double>(ROOT::RDF::TH1DModel("Baryon3Mass", TString::Format("Baryon Mass for %s; %s (GeV/c^{2}); Counts", analysis_info.reaction.getReaction().Data(),analysis_info.reaction.massString({1,3}).Data()).Data(), analysis_info.n_bins, 1.4, 3.5), "BaryonResonanceMass3", "FullWeight");
    ROOT::RDF::RResultPtr<TH2D> angularDistPtr_genMC = moments_df_genMC.Histo2D<double>(ROOT::RDF::TH2DModel("AngularDistribution", TString::Format("Angular Distribution for  %s; cos #theta; #phi", analysis_info.reaction.getReaction().Data()).Data(), analysis_info.n_bins, -1, 1, analysis_info.n_bins, -3.1416, 3.1416), "helCosTheta", "helPhi", "FullWeight");

    MomentsVector mv_data(analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max, maxL, nThreads, &analysis_info.reaction, is_bootstrap);
    MomentsVector mv_accMC(analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max, maxL, nThreads, &analysis_info.reaction, is_bootstrap);
    MomentsVector mv_genMC(analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max, maxL, nThreads, &analysis_info.reaction, is_bootstrap);
    MomentsVector* mv_ptr_data = &mv_data;
    MomentsVector* mv_ptr_accMC = &mv_accMC;
    MomentsVector* mv_ptr_genMC = &mv_genMC;

    //open output file
    moments_df_data.ForeachSlot([mv_ptr_data](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double full_weight){mv_ptr_data->make_moments_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, full_weight);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","FullWeight"});
    mv_data.compile_thread_results();
    std::unique_ptr<TFile> fout_data(TFile::Open(TString::Format("%s/data/moments.root",analysis_info.output_path.Data()),"RECREATE"));
    fout_data->cd();
    if (!is_bootstrap) {
        TH1D* massHist_data = massHistPtr_data.GetPtr();
        TH1D* baryonHist2_data = mass2HistPtr_data.GetPtr();
        TH1D* baryonHist3_data = mass3HistPtr_data.GetPtr();
        TH2D* angularDist_data = angularDistPtr_data.GetPtr();
        baryonHist2_data->Write();
        baryonHist3_data->Write();
        angularDist_data->Write();
        massHist_data->Write();
    }
    for (int i = 0; i < mv_data.moments_histograms.size(); i++) {
        mv_data.moments_histograms[i]->Write();
        mv_data.moments_histograms_raw[i]->Write();
    }
    fout_data->Close();
    fout_data.reset();
    gROOT->cd();

    moments_df_accMC.ForeachSlot([mv_ptr_accMC](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double full_weight){mv_ptr_accMC->make_moments_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, full_weight);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","FullWeight"});
    mv_accMC.compile_thread_results();
    std::unique_ptr<TFile> fout_accMC(TFile::Open(TString::Format("%s/accMC/moments.root",analysis_info.output_path.Data()),"RECREATE"));
    fout_accMC->cd();
    if (!is_bootstrap) {
        TH1D* massHist_accMC = massHistPtr_accMC.GetPtr();
        TH1D* baryonHist2_accMC = mass2HistPtr_accMC.GetPtr();
        TH1D* baryonHist3_accMC = mass3HistPtr_accMC.GetPtr();
        TH2D* angularDist_accMC = angularDistPtr_accMC.GetPtr();
        baryonHist2_accMC->Write();
        baryonHist3_accMC->Write();
        angularDist_accMC->Write();
        massHist_accMC->Write();
    }
    for (int i = 0; i < mv_accMC.moments_histograms.size(); i++) {
        mv_accMC.moments_histograms[i]->Write();
        mv_accMC.moments_histograms_raw[i]->Write();
    }
    fout_accMC->Close();
    fout_accMC.reset();
    gROOT->cd();


    moments_df_genMC.ForeachSlot([mv_ptr_genMC](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double full_weight){mv_ptr_genMC->make_moments_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, full_weight);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","FullWeight"});
    mv_genMC.compile_thread_results();
    std::unique_ptr<TFile> fout_genMC(TFile::Open(TString::Format("%s/genMC/moments.root",analysis_info.output_path.Data()),"RECREATE"));
    fout_genMC->cd();
    TH1D* massHist_genMC = massHistPtr_genMC.GetPtr();
    if (!is_bootstrap) {
        TH1D* baryonHist2_genMC = mass2HistPtr_genMC.GetPtr();
        TH1D* baryonHist3_genMC = mass3HistPtr_genMC.GetPtr();
        TH2D* angularDist_genMC = angularDistPtr_genMC.GetPtr();
        baryonHist2_genMC->Write();
        baryonHist3_genMC->Write();
        angularDist_genMC->Write();
        massHist_genMC->Write();
    }
    massHist_genMC->SetDirectory(nullptr);
    for (int i = 0; i < mv_genMC.moments_histograms.size(); i++) {
        mv_genMC.moments_histograms[i]->Write();
        mv_genMC.moments_histograms_raw[i]->Write();
    }
    fout_genMC->Close();
    fout_genMC.reset();
    gROOT->cd();
    // ========================================================================================
    // Begin making the acceptance matrix
    // ========================================================================================

    
    AcceptanceMatrix acceptance_matrix(analysis_info.n_bins, analysis_info.meson_mass_min, analysis_info.meson_mass_max, maxL, nThreads, massHist_genMC, is_bootstrap);
    AcceptanceMatrix* acceptance_matrix_ptr = &acceptance_matrix;
    moments_df_accMC.ForeachSlot([acceptance_matrix_ptr](unsigned int thread_num, double PolarizationDegree, double helCosTheta, double helPhi, double PolarizationAngleReac, double MesonMass, double full_weight){acceptance_matrix_ptr->make_matrix_by_event(thread_num,PolarizationDegree, helCosTheta, helPhi,PolarizationAngleReac, MesonMass, full_weight);},{"PolarizationDegree","helCosTheta","helPhi","PolarizationAngleReac","MesonResonanceMass","FullWeight"});
    acceptance_matrix.makeMatrix();

    //open output file
    std::unique_ptr<TFile> acceptance_matrix_fout(TFile::Open(TString::Format("%s/accMC/acceptance_matrix.root",analysis_info.output_path.Data()),"RECREATE"));
    acceptance_matrix_fout->cd();
    std::vector<TMatrixD> inverse_acceptance_matrices;
    inverse_acceptance_matrices.reserve(acceptance_matrix.acceptance_matrices.size());
    for (int i = 0; i < acceptance_matrix.acceptance_matrices.size(); i++) {
        acceptance_matrix.acceptance_matrices[i].Write(TString::Format("acceptanceMatrix_%i",i));
        TH2D matrixHist(acceptance_matrix.acceptance_matrices[i]);
        matrixHist.SetDirectory(nullptr);
        matrixHist.Write(TString::Format("acceptanceMatrixHistogram_%i",i));
        TMatrixD inverse_matrix = acceptance_matrix.acceptance_matrices[i];
        inverse_matrix.Invert();
        inverse_matrix.Write(TString::Format("acceptanceMatrixInverse_%i",i));
        TH2D matrixHistInverse(inverse_matrix);
        matrixHistInverse.SetDirectory(nullptr);
        matrixHistInverse.Write(TString::Format("acceptanceMatrixHistogramInverse_%i",i));
        inverse_acceptance_matrices.push_back(std::move(inverse_matrix));
    }
    acceptance_matrix_fout->Close();
    acceptance_matrix_fout.reset();

    // ========================================================================================
    // Begin acceptance correcting
    // ========================================================================================

    MomentsArray moments_index_handler(maxL, 0,0,0,0);
    
    int dimension = 3*(maxL+1)*(maxL+2)/2 - maxL-1;
    std::vector<TMatrixD> moments_vectors_accMC;
    std::vector<TMatrixD> moments_vectors_data;
    std::vector<std::unique_ptr<TH1D>> output_histograms_raw_accMC;
    std::vector<std::unique_ptr<TH1D>> output_histograms_normalized_accMC;

    std::vector<std::unique_ptr<TH1D>> output_histograms_raw_data;
    std::vector<std::unique_ptr<TH1D>> output_histograms_normalized_data;

    ReactionSpecs reaction = analysis_info.reaction;

    std::cout << "Creating moments vectors" << std::endl;

    for (int bin=0; bin<analysis_info.n_bins; bin++) {
        moments_vectors_accMC.emplace_back(dimension,1);
        moments_vectors_data.emplace_back(dimension,1);
    }
    std::vector<int> pns = {2,3};

    std::cout << "creating empty moments histograms" << std::endl;
    for (int alpha=0; alpha < 3; alpha++) {
        for (int L=0; L <= maxL; L++) {
            for (int M=0; M <= L; M++) {
                if (alpha == 2 && M==0) {
                    continue;
                }
                TString unnormalized_title = TString::Format("H^{%i}(%i%i) for %s; %s (GeV/c^{2}) ; H^{%i}(%i%i)",alpha,L, M, analysis_info.reaction.getReaction().Data(), analysis_info.reaction.massString(pns).Data(), alpha,L, M);
                TString normalized_title = TString::Format("#LT_{}H^{%i}(%i%i)#GT for %s; %s (GeV/c^{2}) ; #LT_{}H^{%i}(%i%i)#GT",alpha, L, M, analysis_info.reaction.getReaction().Data(), analysis_info.reaction.massString(pns).Data(),alpha, L, M);
                TString unnormalized_fig_name = TString::Format("H_%i_%i_%i_raw_acceptance_corrected",alpha,L,M);
                TString normalized_fig_name = TString::Format("H_%i_%i_%i_acceptance_corrected",alpha,L,M);
                output_histograms_raw_data.push_back(std::make_unique<TH1D>(unnormalized_fig_name,unnormalized_fig_name,analysis_info.n_bins,analysis_info.meson_mass_min,analysis_info.meson_mass_max));
                output_histograms_raw_data.back()->SetTitle(unnormalized_title);
                output_histograms_raw_data.back()->SetDirectory(nullptr);
                output_histograms_raw_data.back()->SetStats(0);

                output_histograms_normalized_data.push_back(std::make_unique<TH1D>(normalized_fig_name,normalized_fig_name,analysis_info.n_bins,analysis_info.meson_mass_min,analysis_info.meson_mass_max));
                output_histograms_normalized_data.back()->SetTitle(normalized_title);
                output_histograms_normalized_data.back()->SetDirectory(nullptr);
                output_histograms_normalized_data.back()->SetStats(0);

                output_histograms_raw_accMC.push_back(std::make_unique<TH1D>(unnormalized_fig_name,unnormalized_fig_name,analysis_info.n_bins,analysis_info.meson_mass_min,analysis_info.meson_mass_max));
                output_histograms_raw_accMC.back()->SetTitle(unnormalized_title);
                output_histograms_raw_accMC.back()->SetDirectory(nullptr);
                output_histograms_raw_accMC.back()->SetStats(0);

                output_histograms_normalized_accMC.push_back(std::make_unique<TH1D>(normalized_fig_name,normalized_fig_name,analysis_info.n_bins,analysis_info.meson_mass_min,analysis_info.meson_mass_max));
                output_histograms_normalized_accMC.back()->SetTitle(normalized_title);
                output_histograms_normalized_accMC.back()->SetDirectory(nullptr);
                output_histograms_normalized_accMC.back()->SetStats(0);
            }
        }
    }
    std::unique_ptr<TLine> zero_line = std::make_unique<TLine>(analysis_info.meson_mass_min,0,analysis_info.meson_mass_max,0);
    std::cout << "loading moments vectors" << std::endl;
    for (int i=0; i < dimension; i++) {
        for (int bin=0; bin<analysis_info.n_bins; bin++) {
            moments_vectors_data[bin](i,0) = mv_data.moments_histograms_raw[i]->GetBinContent(bin+1);
            moments_vectors_accMC[bin](i,0) = mv_accMC.moments_histograms_raw[i]->GetBinContent(bin+1);
        }
    }
    std::cout << "acceptance correcting moments" << std::endl;
    for (int binNo=0; binNo<analysis_info.n_bins; binNo++) {
        TMatrixD truth_moments_data = inverse_acceptance_matrices[binNo] * moments_vectors_data[binNo];
        TMatrixD truth_moments_accMC = inverse_acceptance_matrices[binNo] * moments_vectors_accMC[binNo];
        for (int alpha = 0; alpha <3; alpha++) {
            for (int L=0; L <= maxL; L++) {
                for (int M=0; M <= L; M++) {
                    if (alpha == 2 && M == 0) {
                        continue;
                    }
                    int index = moments_index_handler.getMomentIndex(alpha, L, M);
                    double moment_value_data = truth_moments_data(index,0);
                    double moment_value_accMC = truth_moments_accMC(index,0);
                    output_histograms_raw_data[index]->SetBinContent(binNo+1, moment_value_data);
                    output_histograms_raw_accMC[index]->SetBinContent(binNo+1, moment_value_accMC);
            }
            }
        }
    }
    std::cout << "saving acceptance corrected moments" << std::endl;
    fout_data.reset(TFile::Open(TString::Format("%s/data/moments.root",analysis_info.output_path.Data()),"UPDATE"));
    fout_data->cd();
    for (int alpha=0; alpha < 3; alpha++) {
        for (int L=0; L <= maxL; L++) {
            for (int M=0; M <= L; M++) {
                if (alpha == 2 && M == 0) {
                    continue;
                }
                int index = moments_index_handler.getMomentIndex(alpha, L, M);
                output_histograms_raw_data[index]->SetMinimum(1.1*min(0.,output_histograms_raw_data[index]->GetMinimum()));
                output_histograms_raw_data[index]->SetMaximum(1.1*max(0.,output_histograms_raw_data[index]->GetMaximum()));
                output_histograms_raw_data[index]->Write(TString::Format("H_%i_%i_%i_raw_acceptance_corrected",alpha, L,M));
                output_histograms_normalized_data[index]->Divide(output_histograms_raw_data[index].get(),output_histograms_raw_data[0].get());
                output_histograms_normalized_data[index]->SetMaximum(1.1*max(0.,output_histograms_normalized_data[index]->GetMaximum()));
                output_histograms_normalized_data[index]->SetMinimum(1.1*min(0.,output_histograms_normalized_data[index]->GetMinimum()));
                output_histograms_normalized_data[index]->Write(TString::Format("H_%i_%i_%i_acceptance_corrected",alpha, L,M));
            }
        }
    }
    fout_data->Close();
    fout_data.reset();

    fout_accMC.reset(TFile::Open(TString::Format("%s/accMC/moments.root",analysis_info.output_path.Data()),"UPDATE"));
    fout_accMC->cd();
    for (int alpha=0; alpha < 3; alpha++) {
        for (int L=0; L <= maxL; L++) {
            for (int M=0; M <= L; M++) {
                if (alpha == 2 && M == 0) {
                    continue;
                }
                int index = moments_index_handler.getMomentIndex(alpha, L, M);
                output_histograms_raw_accMC[index]->SetMinimum(1.1*min(0.,output_histograms_raw_accMC[index]->GetMinimum()));
                output_histograms_raw_accMC[index]->SetMaximum(1.1*max(0.,output_histograms_raw_accMC[index]->GetMaximum()));
                output_histograms_raw_accMC[index]->Write(TString::Format("H_%i_%i_%i_raw_acceptance_corrected",alpha, L,M));
                output_histograms_normalized_accMC[index]->Divide(output_histograms_raw_accMC[index].get(),output_histograms_raw_accMC[0].get());
                output_histograms_normalized_accMC[index]->SetMaximum(1.1*max(0.,output_histograms_normalized_accMC[index]->GetMaximum()));
                output_histograms_normalized_accMC[index]->SetMinimum(1.1*min(0.,output_histograms_normalized_accMC[index]->GetMinimum()));
                output_histograms_normalized_accMC[index]->Write(TString::Format("H_%i_%i_%i_acceptance_corrected",alpha, L,M));
            }
        }
    }
    fout_accMC->Close();
    fout_accMC.reset();
}
