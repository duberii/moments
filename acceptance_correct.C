#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"
#include <TMatrixD.h>

void acceptance_correct(const char* dataset_name, const char* data_type, int maxL, int bootstrap_number) {
    bool is_bootstrap = (bootstrap_number > -1);

    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    if (is_bootstrap) {
        analysis_info.output_path=TString::Format("%s%s/%i/bootstrap_%i/",analysis_info.output_path.Data(),dataset_name, maxL, bootstrap_number);
    } else {
        analysis_info.output_path=TString::Format("%s%s/%i/full/",analysis_info.output_path.Data(),dataset_name,maxL);
    }
    
    MomentsArray moments_index_handler(maxL, 0,0,0,0);
    
    int dimension = 3*(maxL+1)*(maxL+2)/2 - maxL-1;
    std::unique_ptr<TFile> acceptance_matrix_file = std::unique_ptr<TFile>(TFile::Open(analysis_info.output_path + "accMC/acceptance_matrix.root","READ"));
    std::unique_ptr<TFile> uncorrected_moments_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/moments.root", analysis_info.output_path.Data(),data_type),"UPDATE"));
    std::vector<TMatrixD> moments_vectors;
    std::vector<TH1D*> moments_histograms;
    std::vector<TMatrixD> inverse_acceptance_matrices;
    std::vector<std::unique_ptr<TH1D>> output_histograms_raw;
    std::vector<std::unique_ptr<TH1D>> output_histograms_normalized;


    ReactionSpecs reaction = analysis_info.reaction;

    std::cout << "loading matrices" << std::endl;

    for (int bin=0; bin<analysis_info.n_bins; bin++) {
        moments_vectors.emplace_back(dimension,1);
        TMatrixD* inverse_acceptance_matrix = (TMatrixD*)acceptance_matrix_file->Get(TString::Format("acceptanceMatrixInverse_%i",bin));
        inverse_acceptance_matrices.push_back(*inverse_acceptance_matrix);
    }
    double xMin = 0;
    double xMax = 0;
    std::vector<int> pns = {2,3};

    std::cout << "loading moments and creating histograms" << std::endl;
    for (int alpha=0; alpha < 3; alpha++) {
        for (int L=0; L <= maxL; L++) {
            for (int M=0; M <= L; M++) {
                if (alpha == 2 && M==0) {
                    continue;
                }
                moments_histograms.push_back((TH1D*)uncorrected_moments_file->Get(TString::Format("H_%i_%i_%i_raw",alpha,L,M)));
                if (xMin==xMax) {
                    xMin = moments_histograms[0]->GetXaxis()->GetXmin();
                    xMax = moments_histograms[0]->GetXaxis()->GetXmax();
                }
                TString unnormalized_title = TString::Format("H^{%i}(%i%i) for %s; %s (GeV/c^{2}) ; H^{%i}(%i%i)",alpha,L, M, reaction.getReaction().Data(), reaction.massString(pns).Data(), alpha,L, M);
                TString normalized_title = TString::Format("#LT_{}H^{%i}(%i%i)#GT for %s; %s (GeV/c^{2}) ; #LT_{}H^{%i}(%i%i)#GT",alpha, L, M, reaction.getReaction().Data(), reaction.massString(pns).Data(),alpha, L, M);
                TString unnormalized_fig_name = TString::Format("H_%i_%i_%i_raw_acceptance_corrected",alpha,L,M);
                TString normalized_fig_name = TString::Format("H_%i_%i_%i_acceptance_corrected",alpha,L,M);
                output_histograms_raw.push_back(std::make_unique<TH1D>(unnormalized_fig_name,unnormalized_fig_name,analysis_info.n_bins,xMin,xMax));
                output_histograms_raw.back()->SetTitle(unnormalized_title);
                output_histograms_raw.back()->SetDirectory(nullptr);
                output_histograms_raw.back()->SetStats(0);

                output_histograms_normalized.push_back(std::make_unique<TH1D>(normalized_fig_name,normalized_fig_name,analysis_info.n_bins,xMin,xMax));
                output_histograms_normalized.back()->SetTitle(normalized_title);
                output_histograms_normalized.back()->SetDirectory(nullptr);
                output_histograms_normalized.back()->SetStats(0);
            }
        }
    }
    std::unique_ptr<TLine> zero_line = std::make_unique<TLine>(xMin,0,xMax,0);
    std::cout << "loading moments vectors" << std::endl;
    for (int i=0; i < dimension; i++) {
        for (int bin=0; bin<analysis_info.n_bins; bin++) {
            moments_vectors[bin](i,0) = moments_histograms[i]->GetBinContent(bin+1);
        }
    }
    std::cout << "acceptance correcting moments" << std::endl;
    for (int binNo=0; binNo<analysis_info.n_bins; binNo++) {
        TMatrixD truth_moments = inverse_acceptance_matrices[binNo] * moments_vectors[binNo];
        for (int alpha = 0; alpha <3; alpha++) {
            for (int L=0; L <= maxL; L++) {
                for (int M=0; M <= L; M++) {
                    if (alpha == 2 && M == 0) {
                        continue;
                    }
                    int index = moments_index_handler.getMomentIndex(alpha, L, M);
                    double moment_value = truth_moments(index,0);
                    output_histograms_raw[index]->SetBinContent(binNo+1, moment_value);
            }
            }
        }
    }
    std::cout << "saving acceptance corrected moments" << std::endl;
    uncorrected_moments_file->cd();
    //TCanvas canvas("canvas","canvas",700,500);
    for (int alpha=0; alpha < 3; alpha++) {
        for (int L=0; L <= maxL; L++) {
            for (int M=0; M <= L; M++) {
                if (alpha == 2 && M == 0) {
                    continue;
                }
                int index = moments_index_handler.getMomentIndex(alpha, L, M);
                output_histograms_raw[index]->SetMinimum(1.1*min(0.,output_histograms_raw[index]->GetMinimum()));
                output_histograms_raw[index]->SetMaximum(1.1*max(0.,output_histograms_raw[index]->GetMaximum()));
                output_histograms_raw[index]->Write(TString::Format("H_%i_%i_%i_raw_acceptance_corrected",alpha, L,M));
                output_histograms_normalized[index]->Divide(output_histograms_raw[index].get(),output_histograms_raw[0].get());
                output_histograms_normalized[index]->SetMaximum(1.1*max(0.,output_histograms_normalized[index]->GetMaximum()));
                output_histograms_normalized[index]->SetMinimum(1.1*min(0.,output_histograms_normalized[index]->GetMinimum()));
                output_histograms_normalized[index]->Write(TString::Format("H_%i_%i_%i_acceptance_corrected",alpha, L,M));
            }
        }
    }
    uncorrected_moments_file->Close();
    acceptance_matrix_file->Close();
}
