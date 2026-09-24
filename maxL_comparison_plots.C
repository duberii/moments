#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

void maxL_comparison_plots(const char* dataset_name, std::vector<int> maxLs_to_compare) {
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    std::unique_ptr<TFile> output_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/moments.root", analysis_info.output_path.Data(), dataset_name),"UPDATE"));
    std::vector<std::unique_ptr<TFile>> input_files;
    std::vector<TString> labels;
    int largest_L = 0;
    std::cout << "Opening data files" << std::endl;
    for (int maxL : maxLs_to_compare) {
        input_files.emplace_back(TFile::Open(TString::Format("%s%s/%i/full/data/moments.root", analysis_info.output_path.Data(), dataset_name, maxL),"UPDATE"));
        labels.push_back(TString::Format("L_{max}=%i",maxL));
        if (largest_L < maxL) {
            largest_L = maxL;
        }
    }
    std::cout << "Done opening data files" << std::endl;
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("temp_canv","temp_canv",700,500);
    //compare acceptance corrected moments between maxLs
    TString output_path = TString::Format("%s%s/maxL_comparison_plots/",analysis_info.output_path.Data(),dataset_name);
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L =0; L <= largest_L; L++) {
            for (int M =0; M <= L; M++) { 
                if (alpha==2 && M==0) {
                    continue;
                }
                std::vector<TH1D*> histograms;
                std::vector<TH1D*> histograms_raw;
                std::vector<TH1D*> histograms_statistical;
                std::vector<TH1D*> histograms_raw_statistical;
                std::vector<TString> labels_to_use;
                TString fig_name_in_raw = TString::Format("H_%i_%i_%i_raw_acceptance_corrected", alpha, L, M);
                TString fig_name_out_raw = TString::Format("H_%i_%i_%i_raw", alpha, L, M);
                TString fig_name_in = TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M);
                TString fig_name_out = TString::Format("H_%i_%i_%i", alpha, L, M);

                //Only raw
                std::cout << "Making alpha=" << alpha << " L=" << L << " M=" << M << " raw histograms"<< std::endl;
                for (unsigned int i=0; i < maxLs_to_compare.size(); i++) {
                    if (L <= maxLs_to_compare[i]) {
                        TH1D* hist_raw = (TH1D*)input_files[i]->Get(fig_name_in_raw);
                        TH1D* hist_raw_statistical = (TH1D*)input_files[i]->Get(fig_name_in_raw + "_statistical");
                        hist_raw->SetDirectory(nullptr);
                        hist_raw_statistical->SetDirectory(nullptr);
                        histograms_raw.push_back(hist_raw);
                        histograms_raw_statistical.push_back(hist_raw_statistical);
                        labels_to_use.push_back(labels[i]);
                    }
                }
                draw_on_same_canvas_with_statistical_errors(histograms_raw, histograms_raw_statistical,labels_to_use,fig_name_out_raw, output_file.get(), canvas.get(), output_path); 
                std::cout << "Done making alpha=" << alpha << " L=" << L << " M=" << M << " raw histograms"<< std::endl;
                if (alpha == 0 && L == 0 && M == 0) {
                    continue;
                }

                std::cout << "Making alpha=" << alpha << " L=" << L << " M=" << M << " histograms"<< std::endl;

                //Only normalized
                for (unsigned int i=0; i < maxLs_to_compare.size(); i++) {
                    if (L <= maxLs_to_compare[i]) {
                        TH1D* hist = (TH1D*)input_files[i]->Get(fig_name_in);
                        TH1D* hist_statistical = (TH1D*)input_files[i]->Get(fig_name_in + "_statistical");
                        hist->SetDirectory(nullptr);
                        hist_statistical->SetDirectory(nullptr);
                        histograms.push_back(hist);
                        histograms_statistical.push_back(hist_statistical);
                        labels_to_use.push_back(labels[i]);
                    }
                }
                draw_on_same_canvas_with_statistical_errors(histograms, histograms_statistical,labels_to_use,fig_name_out, output_file.get(), canvas.get(), output_path); 
                        
                std::cout << "Done making alpha=" << alpha << " L=" << L << " M=" << M << " histograms"<< std::endl;
            }
        }
    }
    for (unsigned int i= 0; i < input_files.size(); i++) {
        input_files[i]->Close();
    }
    output_file->Close();
}
