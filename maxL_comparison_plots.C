#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

void maxL_comparison_plots(const char* dataset_name, std::vector<int> maxLs_to_compare) {
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    std::unique_ptr<TFile> output_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/moments.root", analysis_info.output_path.Data(), dataset_name),"UPDATE"));
    std::vector<std::unique_ptr<TFile>> input_files;
    std::vector<TString> labels;
    int largest_L = 0;
    for (int maxL : maxLs_to_compare) {
        input_files.emplace_back(TFile::Open(TString::Format("%s%s/%i/full/data/moments.root", analysis_info.output_path.Data(), dataset_name, maxL),"UPDATE"));
        labels.push_back(TString::Format("L_{max}=%i",maxL));
        if (largest_L < maxL) {
            largest_L = maxL;
        }
    }
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("temp_canv","temp_canv",700,500);
    //compare acceptance corrected moments between maxLs
    TString output_path = TString::Format("%s%s/",analysis_info.output_path.Data(),dataset_name);
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L =0; L <= largest_L; L++) {
            for (int M =0; M <= L; M++) { 
                if (alpha==2 && M==0) {
                    continue;
                }
                std::vector<TH1D*> histograms;
                std::vector<TString> labels_to_use;
                TString fig_name_in;
                TString fig_name_out;
                if (alpha == 0 && L == 0 && M == 0) {
                    fig_name_in = TString::Format("H_%i_%i_%i_raw_acceptance_corrected", alpha, L, M);
                    fig_name_out = TString::Format("H_%i_%i_%i_raw", alpha, L, M);
                } else {
                    fig_name_in = TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M);
                    fig_name_out = TString::Format("H_%i_%i_%i", alpha, L, M);
                }
                for (unsigned int i=0; i < maxLs_to_compare.size(); i++) {
                    if (L <= maxLs_to_compare[i]) {
                        TH1D* hist = (TH1D*)input_files[i]->Get(fig_name_in);
                        hist->SetDirectory(nullptr);
                        histograms.push_back(hist);
                        labels_to_use.push_back(labels[i]);
                    }                    
                }
                draw_on_same_canvas(histograms,labels_to_use,fig_name_out,output_file.get(),canvas.get(), output_path);                
            }
        }
    }
    for (unsigned int i= 0; i < input_files.size(); i++) {
        input_files[i]->Close();
    }
    output_file->Close();
}
