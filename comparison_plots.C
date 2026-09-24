#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

void comparison_plots(const char* dataset_name, int maxL) {
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    std::unique_ptr<TFile> data_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/data/moments.root", analysis_info.output_path.Data(), dataset_name, maxL),"READ"));
    std::unique_ptr<TFile> accMC_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/accMC/moments.root", analysis_info.output_path.Data(), dataset_name, maxL),"READ"));
    std::unique_ptr<TFile> genMC_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/genMC/moments.root", analysis_info.output_path.Data(), dataset_name, maxL),"READ"));
    std::unique_ptr<TFile> output_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/moments.root", analysis_info.output_path.Data(), dataset_name, maxL),"UPDATE"));
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("temp_canv","temp_canv",700,500);
    //compare acceptance corrected accMC to genMC (only H^0(0)))
    {
        TH1D* accMC_hist = (TH1D*)accMC_file->Get(TString("H_0_0_0_raw_acceptance_corrected"));
        TH1D* genMC_hist = (TH1D*)genMC_file->Get(TString("H_0_0_0_raw"));
        draw_on_same_canvas({genMC_hist, accMC_hist}, {"genMC", "acceptance-corrected accMC"}, TString("H_0_0_0_accTest"), output_file.get(), canvas.get(), TString::Format("%s%s/%i/acceptance_correcting_comparison_plots/", analysis_info.output_path.Data(), dataset_name, maxL));
    }
    //compare uncorrected moments to corrected ones
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L =0; L <= maxL; L++) {
            for (int M =0; M <= L; M++) { 
                if (alpha==2 && M==0) {
                    continue;
                }
                TH1D* uncorrected_hist_raw = (TH1D*)data_file->Get(TString::Format("H_%i_%i_%i_raw", alpha, L, M));
                TH1D* corrected_hist_raw = (TH1D*)data_file->Get(TString::Format("H_%i_%i_%i_raw_acceptance_corrected", alpha, L, M));
                TH1D* corrected_hist_raw_statistical = (TH1D*)data_file->Get(TString::Format("H_%i_%i_%i_raw_acceptance_corrected_statistical", alpha, L, M));
                draw_on_same_canvas_with_statistical_errors({uncorrected_hist_raw,  corrected_hist_raw},{uncorrected_hist_raw,  corrected_hist_raw_statistical}, {"uncorrected moment","acceptance-corrected moment"}, TString::Format("H_%i_%i_%i_raw", alpha, L, M), output_file.get(), canvas.get(), TString::Format("%s%s/%i/acceptance_correcting_comparison_plots/", analysis_info.output_path.Data(), dataset_name, maxL));
                
                if (alpha != 0 || L != 0 || M != 0) {
                    TH1D* uncorrected_hist = (TH1D*)data_file->Get(TString::Format("H_%i_%i_%i", alpha, L, M));
                    TH1D* corrected_hist = (TH1D*)data_file->Get(TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M));
                    TH1D* corrected_hist_statistical = (TH1D*)data_file->Get(TString::Format("H_%i_%i_%i_acceptance_corrected_statistical", alpha, L, M));
                    draw_on_same_canvas_with_statistical_errors({uncorrected_hist, corrected_hist},{uncorrected_hist, corrected_hist_statistical}, {"uncorrected moment",  "acceptance-corrected moment"}, TString::Format("H_%i_%i_%i", alpha, L, M), output_file.get(), canvas.get(), TString::Format("%s%s/%i/acceptance_correcting_comparison_plots/", analysis_info.output_path.Data(), dataset_name, maxL));
                } 
            }
        }
    }
    data_file->Close();
    accMC_file->Close();
    genMC_file->Close();
    output_file->Close();
}
