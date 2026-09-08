#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"
#include <TMatrixD.h>


void aggregate_results(const char* dataset_name, int maxL, int num_bootstraps) {
    
    
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    std::vector<TString> figure_names;
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas","canvas",700,500);
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L =0; L <= maxL; L++) {
            for (int M =0; M <= L; M++) { 
                if (alpha==2 && M==0) {
                    continue;
                }
                figure_names.push_back(TString::Format("H_%i_%i_%i", alpha, L, M));
                figure_names.push_back(TString::Format("H_%i_%i_%i_raw", alpha, L, M));
                figure_names.push_back(TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M));
                figure_names.push_back(TString::Format("H_%i_%i_%i_raw_acceptance_corrected", alpha, L, M));
            }
        }
    }
    std::vector<TString> data_files =  {"data","accMC","genMC"};
    for (int data_file_number = 0; data_file_number < int(data_files.size()); data_file_number++) {
        std::unique_ptr<TFile> main_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/%s/moments.root", analysis_info.output_path.Data(), dataset_name, maxL, data_files[data_file_number].Data()),"UPDATE"));
        std::vector<std::vector<double>> average(analysis_info.n_bins, std::vector<double>(figure_names.size()));
        std::vector<std::vector<double>> square_average(analysis_info.n_bins, std::vector<double>(figure_names.size()));
        std::vector<std::vector<double>> std_devs(analysis_info.n_bins, std::vector<double>(figure_names.size()));
        std::vector<bool> figure_names_in_current_file(figure_names.size());
        for (int figure_number = 0; figure_number < int(figure_names.size()); figure_number++) { 
            TH1D* figure = dynamic_cast<TH1D*>(main_file->Get(figure_names[figure_number]));
            figure_names_in_current_file[figure_number] = (figure != nullptr);
        }
        for (int bootstrap_number = 0; bootstrap_number < num_bootstraps; bootstrap_number++) {
            std::unique_ptr<TFile> bootstrap_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/bootstrap_%i/%s/moments.root", analysis_info.output_path.Data(), dataset_name, maxL, bootstrap_number, data_files[data_file_number].Data()),"READ"));
            for (int figure_number = 0; figure_number < int(figure_names_in_current_file.size()); figure_number++) {
                if (!figure_names_in_current_file[figure_number]) {
                    continue;
                }
                TH1D* bootstrap_figure = dynamic_cast<TH1D*>(bootstrap_file->Get(figure_names[figure_number]));
                for (int bin=0; bin < analysis_info.n_bins; bin++) {
                    double bin_val = bootstrap_figure->GetBinContent(bin+1);
                    average[bin][figure_number] += bin_val;
                    square_average[bin][figure_number] += bin_val*bin_val;
                }
            }
            bootstrap_file->Close();
        }
        
        for (int figure_number = 0; figure_number < int(figure_names_in_current_file.size()); figure_number++) {
            if (!figure_names_in_current_file[figure_number]) {
                continue;
            }
            TH1D* figure = dynamic_cast<TH1D*>(main_file->Get(figure_names[figure_number]));
            double y_min = 0.0;
            double y_max = 0.0;
            for (int bin=0; bin < analysis_info.n_bins; bin++) {
                average[bin][figure_number] /= num_bootstraps;
                square_average[bin][figure_number] /= num_bootstraps;
                double diff = square_average[bin][figure_number] - average[bin][figure_number]*average[bin][figure_number];
                if (diff <= 0) {
                    std_devs[bin][figure_number] = 0.;
                } else {
                    std_devs[bin][figure_number] = sqrt(diff);
                }
                figure->SetBinError(bin+1, std_devs[bin][figure_number]);
                y_min = min(figure->GetBinContent(bin+1)-std_devs[bin][figure_number], y_min);
                y_max = max(figure->GetBinContent(bin+1)+std_devs[bin][figure_number], y_max);
            }
            figure->SetMaximum(1.1*y_max);
            figure->SetMinimum(1.1*y_min);
            main_file->cd();
            figure->Write();
            //canvas->cd();
            //std::unique_ptr<TLine> zero_line = std::make_unique<TLine>(figure->GetXaxis()->GetXmin(),0,figure->GetXaxis()->GetXmax(),0);
            //figure->Draw("E1");
            //zero_line->Draw("SAME");
            //canvas->SaveAs(TString::Format("%sall/moments/%s/%s.png", analysis_info.output_path, data_files[data_file_number].Data(),figure_names[figure_number].Data()));
        }
        main_file->Close();
    }

}