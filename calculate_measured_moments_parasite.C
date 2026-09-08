#include <stdexcept>
#include <TString.h>
#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

void calculate_measured_moments_parasite(const char* dataset_name, const char* data_type, int maxL, int bootstrap_number, int largerMaxL) {
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    bool is_bootstrap = (bootstrap_number>-1);

    // Sets the number of threads based on the SLURM environment variable
    int nThreads = std::stoi(std::getenv("SLURM_CPUS_PER_TASK"));
    ROOT::EnableImplicitMT(nThreads);

    //Get analysis-specific information
    std::unique_ptr<TFile> fout;
    std::unique_ptr<TFile> fin;
    if (is_bootstrap) {
        fout = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/bootstrap_%i/%s/moments.root",analysis_info.output_path.Data(), dataset_name, maxL, bootstrap_number, data_type),"UPDATE"));
        fin = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/bootstrap_%i/%s/moments.root",analysis_info.output_path.Data(), dataset_name, largerMaxL, bootstrap_number, data_type),"READ"));
    } else {
        fout = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/%s/moments.root",analysis_info.output_path.Data(), dataset_name, maxL, data_type),"UPDATE"));
        fin = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/%s/moments.root",analysis_info.output_path.Data(), dataset_name, largerMaxL, data_type),"READ"));
    }
    std::cout << "done opening files" << std::endl;

    for (int alpha = 0; alpha < 3; alpha++) {
        for (int L=0; L <= maxL; L++) {
            for (int M = 0; M <= L; M++) {
                if (alpha == 2 && M==0) {
                    continue;
                }
                TString raw_file_name = TString::Format("H_%i_%i_%i_raw",alpha, L, M);
                TString file_name = TString::Format("H_%i_%i_%i",alpha, L, M);
                TH1D* raw_moment = (TH1D*)fin->Get(raw_file_name);
                TH1D* moment = (TH1D*)fin->Get(file_name);
                
                fout->cd();
                moment->Write(file_name);
                raw_moment->Write(raw_file_name);
                moment->SetDirectory(nullptr);
                raw_moment->SetDirectory(nullptr);
            }
        }
    }
    std::cout << "finished saving" << std::endl;
    fin->Close();
    fout->Close();
}
