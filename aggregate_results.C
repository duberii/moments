#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"
#include <TMatrixD.h>


void aggregate_results(const char* dataset_name, int maxL, int num_bootstraps) {
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    std::vector<TString> figure_names;
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas","canvas",700,500);
    int dimension = 3*(maxL+1)*(maxL+2)/2 - maxL-1;
    MomentsArray moments_index_handler(maxL);

    //===================================================
    // Making fully bootstrapped figures
    //===================================================

    // First, read the nominal values from the "full" file

    std::unique_ptr<TFile> main_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/data/moments.root", analysis_info.output_path.Data(), dataset_name, maxL),"UPDATE"));
    
    std::vector<std::vector<double>> means = std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::vector<std::vector<double>> means_raw = std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::vector<std::vector<double>> means_acceptance_corrected = std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::vector<std::vector<double>> means_raw_acceptance_corrected = std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::cout<< "Beginning to extract means from main file" << std::endl;
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L =0; L <= maxL; L++) {
            for (int M =0; M <= L; M++) { 
                if (alpha==2 && M==0) {
                    continue;
                }
                int moment_number = moments_index_handler.getMomentIndex(alpha,L,M);
                TH1D* moments_histogram = (TH1D*)main_file->Get(TString::Format("H_%i_%i_%i", alpha, L, M));
                Double_t* moments_values = moments_histogram->GetArray();
                TH1D* moments_histogram_raw = (TH1D*)main_file->Get(TString::Format("H_%i_%i_%i_raw", alpha, L, M));
                Double_t* moments_values_raw = moments_histogram_raw->GetArray();
                TH1D* moments_histogram_acceptance_corrected = (TH1D*)main_file->Get(TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M));
                Double_t* moments_values_acceptance_corrected = moments_histogram_acceptance_corrected->GetArray();
                TH1D* moments_histogram_raw_acceptance_corrected = (TH1D*)main_file->Get(TString::Format("H_%i_%i_%i_raw_acceptance_corrected", alpha, L, M));
                Double_t* moments_values_raw_acceptance_corrected = moments_histogram_raw_acceptance_corrected->GetArray();
                for (int bin = 1; bin <= analysis_info.n_bins; bin++) {
                    means[moment_number][bin-1] = moments_values[bin];
                    means_raw[moment_number][bin-1] = moments_values_raw[bin];
                    means_acceptance_corrected[moment_number][bin-1] = moments_values_acceptance_corrected[bin];
                    means_raw_acceptance_corrected[moment_number][bin-1] = moments_values_raw_acceptance_corrected[bin];
                }
            }
        }
    }
    std::cout<< "Done extracting means from main file" << std::endl;

    // Now find the std dev for all the bootstrapped figures
    std::vector<std::vector<double>> std_devs = std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::vector<std::vector<double>> std_devs_raw = std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::vector<std::vector<double>> std_devs_acceptance_corrected = std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::vector<std::vector<double>> std_devs_raw_acceptance_corrected = std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    
    std::cout<< "Beginning to calculate std dev from bootstrapped files" << std::endl;
    for (int bootstrap_number = 0; bootstrap_number < num_bootstraps; bootstrap_number++) {
        std::unique_ptr<TFile> bootstrap_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/bootstraps/%i/data/moments.root", analysis_info.output_path.Data(), dataset_name, maxL, bootstrap_number),"READ"));
        for (int alpha = 0; alpha <3; alpha++) {
            for (int L =0; L <= maxL; L++) {
                for (int M =0; M <= L; M++) { 
                    if (alpha==2 && M==0) {
                        continue;
                    }
                    int moment_number = moments_index_handler.getMomentIndex(alpha,L,M);
                    TH1D* moments_histogram = (TH1D*)bootstrap_file->Get(TString::Format("H_%i_%i_%i", alpha, L, M));
                    Double_t* moments_values = moments_histogram->GetArray();
                    TH1D* moments_histogram_raw = (TH1D*)bootstrap_file->Get(TString::Format("H_%i_%i_%i_raw", alpha, L, M));
                    Double_t* moments_values_raw = moments_histogram_raw->GetArray();
                    TH1D* moments_histogram_acceptance_corrected = (TH1D*)bootstrap_file->Get(TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M));
                    Double_t* moments_values_acceptance_corrected = moments_histogram_acceptance_corrected->GetArray();
                    TH1D* moments_histogram_raw_acceptance_corrected = (TH1D*)bootstrap_file->Get(TString::Format("H_%i_%i_%i_raw_acceptance_corrected", alpha, L, M));
                    Double_t* moments_values_raw_acceptance_corrected = moments_histogram_raw_acceptance_corrected->GetArray();
                    for (int bin = 1; bin <= analysis_info.n_bins; bin++) {
                        std_devs[moment_number][bin-1] += pow(moments_values[bin]-means[moment_number][bin-1],2);
                        std_devs_raw[moment_number][bin-1] += pow(moments_values_raw[bin]-means_raw[moment_number][bin-1],2);
                        std_devs_acceptance_corrected[moment_number][bin-1] += pow(moments_values_acceptance_corrected[bin]-means_acceptance_corrected[moment_number][bin-1],2);
                        std_devs_raw_acceptance_corrected[moment_number][bin-1] += pow(moments_values_raw_acceptance_corrected[bin]-means_raw_acceptance_corrected[moment_number][bin-1],2);
                    }
                }
            }
        }
    }
    std::cout<< "Finished calculating std dev from bootstrapped files" << std::endl;

    //Actually make and save the new figures with the right uncertainties

    main_file->cd();
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L =0; L <= maxL; L++) {
            for (int M =0; M <= L; M++) { 
                if (alpha==2 && M==0) {
                    continue;
                }
                int moment_number = moments_index_handler.getMomentIndex(alpha,L,M);
                TH1D* moments_histogram = (TH1D*)main_file->Get(TString::Format("H_%i_%i_%i", alpha, L, M));
                TH1D* moments_histogram_raw = (TH1D*)main_file->Get(TString::Format("H_%i_%i_%i_raw", alpha, L, M));
                TH1D* moments_histogram_acceptance_corrected = (TH1D*)main_file->Get(TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M));
                TH1D* moments_histogram_raw_acceptance_corrected = (TH1D*)main_file->Get(TString::Format("H_%i_%i_%i_raw_acceptance_corrected", alpha, L, M));
                for (int bin = 1; bin <= analysis_info.n_bins; bin++) {
                    moments_histogram->SetBinError(bin,sqrt(std_devs[moment_number][bin-1]/num_bootstraps));
                    moments_histogram_raw->SetBinError(bin,sqrt(std_devs_raw[moment_number][bin-1]/num_bootstraps));
                    moments_histogram_acceptance_corrected->SetBinError(bin,sqrt(std_devs_acceptance_corrected[moment_number][bin-1]/num_bootstraps));
                    moments_histogram_raw_acceptance_corrected->SetBinError(bin,sqrt(std_devs_raw_acceptance_corrected[moment_number][bin-1]/num_bootstraps));
                }
                moments_histogram->Write();
                moments_histogram_raw->Write();
                moments_histogram_acceptance_corrected->Write();
                moments_histogram_raw_acceptance_corrected->Write();
            }
        }
    }

    //===================================================
    // Making statistical uncertainty figures
    //===================================================

    // Now find the std dev for all the bootstrapped figures
    std::vector<std::vector<double>> std_devs_acceptance_corrected_statistical= std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::vector<std::vector<double>> std_devs_raw_acceptance_corrected_statistical= std::vector<std::vector<double>>(dimension,std::vector<double>(analysis_info.n_bins));
    std::cout<< "Beginning to calculate std dev (statstical error) from bootstrapped files" << std::endl;
    std::vector<TMatrixD> inverse_acceptance_matrices;
    inverse_acceptance_matrices.reserve(analysis_info.n_bins);
    TFile* acceptance_matrix_fin = TFile::Open(TString::Format("%s%s/%i/full/accMC/acceptance_matrix.root", analysis_info.output_path.Data(), dataset_name, maxL),"READ");
    for (int bin = 0; bin < analysis_info.n_bins; bin++) {
        TMatrixD* mtx = (TMatrixD*)acceptance_matrix_fin->Get(TString::Format("acceptanceMatrixInverse_%i",bin));
        inverse_acceptance_matrices.emplace_back(*mtx);
    }
    for (int bootstrap_number = 0; bootstrap_number < num_bootstraps; bootstrap_number++) {
        std::unique_ptr<TFile> bootstrap_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/bootstraps/%i/data/moments.root", analysis_info.output_path.Data(), dataset_name, maxL, bootstrap_number),"READ"));
        std::vector<TMatrixD> moments_vectors_raw;
        for (int bin = 0; bin < analysis_info.n_bins; bin++) {
            moments_vectors_raw.emplace_back(dimension,1);
        }
        for (int alpha = 0; alpha <3; alpha++) {
            for (int L =0; L <= maxL; L++) {
                for (int M =0; M <= L; M++) { 
                    if (alpha==2 && M==0) {
                        continue;
                    }
                    int moment_number = moments_index_handler.getMomentIndex(alpha,L,M);
                    TH1D* moments_histogram_raw = (TH1D*)bootstrap_file->Get(TString::Format("H_%i_%i_%i_raw", alpha, L, M));
                    Double_t* moments_values_raw = moments_histogram_raw->GetArray();
                    TH1D* moments_histogram = (TH1D*)bootstrap_file->Get(TString::Format("H_%i_%i_%i", alpha, L, M));
                    Double_t* moments_values = moments_histogram->GetArray();
                    for (int bin = 1; bin <= analysis_info.n_bins; bin++) {
                        moments_vectors_raw[bin-1](moment_number,0) = moments_values_raw[bin];
                    }
                }
            }
        }
        for (int bin = 1; bin <= analysis_info.n_bins; bin++) {
            TMatrixD truth_moments_vector = inverse_acceptance_matrices[bin-1] * moments_vectors_raw[bin-1];
            for (int moment_number = 0; moment_number < dimension; moment_number++) {
                std_devs_raw_acceptance_corrected_statistical[moment_number][bin-1] += pow(truth_moments_vector(moment_number,0)-means_raw_acceptance_corrected[moment_number][bin-1],2);
                std_devs_acceptance_corrected_statistical[moment_number][bin-1] += pow(truth_moments_vector(moment_number,0)/truth_moments_vector(0,0)-means_acceptance_corrected[moment_number][bin-1],2);
            }
        }
    }

    //Actually make and save the new figures with the right (statistical) uncertainties

    main_file->cd();
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L =0; L <= maxL; L++) {
            for (int M =0; M <= L; M++) { 
                if (alpha==2 && M==0) {
                    continue;
                }
                int moment_number = moments_index_handler.getMomentIndex(alpha,L,M);
                TH1D* moments_histogram_acceptance_corrected = (TH1D*)((TH1D*)main_file->Get(TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M)))->Clone(TString::Format("H_%i_%i_%i_acceptance_corrected_statistical", alpha, L, M));
                TH1D* moments_histogram_raw_acceptance_corrected = (TH1D*)((TH1D*)main_file->Get(TString::Format("H_%i_%i_%i_raw_acceptance_corrected", alpha, L, M)))->Clone(TString::Format("H_%i_%i_%i_raw_acceptance_corrected_statistical", alpha, L, M));
                for (int bin = 1; bin <= analysis_info.n_bins; bin++) {
                    moments_histogram_acceptance_corrected->SetBinError(bin,sqrt(std_devs_acceptance_corrected_statistical[moment_number][bin-1]/num_bootstraps));
                    moments_histogram_raw_acceptance_corrected->SetBinError(bin,sqrt(std_devs_raw_acceptance_corrected_statistical[moment_number][bin-1]/num_bootstraps));
                }
                moments_histogram_acceptance_corrected->Write();
                moments_histogram_raw_acceptance_corrected->Write();
            }
        }
    }
}