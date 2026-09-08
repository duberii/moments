#include <stdexcept>
#include <TString.h>
#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

void acceptance_matrix_parasite(const char* dataset_name, int maxL, int bootstrap_number, int largerMaxL) {
    bool is_bootstrap = (bootstrap_number>-1);

    //Get analysis-specific information

    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    std::unique_ptr<TFile> input_acceptance_matrix_file;
    std::unique_ptr<TFile> output_acceptance_matrix_file;
    if (is_bootstrap) {
        input_acceptance_matrix_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/bootstrap_%i/accMC/acceptance_matrix.root",analysis_info.output_path.Data(), dataset_name,largerMaxL, bootstrap_number), "READ"));
        output_acceptance_matrix_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/bootstrap_%i/accMC/acceptance_matrix.root",analysis_info.output_path.Data(), dataset_name,maxL, bootstrap_number), "UPDATE"));
    } else {
        input_acceptance_matrix_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/accMC/acceptance_matrix.root",analysis_info.output_path.Data(), dataset_name,largerMaxL), "READ"));
        output_acceptance_matrix_file = std::unique_ptr<TFile>(TFile::Open(TString::Format("%s%s/%i/full/accMC/acceptance_matrix.root",analysis_info.output_path.Data(), dataset_name,maxL), "UPDATE"));
    }

    MomentsArray larger_moment_array(largerMaxL);
    MomentsArray smaller_moment_array(maxL);
    output_acceptance_matrix_file->cd();
    for (int bin =0; bin<analysis_info.n_bins; bin++) {
        TMatrixD output_matrix(smaller_moment_array.dimension,smaller_moment_array.dimension);
        TMatrixD* input_matrix = (TMatrixD*)input_acceptance_matrix_file->Get(TString::Format("acceptanceMatrix_%i",bin));
        for (int alpha = 0; alpha < 3; alpha++) {
            int smaller_lower_bound_x;
            int larger_lower_bound_x;
            if (alpha == 0) {
                smaller_lower_bound_x=0;
                larger_lower_bound_x = 0;
            } else {
                smaller_lower_bound_x= smaller_moment_array.getMomentIndex(alpha-1,maxL,maxL)+1;
                larger_lower_bound_x = larger_moment_array.getMomentIndex(alpha-1,largerMaxL,largerMaxL)+1;
            } 
            int larger_upper_bound_x = larger_moment_array.getMomentIndex(alpha, maxL, maxL);
            for (int alphaPrime = 0; alphaPrime < 3; alphaPrime++) {
                int smaller_lower_bound_y;
                int larger_lower_bound_y;
                if (alphaPrime == 0) {
                    smaller_lower_bound_y=0;
                    larger_lower_bound_y = 0;
                } else {
                    smaller_lower_bound_y= smaller_moment_array.getMomentIndex(alphaPrime-1,maxL,maxL)+1;
                    larger_lower_bound_y = larger_moment_array.getMomentIndex(alphaPrime-1,largerMaxL,largerMaxL)+1;
                } 
                int larger_upper_bound_y = larger_moment_array.getMomentIndex(alphaPrime, maxL, maxL);
                output_matrix.SetSub(smaller_lower_bound_x, smaller_lower_bound_y, input_matrix->GetSub(larger_lower_bound_x, larger_upper_bound_x, larger_lower_bound_y, larger_upper_bound_y));
            }
        }
        output_matrix.Write(TString::Format("acceptanceMatrix_%i",bin));
        TH2D matrixHist(output_matrix);
        matrixHist.SetDirectory(nullptr);
        matrixHist.Write(TString::Format("acceptanceMatrixHistogram_%i",bin));
        TMatrixD inverse_matrix = output_matrix;
        inverse_matrix.Invert();
        inverse_matrix.Write(TString::Format("acceptanceMatrixInverse_%i",bin));
        TH2D matrixHistInverse(inverse_matrix);
        matrixHistInverse.SetDirectory(nullptr);
        matrixHistInverse.Write(TString::Format("acceptanceMatrixHistogramInverse_%i",bin));
    }
    input_acceptance_matrix_file->Close();
    output_acceptance_matrix_file->Close();
}