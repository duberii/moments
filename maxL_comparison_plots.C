#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

void maxL_comparison_plots(const char* dataset_name, std::vector<int> maxLs_to_compare) {
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    std::vector<TString> file_paths;
    std::vector<TString> names;
    for (int maxL : maxLs_to_compare) {
        file_paths.push_back(TString::Format("%s%s/%i/full/data/moments.root", analysis_info.output_path.Data(), dataset_name, maxL));
        names.push_back(TString::Format("L_{max}=%i",maxL));
    }
    compare_results_files(file_paths, names, TString::Format("%s%s/maxL_comparison_plots",analysis_info.output_path.Data(), dataset_name), *std::max_element(maxLs_to_compare.begin(), maxLs_to_compare.end()));

}
