#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"

int main(const char* dataset_name) {
    AnalysisInfo analysis_info = getAnalysisInfo(dataset_name);
    std::cout << analysis_info.output_path.Data() << std::endl;
    return 0;
}
