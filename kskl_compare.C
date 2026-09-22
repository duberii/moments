#include "/N/u/rdube/Quartz/work/moments_workflow/shared_funcs.C"
#include <TMatrixD.h>


void kskl_compare() {
    std::vector<TString> figure_names;
    std::unique_ptr<TCanvas> canvas = std::make_unique<TCanvas>("canvas","canvas",700,500);
    figure_names.push_back("H_0_0_0_raw_acceptance_corrected");
    for (int alpha = 0; alpha <3; alpha++) {
        for (int L =0; L <= 6; L++) {
            for (int M =0; M <= L; M++) { 
                if (alpha==2 && M==0) {
                    continue;
                }
                if (alpha == 0 && L==0 && M==0) {
                    continue;
                }
                figure_names.push_back(TString::Format("H_%i_%i_%i_acceptance_corrected", alpha, L, M));
            }
        }
    }
    std::unique_ptr<TFile> spring_2020_file = std::unique_ptr<TFile>(TFile::Open("/N/u/rdube/Quartz/work/moments_workflow/results/kskl_spring_2020/6/full/data/moments.root","READ"));
    std::unique_ptr<TFile> gluex_I_file = std::unique_ptr<TFile>(TFile::Open("/N/u/rdube/Quartz/work/moments_workflow/results/kskl_GlueX_I/6/full/data/moments.root","READ"));
    for (int figure_number = 0; figure_number < int(figure_names.size()); figure_number++) {
        TLegend legend(0.2,0.15);
        TH1D* spring_2020_figure = dynamic_cast<TH1D*>(spring_2020_file->Get(figure_names[figure_number]));
        TH1D* gluex_I_figure = dynamic_cast<TH1D*>(gluex_I_file->Get(figure_names[figure_number]));
        double y_min = min(spring_2020_figure->GetMinimum(),gluex_I_figure->GetMinimum());
        double y_max = max(spring_2020_figure->GetMaximum(),gluex_I_figure->GetMaximum());
        gluex_I_figure->SetMaximum(1.1*y_max);
        gluex_I_figure->SetMinimum(1.1*y_min);
        gluex_I_figure->SetMarkerColor(kRed);
        gluex_I_figure->SetLineColor(kRed);
        spring_2020_figure->SetMarkerColor(kBlue);
        spring_2020_figure->SetLineColor(kBlue);
        legend.AddEntry(gluex_I_figure,"GlueX-I","lp");
        legend.AddEntry(spring_2020_figure,"Spring 2020","lp");
        gluex_I_figure->Draw();
        spring_2020_figure->Draw("SAME");
        legend.Draw();
        canvas->SaveAs(TString::Format("/N/u/rdube/Quartz/work/moments_workflow/kskl_comparison/%s.png",figure_names[figure_number].Data()));
        if (figure_number == 0) {
            spring_2020_figure->Scale(6.28318530718/132.4);
            gluex_I_figure->Scale(6.28318530718/125.0);
            double y_min = min(0.0,min(spring_2020_figure->GetMinimum(),gluex_I_figure->GetMinimum()));
            double y_max = max(spring_2020_figure->GetMaximum(),gluex_I_figure->GetMaximum());
            gluex_I_figure->SetMaximum(1.1*y_max);
            gluex_I_figure->SetMinimum(1.1*y_min);
            gluex_I_figure->SetMarkerColor(kRed);
            gluex_I_figure->SetLineColor(kRed);
            spring_2020_figure->SetMarkerColor(kBlue);
            spring_2020_figure->SetLineColor(kBlue);
            gluex_I_figure->SetTitle("#gamma p -> p K^{0}_{S} K^{0}_{L} Cross Sections; K^{0}_{S} K^{0}_{L} Mass (GeV/c^{2}); Cross Section (pb)");
            gluex_I_figure->Draw();
            spring_2020_figure->Draw("SAME");
            legend.Draw();
            canvas->SaveAs("/N/u/rdube/Quartz/work/moments_workflow/kskl_comparison/cross_sections.png");
        }
    }
    spring_2020_file->Close();
    gluex_I_file->Close();
}
