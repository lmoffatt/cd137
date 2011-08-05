#ifndef MEDIA_H_INCLUDED
#define MEDIA_H_INCLUDED

class APC_cells;
class NK_cells;
class LT_cells;

///Culture Media
class Media
{
    public:
        /// Interpheron gamma concentration in the media
        double& IFNgamma();
        const double& IFNgamma()const;

        /// Tumor Necrosis Factor alpha concentration in the media
        double& TNF();
        const double& TNF()const;

        /// TNF degradation
        double& TNF_degradation ();
        const double& TNF_degradation () const;

        /// IFN degradation
        double& IFN_degradation ();
        const double& IFN_degradation () const;

        /// Maximum number of cells tolerated by the media
        /// Proliferation occurs when the number of cells is lower than this number
        double& Max_num_cells();
        const double& Max_num_cells() const;

        /// Number of cells contained by the media
        double& num_cells();
        const double& num_cells() const;

        /// Antigen concentration
        double& Ag();
        const double& Ag()const;

        /// blocking Ab concentration
        double& Ab();
        const double& Ab()const;

        void update(double time_step,
                    const APC_cells& APC_ ,
                    const NK_cells& NK,
                    const LT_cells& LT_);

        Media(double max_num_cells_,
              double init_num_cells,
              double Ag_,
              double Ab_,
              double IFNgamma_init,
              double TNF_init,
              double TNF_deg_init,
              double IFN_deg_init);
     //         double internalization_init);

        Media(){}




    private:
        double IFNgamma_d;
        double TNF_d;
        double max_num_cells_d;
        double num_cells_d;
        double Ag_d;
        double Ab_d;
        double TNF_deg;
        double IFN_deg;
   //     double Ag_internalization_rate;

};



#endif // MEDIA_H_INCLUDED
