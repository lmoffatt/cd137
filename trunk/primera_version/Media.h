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

        void update(double time_step,const APC_cells& APC_ ,const NK_cells& NK,const LT_cells& LT_);

        Media(double max_num_cells_,
              double init_num_cells,
              double Ag_=0,
              double IFNgamma_init=0,
              double TNF_init=0,
              double TNF_deg_init=0,
              double internalization_init=0):

              IFNgamma_d(IFNgamma_init),
              TNF_d(TNF_init),
              max_num_cells_d(max_num_cells_),
              num_cells_d(init_num_cells),
              Ag_d(Ag_),
              TNF_deg(TNF_deg_init),
              Ag_internalization_rate (internalization_init) {}

        Media(){}




    private:
        double IFNgamma_d;
        double TNF_d;
        double max_num_cells_d;
        double num_cells_d;
        double Ag_d;
        double Ag_internalization_rate;
        double TNF_deg;
};



#endif // MEDIA_H_INCLUDED
