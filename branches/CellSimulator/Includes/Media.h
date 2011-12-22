#ifndef MEDIA_H_INCLUDED
#define MEDIA_H_INCLUDED
#include "SimParameters.h"
#include "Treatment.h"
#include "Includes/Parameters.h"

class APC_cells;
class NK_cells;
class LT_cells;

///Culture Media
class Media
{
    public:
    ~Media(){}
        /// Interpheron gamma concentration in the media
    /*1*/    double& IFNgamma();
             const double& IFNgamma()const;

        /// Tumor Necrosis Factor alpha concentration in the media
    /*2*/    double& TNF();
             const double& TNF()const;

        /// TNF degradation
    /*3*/    double& TNF_degradation ();
             const double& TNF_degradation () const;

        /// IFN degradation
    /*4*/    double& IFN_degradation ();
             const double& IFN_degradation () const;

        /// Number of cells contained by the media
    /*5*/    double& num_cells();
             const double& num_cells() const;

        /// Antigen concentration
    /*6*/    double& Ag();
             const double& Ag()const;

        /// blocking Ab concentration
    /*7*/    double& Ab();
             const double& Ab()const;

    /*8*/    double& TymidineTriteate();
             const double& TymidineTriteate() const;

    /*9*/    double& Prol_TymTr();
             const double& Prol_TymTr() const;

    /*10*/   double& Tymidine_incorporated();
             const double& Tymidine_incorporated()const;


        void update(double& time_step,double t_run,const APC_cells& APC ,const NK_cells& NK,const LT_cells& LT);

        Media(double IFNgamma_init,
              double TNF_init,
              double Tymidine_incorprated_init,
              double TNF_deg_init,
              double IFN_deg_init,
              double init_num_cells,
              double init_Ag,
              double init_Ab,
             double Prol_TymTr_init
              );

        Media(const Parameters& sp,
              const Treatment& tr);

        Media(const Media& other);


        Media(){}

        Media&
        operator=(const Media& other);
        friend void swap(Media& one, Media& other);

	friend std::ostream& operator<<(std::ostream& s, const Media& c);

    private:
        double IFNgamma_d;
        double TNF_d;
        double TNF_deg_d;
        double IFN_deg_d;
        double num_cells_d;
        double Ag_d;
        double Ab_d;
        double TymidineTriteate_d;
        double Prol_TymTr_d;/// conversion rate Tym/prol
        double Tymidine_incorporated_d;



};



#endif // MEDIA_H_INCLUDED
