#ifndef MEDIA_H_INCLUDED
#define MEDIA_H_INCLUDED
#include "SimParameters.h"
#include "Treatment.h"
#include "Includes/Parameters.h"
#include "RungeKutta4.h"
class APC_cells;
class NK_cells;
class LT_cells;

///Culture Media
class Media
{
public:
    ~Media(){}
    /// 1) Variables (5)
    /// Interpheron gamma concentration in the media
    /*1*/    double& IFNgamma();
    const double& IFNgamma()const;
    /// Tumor Necrosis Factor alpha concentration in the media
    /*2*/    double& TNF();
    const double& TNF()const;
    /// 2) Number of cells contained by the media
    /*3*/    double& num_cells();
    const double& num_cells() const;
    /// Ag concentration
    /*4*/    double& Ag();
    const double& Ag()const;
    /// Tymidine incorporated by cells in media
    /*5*/   double& Tymidine_incorporated();
    const double& Tymidine_incorporated()const;
    /// Tymidine Triate presence in media (0 or 1)
    /**/    double& TymidineTriteate();
    const double& TymidineTriteate() const;
    /// blocking Ab presence in media (0 or 10)
    /**/    double& Ab();
    const double& Ab()const;

    /// 2) Parameters 4
    /// TNF degradation
    /*1*/    double& TNF_degradation ();
    const double& TNF_degradation () const;
    /// IFN degradation
    /*2*/    double& IFN_degradation ();
    const double& IFN_degradation () const;
    /// 3) Ag Parameters
    /*3*/    double& Ag_degradation ();
    const double& Ag_degradation () const;
    /// Timidine Triate/Proliferation conversion factor
    /*4*/    double& Prol_TymTr();
    const double& Prol_TymTr() const;

    void update(double& time_step,double t_run,const APC_cells& APC ,const NK_cells& NK,const LT_cells& LT);

    Media(double IFNgamma_init,
          double TNF_init,
          double Tymidine_incorprated_init,
          double TNF_deg_init,
          double IFN_deg_init,
          double Ag_deg_init,
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


    std::vector<double> Derivative(const APC_cells& APC ,const NK_cells& NK,const LT_cells& LT);

    std::vector<double> getState() const;

    void setState(std::vector<double> y,double t_run,const APC_cells& APC ,const NK_cells& NK,const LT_cells& LT);


private:
    /// variables (5)
    /*1*/ double IFNgamma_d;
    /*2*/ double TNF_d;
    /*3*/ double num_cells_d;
    /*4*/ double Ag_d;
    /*5*/ double Tymidine_incorporated_d;

    /**/ double TymidineTriteate_d; /// presence in media
    /**/ double Ab_d;

    /// Parameters (4)
    /*1*/ double TNF_deg_d;
    /*2*/ double IFN_deg_d;
    /*3*/ double Ag_deg_d;
    /*4*/ double Prol_TymTr_d;/// conversion rate Tym/prol
};


#endif // MEDIA_H_INCLUDED
