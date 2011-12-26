#ifndef MATERIAL_MODEL_H
#define MATERIAL_MODEL_H

void elasticity_3D(double *D, double lambda, double mu);
void thermoelasticity_3D(double *Db, double lambda, double mu, double alphat);

void cubic_elasticity_3D(double *D, double lambda, double mu,
                         double alpha, double *a, double *b);
void cubic_thermoelasticity_3D(double *Db, double lambda, double mu, 
                               double alpha, double *a, double *b, 
                               double alphat);

void elasticity_2D_strain(double *D, double lambda, double mu);
void thermoelasticity_2D_strain(double *Db, double lambda, double mu, 
                                double alphat);
void elasticity_2D_stress(double *D, double lambda, double mu);
void thermoelasticity_2D_stress(double *Db, double lambda,
                                double mu, double alphat);
void elasticity_axis(double *D, double lambda, double mu);
void thermoelasticity_axis(double *Db, double lambda, double mu, double alphat);


void cubic_elasticity_2D_strain(double *D, double lambda, double mu,
                                double alpha, double *a, double *b);
void cubic_thermoelasticity_2D_strain(double *Db, double lambda, double mu,
                                      double alpha, double *a,
                                      double *b, double alphat);
void cubic_elasticity_2D_stress(double *D, double lambda, double mu,
                                double alpha, double *a, double *b);
void cubic_thermoelasticity_2D_stress(double *Db, double lambda,
                                      double mu, double alpha,
                                      double *a, double *b,
                                      double alphat);

void hex_elasticity_3D(double *D, double *coeff, double *a, double *b);
void hex_elasticity_2D_strain(double *Cs, double *coeff, double *a, double *b);
void hex_elasticity_2D_stress(double *Cs, double *coeff, double *a, double *b);

void piezo_elasticity_3D(double *Cb, double lambda, double mu,
                         double *pz, double kds);
void piezo_elasticity_3D(double *Cb, double *Cmech, 
                         double *pz, double *kds);

void piezo_elasticity_2D_strain(double *Cb, double lambda, double mu,
                                double *pz, double kds);
void piezo_elasticity_2D_strain(double *Cb, double *Cmech, double *pz, 
                                double *kds);
void piezo_elasticity_2D_stress(double *Cb, double lambda, 
                                double mu, double *pz, double kds);

void piezo_elasticity_2D_stress(double *Cb, double *Cmech, double *pz, 
                                double *kds);
void piezo_elasticity_2HD_stress(double *Cb, double lambda, double mu,
                                 double *pz, double kds);

void piezo_elasticity_2HD_stress(double *Cb, double *Cmech, double *pz, 
                                 double *kds);

void hex_piezo_3D(double *pz, double *coeff, double *a, double *b);
void hex_dielectric_3D(double *kds, double *coeff, double *a, double *b);


#endif                          /* MATERIAL_MODEL_H */
