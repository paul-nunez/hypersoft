
// Compile as: g++ hypersim.c -o hypersim -lm -lcfitsio -lgsl -lgslcblas -lfftw3
//Modified Nov 9 2013

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

# include <fftw3.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define N 1024
#define STRING_LENGTH 100
#define FLAG_LENGTH 6
#define PI 3.141592
#define TINY 1e-30
#define SMALL 1e-5
#define EPSILON 0.0005
#define RAD_TO_MAS 206264849.15927

#include "../Tools/fits.c" //Includes functions to generate fits files
#include "../Tools/misc.c"
#include "../Tools/spot_find.c"

double **sft(int n, int sign, double **re, double **im, double **image); //not used?
double **truncate_fft(int n1, int n2, double **image1, double **image2);
double **pad(double **image_small, double **image_big, int Nmax, int Nmin);
double find_mtf_grid(int n, struct ppl *ppl);
double find_mtf_spacing(int n, struct mtfunc *mtf);
double find_ppl_spacing(int n, struct ppl *ppl);
void find_difference(int n, double **a, double **b);
double **scale_image(int n, double desired_size, double input_size, double **prist, double **scaled);
double **adjust_image_size(int n, double desired_size, double input_size, double **prist);
double **read_pgm(int n, double fov, char *file_name);
fftw_complex *fftw_format(int n, double **real, double **imaginary, fftw_complex *a);
double **fftw_matrix(int n, fftw_complex *a, double **matrix, int flag);
void print_file(int n, int offset, double step,  double **image, char image_file[]);
double **modulus(int n, double **re, double **im);
//double **change_to_2d(int n, double *image, double **data, int q);
double **shift(int n, double **image);
double **transpose(int n, double **image);
double **multiply(int n, double **image1, double **image2);
fftw_complex *autocorrelate_fftw(int n, double **image, double **phase);
double **random_phase(int n, double sigma, double **phase, struct ppl *ppl);
void display(char *name, int n, int offset, double range, double **image);
double **remove_offset(int n, double **image);
double **initialize_matrix(int n, double **image, double init);
double **normalize(int n, double **image);
double find_average(int n, double **image);
double **anti_symmetrize(int n, double **image); //not used?
double **symmetrize(int n, double **image); //not used?
double **generate_pupil(int n, struct ppl *ppl);
void read_pupil(struct ppl *ppl);
void ppl_uncertainty(double sigma, struct ppl *ppl);
double **densify(int n, double gamma, struct mtfunc *mtfunc, double **image); //TODO: STILL ISSUES HERE?
void redefine_mtf(int n, int flag, double gamma, struct mtfunc *mtfunc);
double **poisson(int n, double **image);
double **include_photons(int n, double avgpp, double **image, double average);
void main_spot_find(int n, double threshold, struct ppl *ppl, double **image);
double find_next_gamma(int n, double gamma, struct mtfunc *mtfunc);
void find_mtf_centers(int n, struct ppl *ppl, struct mtfunc *mtf);
double find_error(int n, double *conv_img, double *conv_guess);
int readpilot(struct plt *plt); //in Tools/pilot_read.c
void initialize_pilot(struct plt *plt);

struct plt{
  int psf_flag;
  int fits_flag;
  int fits_npix;
  int pgm_flag;
  int pgm_npix;
  char *pilotfile;
  double lambda;
  char *pristine_name;
  char *res_pristine_name;  
  char *pupil_pgm_name;
  char *gen_par_name;
  int gen_par_flag;
  double pupil_size;
  int pupil_pgm_flag;
  int pupil_fits_flag;
  char *pupil_par_name;
  char *pupil_fits_name;
  char *pupil_fits_read;
  int ppl_uncertain;
  double ppl_err;
  int ppl_phase;
  double ppl_phase_error;
  char *mod_fits_name;
  char *mtf_fits_name;
  int photon;
  double tphoton;
  double nphoton;
  double gamma;
  int dfcorr;
  char *fft_fizeau_name;
  char *fizeau_fits_name;
  char *fft_densified_name;
  char *densified_fits_name;
  char *rl_fits_name;
  double image_size;
};

struct ppl{
  char *file_name;
  int n_sub;
  double pscale;
  int *pos_x;
  int *pos_y;
  double step;
  double size;
  int *diameter;
  int *diameter_c;
  double **pupil;
  double **pupil_centers;
};

struct mtfunc{
  double *pos_x;
  double *pos_y;
  double **centers;
  double *radius;
  int *dpos_x;
  int *dpos_y;
  double *fdpos_x;
  double *fdpos_y;
  int n_sub;
};


int main(int n, char **arg){

  struct plt plt;
  struct ppl ppl;
  struct mtfunc mtfunc;


  if(arg[1]==NULL || strcmp(arg[1], "-h")==0 || n<1){
        fprintf(stderr, "\nPlease use a pilot file with the -p flag\n\n");    
	fprintf(stderr, "Example:\n ./hypersim -p simulation.pilot\n");
	exit(1);
  }    

  initialize_pilot(&plt); //Initialize pilot variables

  if (strcmp(arg[1], "-p")==0){
    strcpy(plt.pilotfile, arg[2]);
    readpilot(&plt);
  }
  else{
    fprintf(stderr, "Pilot file not given?\n");
    exit(1);
  }

  if(plt.pupil_pgm_flag==1) fprintf(stderr, "Pupil .pgm file provided\n");
  if(plt.pupil_fits_flag==1) fprintf(stderr, "Pupil .fits file provided\n");
  if(plt.pupil_pgm_flag==0 && plt.pupil_fits_flag==0) fprintf(stderr, "Pupil .par file provided\n");

  //Allocate memory for pupil
  ppl.pupil=allocate_matrix(N, N, ppl.pupil);
  ppl.pupil_centers=allocate_matrix(N, N, ppl.pupil_centers);
  ppl.pscale=1; //pupil scaling factor can be generalized in case pupil is scaled

  //Read pupil parameter file unless a pupil pgm/fits file is given
  ppl.file_name=new char[STRING_LENGTH];
  strcpy(ppl.file_name, plt.pupil_par_name); 
  if(plt.pupil_pgm_flag==0 && plt.pupil_fits_flag==0){
    read_pupil(&ppl); //This can also scale pupil. TODO: if scaling pupil, would be nice to center it as well
    ppl.pupil=transpose(N, ppl.pupil);  //Transpose pupil (just to be consistent with David Mary's code)
  }
  if(plt.pupil_size!=-1){
    ppl.size=plt.pupil_size;
    ppl.step=ppl.size/N; //TODO: REMOVE THIRD LINE FROM PAR FILES?
  }

  //Generate pupil from parameter .par file unless a pgm file is given
  if(plt.pupil_pgm_flag==0 && plt.pupil_fits_flag==0) ppl.pupil=generate_pupil(N, &ppl);

  //Read pupil FITS file if supplied in pilot file
  if(plt.pupil_fits_flag==1){
    fprintf(stderr, "\nPupil fits read\n");
    double *ppl1=new double[2*N*N];
    ppl1=readfits(plt.pupil_fits_read, N);
    ppl.pupil=change_to_2d(N, ppl1, ppl.pupil, 0);
    ppl.pupil=normalize(N, ppl.pupil);
  }
  
  //Read pupil PGM file if supplied in pilot file
  if(plt.pupil_pgm_flag==1){
    ppl.pupil=read_pgm(N, ppl.size, plt.pupil_pgm_name);
    ppl.pupil=normalize(N, ppl.pupil);   
    ppl.pupil=transpose(N, ppl.pupil); //Transpose pupil (just to be consistent with David Mary's code)
  }

  //Find pupil centers and diameters if requested by pilot
  if(plt.gen_par_flag==1){
    main_spot_find(N, 0.9, &ppl, ppl.pupil);
    //Create pupil parameter file from pgm or fits information
    FILE *test_par=fopen(plt.gen_par_name, "w");
    fprintf(test_par, "%d\n%lf\n", ppl.n_sub, ppl.size);
    for(int i=0; i<ppl.n_sub; i++) fprintf(test_par, "%d %d %d\n", ppl.pos_x[i], ppl.pos_y[i], ppl.diameter[i] );
    fclose(test_par);    
  }
  
  //Set physical pupil size if it is supplied in the pilot file (otherwise 1st line in .par file)
  if(plt.pupil_size!=-1) ppl.size=plt.pupil_size;
  ppl.step=ppl.size/N;  //Pupil pixel size
  fprintf(stderr, "Image size %lf mas\nPupil size %lf \nPupil pixel size %lf\n", plt.image_size, plt.pupil_size, ppl.step);
  plt.lambda=plt.lambda*1E-9; //nanometers

  //NOTE: Maximal size of pupil plane can be defined by the pristine image pixel size and lambda.
  //The pupil size is defined by the user. Sub-aperture pixel positons must be scaled 
  //by the ratio of the maximal pupil size and the user-specified pupil size.  

  double maximal_pupil_size=(plt.lambda*(N/plt.image_size)*RAD_TO_MAS);
  ppl.pscale=1; 
  fprintf(stderr, "Maximal pupil size is : %lf (as defined by pristine pixel size)\n", maximal_pupil_size);

  //sub-pupil location uncertainty 
  int ppl_loc_error=int(floor(plt.ppl_err*N/ppl.size+0.5)); 
  if(plt.ppl_uncertain==1) ppl_uncertainty(ppl_loc_error, &ppl);

  //Pupil is complex in general. One can add phase errors (see lgs1.c) 
  double **pupil_phase=allocate_matrix(N, N, pupil_phase);
  pupil_phase=initialize_matrix(N, pupil_phase, 0);
  if(plt.ppl_phase==1) pupil_phase=random_phase(N, plt.ppl_phase_error, pupil_phase, &ppl);
  
  //Create pupil fits file
  fitsout(N, ppl.step, ppl.pupil,  plt.pupil_fits_name); 

  //Allocate memory for pristine image, and pristine FFT
  char *prist_file=new char[STRING_LENGTH];
  strcpy(prist_file, plt.pristine_name);
  double image_size=plt.lambda/(ppl.step)*RAD_TO_MAS;
  fprintf(stderr, "Pristine image size should be lambda/dx= %0.3f mas,\nas defined by the sampling in the pupil plane \n\n", image_size);

  //Read pgm image  
  double image_step;
  double **prist;  
  double **prist_big;
  int Nmax=0; //Maximal dimension of padded pristine image
  if(plt.psf_flag==0){
    if(plt.fits_flag==1){//Pristine FITS option
      double *prist1=new double[2*plt.fits_npix*plt.fits_npix];
      prist1=readfits(prist_file, plt.fits_npix);
      prist=allocate_matrix(plt.fits_npix, plt.fits_npix, prist);
      prist=change_to_2d(plt.fits_npix, prist1, prist, 0);
      //The image must be padded with zeroes up to Nmax
      double dtheta=plt.image_size/(plt.fits_npix*RAD_TO_MAS);
      image_step=dtheta*RAD_TO_MAS;
      Nmax=int(floor(plt.lambda/(ppl.step*dtheta)+0.5));     
      prist_big=allocate_matrix(Nmax, Nmax, prist_big);
      prist_big=pad(prist, prist_big, Nmax, plt.fits_npix);
      prist_big=shift(Nmax, prist_big); //To have center in 0,0 (not Nmax/2, Nmax/2) 
    }
    if(plt.pgm_flag==1){//Pristine PGM option
      prist=allocate_matrix(plt.pgm_npix, plt.pgm_npix, prist);
      prist=read_pgm(plt.pgm_npix, plt.image_size, prist_file);  
      //The image must be padded with zeroes up to Nmax
      double dtheta=plt.image_size/(plt.pgm_npix*RAD_TO_MAS);
      image_step=dtheta*RAD_TO_MAS;
      Nmax=int(floor(plt.lambda/(ppl.step*dtheta)+0.5));     
      prist_big=allocate_matrix(Nmax, Nmax, prist_big);
      prist_big=pad(prist, prist_big, Nmax, plt.pgm_npix);
      prist_big=shift(Nmax, prist_big); //To have center in 0,0 (not Nmax/2, Nmax/2) 
    }
  }
  //If psf is requested, the pristine image is a point source
  if(plt.psf_flag==1){
    fprintf(stderr, "PSF requested\n");
    prist=allocate_matrix(N, N, prist);
    prist=initialize_matrix(N, prist, 0);
    prist[0][0]=1.0; //Pristine centered at (0,0). Not (512,512)
  }
  
  //NOTE: To resample pristine image I pad with zeros and then truncate the fft  
  //Take fftw  
  if(plt.psf_flag==1) Nmax=N;
  double **im_prist=allocate_matrix(Nmax, Nmax, im_prist);
  im_prist=initialize_matrix(Nmax, im_prist, 0);
  fftw_complex *prist_re;
  if(plt.psf_flag==0) prist_re=fftw_format(Nmax, prist_big, im_prist, prist_re);//Creates input for fftw
  if(plt.psf_flag==1) prist_re=fftw_format(Nmax, prist, im_prist, prist_re);//Creates input for fftw
  fftw_complex *prist_fftw = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * Nmax * Nmax );
  fftw_plan plan_forward = fftw_plan_dft_2d (Nmax, Nmax, prist_re, prist_fftw,FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute ( plan_forward );  
  double **re_fft_prist_big=allocate_matrix(Nmax, Nmax, re_fft_prist_big);
  re_fft_prist_big=fftw_matrix(Nmax, prist_fftw, re_fft_prist_big, 0);
  double **im_fft_prist_big=allocate_matrix(Nmax, Nmax, im_fft_prist_big);
  im_fft_prist_big=fftw_matrix(Nmax, prist_fftw, im_fft_prist_big, 1);
  //Free some memory
  if(plt.psf_flag==0){
    for(int i=0; i<Nmax; i++) delete [] prist_big[i];
    delete [] prist_big;
  }
  
  //Now truncate the fft
  double **re_fft_prist=allocate_matrix(N, N, re_fft_prist);
  double **im_fft_prist=allocate_matrix(N, N, im_fft_prist);
  re_fft_prist=truncate_fft(Nmax, N, re_fft_prist_big, re_fft_prist);
  im_fft_prist=truncate_fft(Nmax, N, im_fft_prist_big, im_fft_prist);   

  //Free some memory
  for(int i=0; i<Nmax; i++){
    delete [] re_fft_prist_big[i];
    delete [] im_fft_prist_big[i];
  }
  delete [] re_fft_prist_big;
  delete [] im_fft_prist_big;
  fftw_free(prist_re);
  fftw_free(prist_fftw);
  fftw_destroy_plan(plan_forward);
  
  //Now display magnitude of the FFT
  double **mod_fft_prist=modulus(N, re_fft_prist, im_fft_prist);
  mod_fft_prist=shift(N, mod_fft_prist);
  fitsout(N, ppl.step, mod_fft_prist, plt.mod_fits_name); 

  //Take inverse fftw to view modified pristine image
  fftw_complex *scaled_prist=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  scaled_prist=fftw_format(N, re_fft_prist, im_fft_prist, scaled_prist);
  fftw_plan p=fftw_plan_dft_2d ( N, N, scaled_prist, scaled_prist, FFTW_BACKWARD,  FFTW_ESTIMATE );
  fftw_execute(p);
  double **new_prist=allocate_matrix(N, N, new_prist);
  new_prist=fftw_matrix(N, scaled_prist, new_prist, 0);
  new_prist=shift(N, new_prist);
  new_prist=normalize(N, new_prist);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(new_prist[i][j]<0) new_prist[i][j]=0;
    }
  }
  //fitsout(N, image_step, new_prist, "Images/Secondary/new_prist.fits");
  fitsout(N, image_step, new_prist, plt.res_pristine_name);
  fftw_free(scaled_prist);
  fftw_destroy_plan(p);
  for(int i=0; i<N; i++) delete [] new_prist[i];
  delete [] new_prist;

  //find mtf centers with pupil information  
  find_mtf_centers(N, &ppl, &mtfunc);

  /*for(int i=0; i<mtfunc.n_sub; i++){
  fprintf(stdout, "%lf %lf\n", mtfunc.pos_x[i], mtfunc.pos_y[i]);
  }*/
    
  //Create OTF from input pupil
  double *otf1=new double[2*N*N];
  double **re_otf=allocate_matrix(N, N, re_otf);
  double **im_otf=allocate_matrix(N, N, im_otf);  
  printf("Autocorrelating pupil\n");
  fftw_complex *pup_auto=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  pup_auto=autocorrelate_fftw(N, ppl.pupil, pupil_phase);
  re_otf=fftw_matrix(N, pup_auto, re_otf, 0);
  re_otf=shift(N, re_otf);
  //im_otf=initialize_matrix(N, im_otf, 0);
  im_otf=fftw_matrix(N, pup_auto, im_otf, 1);
  im_otf=shift(N, im_otf); //Useful only for complex pupils
  double **mtf_diluted=allocate_matrix(N, N, mtf_diluted);  
  mtf_diluted=modulus(N, re_otf, im_otf);
  mtf_diluted=normalize(N, mtf_diluted);

  //Display mtf
  fitsout(N, 1, mtf_diluted, plt.mtf_fits_name);
  double **re_fft_fizeau=allocate_matrix(N,N, re_fft_fizeau);
  double **im_fft_fizeau=allocate_matrix(N,N, im_fft_fizeau);
  re_fft_prist=shift(N, re_fft_prist);
  im_fft_prist=shift(N, im_fft_prist);

  //Multiply by OTF
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      re_fft_fizeau[i][j]=re_fft_prist[i][j]*re_otf[i][j]-im_fft_prist[i][j]*im_otf[i][j];
      im_fft_fizeau[i][j]=re_fft_prist[i][j]*im_otf[i][j]+im_fft_prist[i][j]*re_otf[i][j];
    }
  }
  re_fft_fizeau=shift(N, re_fft_fizeau);
  im_fft_fizeau=shift(N, im_fft_fizeau);


  double **mod_fft_fizeau=allocate_matrix(N, N, mod_fft_fizeau);
  mod_fft_fizeau=modulus(N, re_fft_fizeau, im_fft_fizeau);
  mod_fft_fizeau=normalize(N, mod_fft_fizeau);
  mod_fft_fizeau=shift(N, mod_fft_fizeau);

  //Display fft of fizeau image
  fitsout(N, 1, mod_fft_fizeau, plt.fft_fizeau_name);

  //Densify: To find corrected gamma, we must find minumum mtf spacing
  double dmin=0;
  double gamma;
  if(plt.dfcorr==1){
    dmin=find_ppl_spacing(N, &ppl);
    //double dmin=find_mtf_spacing(N, &mtfunc);
    //double grid=find_mtf_grid(N, &ppl);
    fprintf(stderr, "Minimum mtf spacing: %0.3lf\n", dmin);
    double dprime=floor(dmin-dmin*plt.gamma);
    fprintf(stderr, "Undersampling by %0.2lf\n", dprime);
    gamma=(dmin-dprime)/dmin;
    fprintf(stderr, "Acceptable gamma: %0.8lf\n\n", gamma);
  }
  else gamma=plt.gamma;    

  re_fft_fizeau=shift(N, re_fft_fizeau);//Puts zero freq in N/2,N/2 (easier for densifying)
  im_fft_fizeau=shift(N, im_fft_fizeau);
  double **re_fft_den=allocate_matrix(N, N, re_fft_den);
  double **im_fft_den=allocate_matrix(N, N, im_fft_den);
  for(int j=0; j<N; j++){
    for(int i=0; i<N; i++){
      re_fft_den[i][j]=re_fft_fizeau[i][j];
      im_fft_den[i][j]=im_fft_fizeau[i][j];
    }
  }
  re_fft_den=densify(N, gamma, &mtfunc, re_fft_fizeau); //TODO: problems with ovelapping frequencies 
  im_fft_den=densify(N, gamma, &mtfunc, im_fft_fizeau);

  //Puts zero freq at 0,0 (for computation of fft)
  re_fft_fizeau=shift(N, re_fft_fizeau);
  im_fft_fizeau=shift(N, im_fft_fizeau);
  re_fft_den=shift(N, re_fft_den);
  im_fft_den=shift(N, im_fft_den);

  double **mod_fft_densified=allocate_matrix(N, N, mod_fft_densified);
  mod_fft_densified=modulus(N, re_fft_den, im_fft_den);
  mod_fft_densified=normalize(N, mod_fft_densified);

  //Create fft densified fits file
  mod_fft_densified=shift(N, mod_fft_densified);
  fitsout(N, 1, mod_fft_densified, plt.fft_densified_name);  

  //Inverse fftw (to obtain fizeau)
  fftw_complex *fftw_fizeau=fftw_format(N, re_fft_fizeau, im_fft_fizeau, fftw_fizeau);//Creates input for fftw (matrix to fftw format)
  fftw_complex *fizeau = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );//Allocate output memory in fftw format
  fftw_plan plan_backward = fftw_plan_dft_2d ( N, N, fftw_fizeau, fizeau, FFTW_BACKWARD,  FFTW_ESTIMATE );//Plan optimal fft algthm
  fftw_execute ( plan_backward );  //take fft
  double **fizeau_image=allocate_matrix(N, N, fizeau_image); //Allocate memory for matrix format
  fizeau_image=fftw_matrix(N, fizeau, fizeau_image, 0); //create image matrix
  fizeau_image=normalize(N, fizeau_image);
  fizeau_image=shift(N, fizeau_image);
  fitsout(N, image_step, fizeau_image, plt.fizeau_fits_name);

  //Inverse fftw (to obtain densified image)
  fftw_complex *fftw_densified=fftw_format(N, re_fft_den, im_fft_den, fftw_densified);
  fftw_complex *densified = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  fftw_plan plan_backward_den = fftw_plan_dft_2d ( N, N, fftw_densified, densified, FFTW_BACKWARD,  FFTW_ESTIMATE );//Plan optimal fft algthm
  fftw_execute ( plan_backward_den );  //take fft
  double **densified_image=allocate_matrix(N, N, densified_image); //Allocate memory for matrix format
  densified_image=fftw_matrix(N, densified, densified_image, 0); //create matrix
  densified_image=shift(N, densified_image);
  densified_image=normalize(N, densified_image);

  //****
  //Rescale densified image to simulated inverted galilean telescopes
  Nmax=N/gamma;
  //Pad with zeroes
  double **densified_big=allocate_matrix(Nmax, Nmax, densified_big);
  densified_big=pad(densified_image, densified_big, Nmax, N);
  densified_big=shift(Nmax, densified_big);
  //Now take fft
  double **im_densified_big=allocate_matrix(Nmax, Nmax, im_densified_big);
  im_densified_big=initialize_matrix(Nmax, im_densified_big, 0);
  fftw_complex *densified_big_re;
  densified_big_re=fftw_format(Nmax, densified_big, im_densified_big, densified_big_re);//Creates input for fftw
  fftw_complex *densified_big_fftw = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * Nmax * Nmax );
  fftw_plan plan_fwd = fftw_plan_dft_2d (Nmax, Nmax, densified_big_re, densified_big_fftw,FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute ( plan_fwd );  
  double **re_fft_densified_big=allocate_matrix(Nmax, Nmax, re_fft_densified_big);
  re_fft_densified_big=fftw_matrix(Nmax, densified_big_fftw, re_fft_densified_big, 0);
  double **im_fft_densified_big=allocate_matrix(Nmax, Nmax, im_fft_densified_big);
  im_fft_densified_big=fftw_matrix(Nmax, densified_big_fftw, im_fft_densified_big, 1);
  //Free some memory
  for(int i=0; i<Nmax; i++) delete [] densified_big[i];
  delete [] densified_big;  

  //Now truncate the fft
  double **re_fft_dens=allocate_matrix(N, N, re_fft_dens);
  double **im_fft_dens=allocate_matrix(N, N, im_fft_dens);
  re_fft_dens=truncate_fft(Nmax, N, re_fft_densified_big, re_fft_dens);
  im_fft_dens=truncate_fft(Nmax, N, im_fft_densified_big, im_fft_dens);  

  //Free some memory
  for(int i=0; i<Nmax; i++){
    delete [] re_fft_densified_big[i];
    delete [] im_fft_densified_big[i];
  }
  delete [] re_fft_densified_big;
  delete [] im_fft_densified_big;
  fftw_free(densified_big_re);
  fftw_free(densified_big_fftw);
  fftw_destroy_plan(plan_fwd);

  //Take inverse fftw to view modified densified image
  fftw_complex *scaled_dense=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  scaled_dense=fftw_format(N, re_fft_dens, im_fft_dens, scaled_dense);
  fftw_plan p1=fftw_plan_dft_2d ( N, N, scaled_dense, scaled_dense, FFTW_BACKWARD,  FFTW_ESTIMATE );//OJO repeticion?
  fftw_execute(p1);
  double **new_dense=allocate_matrix(N, N, new_dense);
  new_dense=fftw_matrix(N, scaled_dense, new_dense, 0);
  new_dense=shift(N, new_dense);
  new_dense=normalize(N, new_dense);
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      if(new_dense[i][j]<0) new_dense[i][j]=0;
    }
  }
  fitsout(N, image_step, new_dense, "./Data/new_dense.fits");
  fftw_free(scaled_dense);
  fftw_destroy_plan(p1);
  //for(int i=0; i<N; i++) delete [] new_dense[i];
  //delete [] new_dense;

  //***


  //Positivity. Remove negative values. These seem to result from improper densification (by +-1 pix) of some mtf spots
  //for(int i=0; i<N; i++) for(int j=0; j<N; j++) if(densified_image[i][j]<0) densified_image[i][j]=0;

  //Add noise
  double avgpp=plt.nphoton; //Average photons per pixel 
  //TODO: include detector sampling
  double **sim_data=allocate_matrix(N, N, sim_data);
  double **thermal_noise=allocate_matrix(N, N, sim_data);
  if(plt.photon==1){ 
    double Tnoise=plt.tphoton/plt.nphoton; //Fraction of thermal photons
    thermal_noise=initialize_matrix(N, thermal_noise, Tnoise);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
	sim_data[i][j]=fizeau_image[i][j]+thermal_noise[i][j];
      }
    }
    sim_data=normalize(N, sim_data);
    double avg=find_average(N, sim_data);
    sim_data=include_photons(N, avgpp, sim_data, avg);//Simulate data with densification
    /*  double Tnoise=plt.tphoton/plt.nphoton; //Percent of average number of photons of source
    thermal_noise=initialize_matrix(N, thermal_noise, Tnoise);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
	sim_data[i][j]=new_dense[i][j]+thermal_noise[i][j];
      }
    }
    double avg=find_average(N, sim_data);
    sim_data=include_photons(N, avgpp, sim_data, avg);//Simulate data with densification
    */
  }
  if(plt.photon!=1) sim_data=normalize(N, densified_image); 

  double step=plt.lambda/(ppl.step*N)*RAD_TO_MAS;//Image step //TODO: erase this?
  fitsout(N, image_step, sim_data, plt.densified_fits_name); //Simulated data
  //print_file(N, 0, N*image_step,  sim_data, "Data/sim_data.txt");


  //For taking the autocorrelation of a speckle pattern
  /*double **auto_corr=allocate_matrix(N, N, auto_corr);
  double **zero=allocate_matrix(N, N, zero);
  zero=initialize_matrix(n, zero, 0);  
  fftw_complex *acorr=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  acorr=autocorrelate_fftw(N, fizeau_image, zero);
  auto_corr=fftw_matrix(N, acorr, auto_corr, 0);
  char *name=new char[50];
  strcpy(name, "./Images/Secondary/corr.fits");
  auto_corr=shift(N, auto_corr);
  fitsout(N, 1, auto_corr, name);*/

  return 0;
    
}

//=====================================================================
//=====================================================================

double find_mtf_grid(int n, struct ppl *ppl){

  int d_min=1;
  int counter=0;
  for(int d=1; d<n; d++){
    for(int k=0; k<ppl->n_sub; k++){
      for(int i=0; i<n; i=i+d){
	for(int j=0; j<n; j=j+d){
	  if(ppl->pupil_centers[i][j]==1) counter++;
	}
      }  
    }
    fprintf(stderr, "%d counter %d %d\n", d, counter, ppl->n_sub);
    counter=0;
  }
  exit(1);
  return d_min;
}

//=====================================================================
//=====================================================================

double find_mtf_spacing(int n, struct mtfunc *mtf){
  double d=0;
  double d_min=1024;
  int k1=0;
  int k2=0;
  for(int k=0; k<mtf->n_sub; k++){
    for(int l=0; l<mtf->n_sub; l++){
      d=sqrt(pow(mtf->pos_x[k]-mtf->pos_x[l],2)+pow(mtf->pos_y[k]-mtf->pos_y[l],2));
      double ratio1=(mtf->pos_x[k]-n/2)/sqrt(pow((mtf->pos_x[k]-n/2),2)+pow((mtf->pos_y[k]-n/2),2));
      double ratio2=(mtf->pos_x[l]-n/2)/sqrt(pow((mtf->pos_x[l]-n/2),2)+pow((mtf->pos_y[l]-n/2),2));
      if(d<d_min && d!=0 && ratio1==ratio2){
	d_min=d;
	k1=k;
	k2=l;
	if(d<mtf->radius[k]){
	  fprintf(stderr, "MTF SEPARATION TOO SMALL. exiting program\n");
	  exit(1);	  
	}
      }
    }    
  }
  fprintf(stderr, "Smallest mtf distance between (%0.2lf, %0.2lf) and (%0.2lf, %0.2lf)\n", mtf->pos_x[k1], mtf->pos_y[k1], mtf->pos_x[k2], mtf->pos_y[k2] );
  return d_min;
}

//=====================================================================
//=====================================================================


double find_ppl_spacing(int n, struct ppl *ppl){
  double d=0;
  double d_min=1024;
  int k1=0;
  int k2=0;
  for(int k=0; k<ppl->n_sub; k++){
    for(int l=0; l<ppl->n_sub; l++){
      d=sqrt(pow(ppl->pos_x[k]-ppl->pos_x[l],2)+pow(ppl->pos_y[k]-ppl->pos_y[l],2));
      double ratio1=(ppl->pos_x[k]-n/2)/sqrt(pow((ppl->pos_x[k]-n/2),2)+pow((ppl->pos_y[k]-n/2),2));
      double ratio2=(ppl->pos_x[l]-n/2)/sqrt(pow((ppl->pos_x[l]-n/2),2)+pow((ppl->pos_y[l]-n/2),2));
      if(d<d_min && d!=0 && ratio1==ratio2){
	//d_min=d;
	if(pow(ppl->pos_x[k]-ppl->pos_x[l],2)<pow(ppl->pos_y[k]-ppl->pos_y[l],2) && pow(ppl->pos_x[k]-ppl->pos_x[l],2)>0) d_min=sqrt(pow(ppl->pos_x[k]-ppl->pos_x[l],2));
	if(pow(ppl->pos_y[k]-ppl->pos_y[l],2)<=pow(ppl->pos_x[k]-ppl->pos_x[l],2) && pow(ppl->pos_y[k]-ppl->pos_y[l],2)>0) d_min=sqrt(pow(ppl->pos_y[k]-ppl->pos_y[l],2));
	//else d_min=sqrt(pow(ppl->pos_y[k]-ppl->pos_y[l],2));
	k1=k;
	k2=l;
	if(d<ppl->diameter[k]){
	  fprintf(stderr, "PPL SEPARATION TOO SMALL. exiting program\n");
	  exit(1);	  
	}
      }
    }    
  }
  fprintf(stderr, "Smallest distance between sub-pupils (%d, %d) and (%d, %d) corresponds to a grid with spacing %d\n", ppl->pos_x[k1], ppl->pos_y[k1], ppl->pos_x[k2], ppl->pos_y[k2], int(floor(d_min+0.5)) );

  return d_min;
}


//=====================================================================
//=====================================================================

double **sft(int n, int sign, double **re, double **im, double **image){
  //Inverse fourier transform
  //VERY SLOW!! useful for debugging
  fprintf(stderr, "Doing slow fourier transform\n");

  // Only do real part (imaginary part should be near zero). for each k and l, a sum over i and j is made 
  float sum, sum2;
  int k, l, i, j;

  k=0; l=0;
  for(k = 0; k < n; ++k) {
    for( i = 0, sum = 0,  sum2=0; i < n; i++) {
      for(j = 0; j < n; j++){
	sum = re[i][j]*cos(sign*2*PI*(k*i+l*j)/n)-im[i][j]*sin(sign*2*PI*(k*i+l*j)/n)+sum;
	//sum2= re[i][j]*sin(sign*2*PI*(k*i+l*j)/n)+im[i][j]*cos(sign*2*PI*(k*i+l*j)/n)+sum2;
      }
    }
    image[k][l] = sum;
    printf("%d  %d  %f\n", k, l , image[k][l]); //This print can be piped out
  }

  //im_image[k][l] = sum2/n;
  //if(im_image[k][l]/image[k][l]>0.0001) printf("%d  %d has imaginary part\n", k, l);
  exit(1);

  for(k = 0; k < n; ++k) {if(k==n/4)fprintf(stderr, "done 25 percent\n");
      for(l = 0; l < n; ++l) {
	for( i = 0, sum = 0,  sum2=0; i < n; i++) {
	  for(j = 0; j < n; j++){
	    sum = re[i][j]*cos(sign*2*PI*(k*i+l*j)/n)-im[i][j]*sin(sign*2*PI*(k*i+l*j)/n)+sum;
	    //sum2= re[i][j]*sin(sign*2*PI*(k*i+l*j)/n)+im[i][j]*cos(sign*2*PI*(k*i+l*j)/n)+sum2;
	  }
	}
	image[k][l] = sum/n;
	printf("%d  %d  %f\n", k, l , image[k][l]); //This print can be piped out
	//im_image[k][l] = sum2/n;
	//if(im_image[k][l]/image[k][l]>0.0001) printf("%d  %d has imaginary part\n", k, l);
      }
      printf("%d ", k);
    }
    
  fprintf(stderr,"done sft\n\n");
  
  return image;
}

//=====================================================================
//=====================================================================


double **anti_symmetrize(int n, double **image){

  for(int j=0; j<n/2; j++){
    for(int i=0; i<n/2; i++){
      //image[n-i-1][n-j-1]=-image[i+1][j+1];
      //image[n-i-1][j+1]=-image[i+1][n-j-1];

      image[n/2-1-i][n/2-1-j]=-image[n/2+i][n/2+j];
      image[n-1-i][j]=-image[i][n-1-j];

    }
  }

  return image;

}

//=====================================================================
//=====================================================================

double **symmetrize(int n, double **image){

  for(int j=0; j<n/2; j++){
    for(int i=0; i<n/2; i++){
      //image[n-i-1][n-j-1]=image[i+1][j+1];
      //image[n-i-1][j+1]=image[i+1][n-j-1];

      image[n/2-1-i][n/2-1-j]=image[n/2+i][n/2+j];
      image[n-1-i][j]=image[i][n-1-j];

    }
  }

  return image;

}


//=====================================================================
//=====================================================================

double find_error(int n, double *ref, double *conv_guess){

  double sum=0;
  double sum1=0;
  for(int i=0; i<n; i++){
    sum=sum+(ref[i]-conv_guess[i])*(ref[i]-conv_guess[i]);
    sum1=sum1+ref[i]*ref[i];
  }
  double error=0;
  if(sum1=!0) error=sum/sum1;
  else{fprintf(stderr, "Division by 0 in find_data_error\n"); exit(1);}

  return error;
}

//=====================================================================
//=====================================================================

void find_difference(int n, double **a, double **b){
  int k=0;
  double max=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      //fprintf(stdout, "%d %d %lf %lf\n", i, j, a[i][j], b[i][j]); 
      if(sqrt((a[i][j]-b[i][j])*(a[i][j]-b[i][j]))>0.00001){ 
	//fprintf(stdout, "%d %d %lf\n", i, j, a[i][j]-b[i][j]); 
	k++;
	//if(sqrt((a[i][j]-b[i][j])*(a[i][j]-b[i][j])>max)) max=sqrt((a[i][j]-b[i][j])*(a[i][j]-b[i][j]))/(a[i][j]+b[i][j])*100;
	if(sqrt((a[i][j]-b[i][j])*(a[i][j]-b[i][j])>max)) max=sqrt((a[i][j]-b[i][j])*(a[i][j]-b[i][j]));
      }
    }
  }
      fprintf(stdout, "Found %d differences as big as %lf\n", k, max);


}


//=====================================================================
//=====================================================================

void find_mtf_centers(int n, struct ppl *ppl, struct mtfunc *mtf){
  fprintf(stderr, "Defining mtf 'spot' positions\n");
  //To find mtf spot positions, we can take the autocorrelation of pupil positions

  //First create an image that is 1 when there is a sub-pupil, and 0 otherwise
  double **pupil_centers=allocate_matrix(n, n, pupil_centers);
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      pupil_centers[i][j]=0; //Initialize
    }
  }
  
  for(int k=0; k<ppl->n_sub; k++){
    for(int j=0; j<n; j++){
      for(int i=0; i<n; i++){
	if(i==ppl->pos_x[k] && j==ppl->pos_y[k]) pupil_centers[i][j]=1;
      }
    }
  }

  //fitsout(n, 1, ppl->pupil, "test.fits");
  //exit(1);
  //Now find the autocorrelation of pupil_centers to find the mtf spot centers
  double **mtf_centers=allocate_matrix(n, n, mtf_centers);
  double **zero=allocate_matrix(n, n, zero);
  zero=initialize_matrix(n, zero, 0);
  fftw_complex *mtfc=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  mtfc=autocorrelate_fftw(n, pupil_centers, zero);
  mtf_centers=fftw_matrix(n, mtfc, mtf_centers, 0);
  mtf_centers=normalize(n, mtf_centers);
  mtf_centers=shift(n, mtf_centers);//This is to find spots with zero at 512

  char *mtf_centers_name=new char[STRING_LENGTH];
  strcpy(mtf_centers_name, "./Images/Secondary/mtf_centers.fits");
  fitsout(N, 1, mtf_centers, mtf_centers_name);

  //Find number of spots
  mtf->n_sub=0;
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      if(mtf_centers[i][j]>SMALL) mtf->n_sub=mtf->n_sub+1;
      //if(mtf_centers[i][j]>0.0121) mtf->n_sub=mtf->n_sub+1;
    }
  }
  fprintf(stderr, "nspot %d\n", mtf->n_sub);
  //Now find mtf spot positions
  mtf->pos_x=new double[mtf->n_sub];
  mtf->pos_y=new double[mtf->n_sub];
  mtf->radius=new double[mtf->n_sub];
  mtf->centers=allocate_matrix(n, n, mtf->centers);
  mtf->centers=initialize_matrix(n, mtf->centers, 0);
  int k=0;
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      if(mtf_centers[i][j]>SMALL){
	mtf->pos_x[k]=i;
	mtf->pos_y[k]=j;
	mtf->centers[i][j]=1;
	k++;
      }
    }
  }
  for(int i=0; i<mtf->n_sub; i++) mtf->radius[i]=ppl->diameter[0]; //TODO: must generalize to mtf radii of variable size!
}

//=====================================================================
//=====================================================================

void redefine_mtf(int n, int flag, double gamma, struct mtfunc *mtfunc){
  for(int k=0; k<mtfunc->n_sub; k++){
    if(flag==1){
      mtfunc->pos_y[k]=mtfunc->pos_y[k]-mtfunc->dpos_y[k];
      mtfunc->pos_x[k]=mtfunc->pos_x[k]-mtfunc->dpos_x[k];
    }
    if(flag==-1){
      mtfunc->pos_y[k]=mtfunc->pos_y[k]+mtfunc->dpos_y[k];
      mtfunc->pos_x[k]=mtfunc->pos_x[k]+mtfunc->dpos_x[k];      
    }
  }
}

//=====================================================================
//=====================================================================

double find_average(int n, double **image){
  double sum=0;
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      sum=sum+image[i][j];
    }
  }
  double average=sum/(n*n);
  return average;
}

//=====================================================================
//=====================================================================

 double **include_photons(int n, double avgpp, double **image, double average){
   fprintf(stderr, "Including %0.2e photons", avgpp);
   //First find sum of image
   double sum=0;
   for(int j=0; j<n; j++) for(int i=0; i<n; i++) sum=sum+image[i][j];     

   double max=0;     
   //Now include number of photons per pixel without noise
   for(int j=0; j<n; j++){
     for(int i=0; i<n; i++){
       //image[i][j]=image[i][j]*avgpp/average;
       image[i][j]=image[i][j]/sum*avgpp;
       if(image[i][j]>max) max=image[i][j];
     }
   }
   fprintf(stderr, "\n\nMAX=%lf \n\n ", max);

   //Now add Poisson noise
   image=poisson(N, image); 

   return image;

 }
/*
double **include_photons(int n, double avgpp, double **image, double average){
  fprintf(stderr, "Including %0.2e photons per pixel", avgpp);
  //First find average of image
  double sum=0;
  image=normalize(n, image);
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      sum=sum+image[i][j];
    }
  }
  fprintf(stderr, "sum %lf\n", sum);
  //double average=sum/(n*n);
  //sum=1;

  //Normalize image so that sum of photons =1 
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      image[i][j]=image[i][j]/sum;
    }
  }

  //Now include number of photons per pixel without noise
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      //image[i][j]=image[i][j]*avgpp/average;
      image[i][j]=image[i][j]*avgpp;
    }
  }

  //Now add Poisson noise
  image=poisson(N, image); 

  return image;

}
*/
//=====================================================================
//=====================================================================


double **poisson(int n, double **image){
  //Taken from http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Distribution-Examples.html
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  long seed;
  seed = time (NULL);    
  gsl_rng_set (r, seed);       

  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      image[i][j]=gsl_ran_poisson (r, image[i][j]);
    }
  }

  return image;

}


//=====================================================================
//=====================================================================

double **scale_image(int n, double desired_size, double input_size, double **prist, double **scaled){
  //TODO: when the desired image is very big, it appears truncated, because it samples up to n/s
  scaled=allocate_matrix(n, n, scaled);
  double **scaled2=allocate_matrix(n,n,scaled2);
  double s=desired_size/input_size;
  if(s<1){
    fprintf(stderr, "\nPlease supply an image with (total) angular extent greater or equal to %0.3f mas\n\n", desired_size);
    //exit(1);
  };
  int d=int(floor((n-n/s)/2+0.5));
  //int d=0;
  fprintf(stderr, "Scaling image by %0.5f\n", s);
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      int is=int(i*s);
      int js=int(j*s);
      if(is<n && js<n && is>=0 && js>=0) scaled[i][j]=prist[is][js];      
      else scaled[i][j]=0;
    }
  }
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      //I don't understand why I can't do the translation in the previous loops
      if(i-d>=0 && i-d<n && j-d>=0 && j-d<n) scaled2[i][j]=scaled[i-d][j-d]; 
      //if(i-d>=0 && j-d>=0) scaled2[i][j]=scaled[i-d][j-d];
      else scaled2[i][j]=0;
    }
  } 

  return scaled2;

}


//=====================================================================
//=====================================================================

double find_next_gamma(int n, double gamma, struct mtfunc *mtfunc){
  fprintf(stderr, "Input gamma= %lf\n", gamma);
  mtfunc->fdpos_x=new double[mtfunc->n_sub];
  mtfunc->fdpos_y=new double[mtfunc->n_sub];
  int flag=0;
  int count=0;
  while(count<mtfunc->n_sub/2-1){
    count=0;
    for(int k=0; k<mtfunc->n_sub; k++){
      double d=sqrt((mtfunc->pos_x[k]-n/2)*(mtfunc->pos_x[k]-n/2)+(mtfunc->pos_y[k]-n/2)*(mtfunc->pos_y[k]-n/2));
      int dx=int(floor(mtfunc->pos_x[k]-n/2+0.5));
      int dy=int(floor(mtfunc->pos_y[k]-n/2+0.5));
      double alpha=(mtfunc->pos_y[k]-n/2)/(mtfunc->pos_x[k]-n/2);//Conserved quantity
      if(abs(dx)>0){ 
	mtfunc->fdpos_x[k]=((((mtfunc->pos_x[k]-n/2)/abs(mtfunc->pos_x[k]-n/2))*(gamma*d/sqrt(1+alpha*alpha))));
	mtfunc->fdpos_y[k]=((((mtfunc->pos_y[k]-n/2)/abs(mtfunc->pos_y[k]-n/2))*(gamma*sqrt(alpha*alpha)*d/sqrt(1+alpha*alpha))));
      }
      //printf("h\n");
      if(dy==0){
	mtfunc->fdpos_x[k]=((((mtfunc->pos_x[k]-n/2)/abs(mtfunc->pos_x[k]-n/2))*(gamma*d/sqrt(1+alpha*alpha))));
	mtfunc->fdpos_y[k]=((gamma*sqrt(alpha*alpha)*d/sqrt(1+alpha*alpha)+0.5));          
      }
      if(dx==0){ 
	mtfunc->fdpos_x[k]=0;
	mtfunc->fdpos_y[k]=((gamma*d));
      }
      if(dx==0 && dy<0){ 
	mtfunc->fdpos_x[k]=0;
	mtfunc->fdpos_y[k]=-((gamma*d));
      }
      if(sqrt(pow((mtfunc->fdpos_x[k]-floor(mtfunc->fdpos_x[k]+0.5)),2))<EPSILON && sqrt(pow((mtfunc->fdpos_y[k]-floor(mtfunc->fdpos_y[k]+0.5)),2))<EPSILON) count++;
      if(sqrt(pow((mtfunc->fdpos_x[k]-floor(mtfunc->fdpos_x[k]+0.5)),2))>EPSILON || sqrt(pow((mtfunc->fdpos_y[k]-floor(mtfunc->fdpos_y[k]+0.5)),2))>EPSILON) gamma=gamma+0.0000001;
    }
    flag ++;
    //fprintf(stderr, "%d count %d, %d\n", flag, count, mtfunc->n_sub-1);
  }
  fprintf(stderr, "%d Output gamma= %lf %d %d\n", flag, gamma, count, mtfunc->n_sub-1);
  return gamma;
}


//=====================================================================
//=====================================================================


double **densify(int n, double gamma, struct mtfunc *mtfunc, double **image){
  //TODO:  must deal with overlapping frequencies
  fprintf(stderr, "Densifying\n");

  double **image_d=allocate_matrix(n, n, image_d);
  image_d=initialize_matrix(n, image_d, 0);
  
  mtfunc->dpos_x=new int[mtfunc->n_sub];
  mtfunc->dpos_y=new int[mtfunc->n_sub];

  //gamma=find_next_gamma(n, gamma, mtfunc);

  for(int k=0; k<mtfunc->n_sub; k++){
    //double a=1.1;
    double a=1.25; //This has been optimized so that gamma of 1 gives exactly the same psf
    double d=sqrt((mtfunc->pos_x[k]-n/2)*(mtfunc->pos_x[k]-n/2)+(mtfunc->pos_y[k]-n/2)*(mtfunc->pos_y[k]-n/2));
    int dx=int(floor(mtfunc->pos_x[k]-n/2+0.5));
    int dy=int(floor(mtfunc->pos_y[k]-n/2+0.5));

    mtfunc->dpos_x[k]=int(floor((mtfunc->pos_x[k]-n/2)-(gamma*(mtfunc->pos_x[k]-n/2))+0.5));
    mtfunc->dpos_y[k]=int(floor((mtfunc->pos_y[k]-n/2)-(gamma*(mtfunc->pos_y[k]-n/2))+0.5));

    int start_y=int(floor(mtfunc->pos_y[k]-a*mtfunc->radius[k]+0.5));
    int start_x=int(floor(mtfunc->pos_x[k]-a*mtfunc->radius[k]+0.5));
    int end_y=int(floor(mtfunc->pos_y[k]+a*mtfunc->radius[k]+0.5));
    int end_x=int(floor(mtfunc->pos_x[k]+a*mtfunc->radius[k]+0.5));
    for(int j=start_y; j<=end_y; j++){
      for(int i=start_x; i<=end_x; i++){
	if(i>0 && i<n && j>0 &&j<n){
	  if(i-mtfunc->dpos_x[k]>=0 && i-mtfunc->dpos_x[k]<n && j-mtfunc->dpos_y[k]>=0 && j-mtfunc->dpos_y[k]<n){
	    image_d[i-mtfunc->dpos_x[k]][j-mtfunc->dpos_y[k]]=image[i][j]; 
	  }
	}
      }
    }

  }
  return image_d;
  
}

//=====================================================================
//=====================================================================

double **generate_pupil(int n, struct ppl *ppl){
  fprintf(stderr, "Generating pupil from .par file\n");
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      ppl->pupil[i][j]=0;
      ppl->pupil_centers[i][j]=0;
    }
  }
  
  for(int i=0; i<ppl->n_sub; i++){
    for(int l=ppl->pos_y[i]-int(ppl->diameter[i]/2); l<=ppl->pos_y[i]+int(ppl->diameter[i]/2); l++){
      for(int k=ppl->pos_x[i]-int(ppl->diameter[i]/2); k<=ppl->pos_x[i]+int(ppl->diameter[i]/2); k++){
	if((k-ppl->pos_x[i])*(k-ppl->pos_x[i])+(l-ppl->pos_y[i])*(l-ppl->pos_y[i])<=int(ppl->diameter[i]/2)*int(ppl->diameter[i]/2)){
	  ppl->pupil[k][l]=1;
	}
      }
    }
  }

  return ppl->pupil;

}

//=====================================================================
//=====================================================================

//=====================================================================
//=====================================================================


double **transpose(int n, double **image){
  double **timage=allocate_matrix(N, N, timage);
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
	timage[i][j]=image[j][i];
      }
    }
    for(int i=0; i<N; i++){
      for(int j=0; j<N; j++){
	image[i][j]=timage[i][j];
      }
    }
    for(int i=0; i<N; i++) delete [] timage[i];
    delete [] timage;    
    
    return image;

}

//=====================================================================
//=====================================================================

double **normalize(int n, double **image){
  double max=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(image[i][j]>max) max=image[i][j];
    }
  }
  
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      image[i][j]=image[i][j]/max;
    }
  }

  return image;

}

//=====================================================================
//=====================================================================

double **remove_offset(int n, double **image){
  double min=1e30;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(image[i][j]<min) min=image[i][j];
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      image[i][j]=image[i][j]-min;
    }
  }

  return image;
}


//=====================================================================
//=====================================================================

void display(char *name, int n, int offset, double range, double **image){
  fprintf(stderr, "\nCalling display function\n");
  FILE *settings;

  char *txt_name=new char[STRING_LENGTH];
  strcpy(txt_name, name);
  strcat(txt_name,".txt");
  print_file(n, offset, range, image, txt_name);

  char *settings_name=new char[STRING_LENGTH];
  strcpy(settings_name, name);
  strcat(settings_name, ".plt");
  
  char *ps_name=new char[STRING_LENGTH];
  strcpy(ps_name, name);
  strcat(ps_name, ".ps");

  char *pdf_name=new char[STRING_LENGTH];
  strcpy(pdf_name, name);
  strcat(pdf_name, ".pdf");

  double step=double(range/N);
  double plot_range=n*step;
  settings=fopen(settings_name, "w");
  //fprintf(settings, "set terminal postscript\nset output '%s'\nset xrange [-%lf:%lf]\nset yrange [-%lf:%lf]\nset pm3d map\nset size ratio 1\nsplot '%s' u 1:2:3\n", ps_name, plot_range, plot_range, plot_range, plot_range, txt_name);
  fprintf(settings, "set terminal postscript\nset output '%s'\nset pm3d map\nset size ratio 1\nsplot '%s' u 1:2:3\n", ps_name, txt_name);
  fclose(settings);
  
  char *plot_command=new char[100];
  sprintf(plot_command, "gnuplot %s", settings_name);
  system(plot_command);

  char *convert_command=new char[100];
  sprintf(convert_command, "ps2pdf %s %s", ps_name, pdf_name);
  system(convert_command);

  char *command=new char[100];
  sprintf(command,"xpdf %s &", pdf_name);
  //system(command);

}

//=====================================================================
//=====================================================================

//double **autocorrelate_fftw(int n, double **image, double **phase){
fftw_complex *autocorrelate_fftw(int n, double **image, double **phase){
  fprintf(stderr, "Calculating autocorrelation\n");
  //Take fftw of image
  double **im=allocate_matrix(n, n, im);  

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      image[i][j]=image[i][j]*cos(phase[i][j]);
      im[i][j]=image[i][j]*sin(phase[i][j]);
    }
  }
  //im=initialize_matrix(n, im, 0); //TODO: GENERALIZE TO COMPLEX PUPIL

  fftw_complex *input=fftw_format(n, image, im, input);
  fftw_complex *output = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * n * n );
  fftw_plan plan_forward = fftw_plan_dft_2d ( n, n, input, output, FFTW_FORWARD,  FFTW_ESTIMATE );
  fftw_execute ( plan_forward );  
  double **re=allocate_matrix(n, n, re);
  re=fftw_matrix(N, output, re, 0);
  im=fftw_matrix(N, output, im, 1);
  double **mod=allocate_matrix(n, n, mod);
  mod=modulus(n, re, im);  

  double **conv=allocate_matrix(n, n, conv);
  conv=multiply(n, mod, mod);
  conv=normalize(n, conv);

  im=initialize_matrix(n, im, 0);
  input=fftw_format(N, conv, im, input);
  fftw_plan plan_backward = fftw_plan_dft_2d ( n, n, input, output, FFTW_BACKWARD,  FFTW_ESTIMATE );
  fftw_execute ( plan_backward);  

  //conv=fftw_matrix(n, output, conv, 0);

  /*double **imaginary=allocate_matrix(n, n, imaginary);
  imaginary=fftw_matrix(n, output, imaginary, 1);
  char *name=new char[50];
  strcpy(name, "imaginary.fits"); 
  fitsout(N, 1, imaginary,  name); */ //To check that imaginary part is zero

  //return conv;

  return output;

}

//=====================================================================
//=====================================================================

double **multiply(int n, double **image1, double **image2){
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      image1[i][j]=image1[i][j]*image2[i][j];
    }
  }
  return image1;
}


//=====================================================================
//=====================================================================
/*
double **change_to_2d(int n, double *image, double **data, int q){
  //q sets the starting point for either real (q=0) or imaginary image (q=1)
  int k=q;
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      data[i][j]=image[k];
      if(isnan(data[i][j])){
	fprintf(stderr, "NAN IMG\n");
	exit(1);
      }
      k=k+2;
    }
  }
  return data;

}
*/

//=====================================================================
//=====================================================================

double **pad(double **image_small, double **image_big, int Nmax, int Nmin){

  printf("Padding up to %d (Original size %d)\n", Nmax, Nmin);
  image_big=initialize_matrix(Nmax, image_big, 0);
  for(int i=0; i<Nmax; i++){
    for(int j=0; j<Nmax; j++){
      if(i<Nmin && j<Nmin && i+(Nmax/2-Nmin/2)<Nmax && j+(Nmax/2-Nmin/2)<Nmax){
	image_big[i+(Nmax/2-Nmin/2)][j+(Nmax/2-Nmin/2)]=image_small[i][j];
      }
    }
  }

  return image_big;

}


//====================================================================
//====================================================================

double **shift(int n, double **image){
  double **centered_image=allocate_matrix(n,n, centered_image);
  for(int j=0; j<n; j++){
      for(int i=0; i<n; i++){	
      if(i<n/2 && j<n/2) centered_image[i][j]=image[i+n/2][j+n/2];  //First quadrant      
      if(i>=n/2 && j<n/2) centered_image[i][j]=image[i-n/2][j+n/2]; //second quadrant
      if(i<n/2 && j>=n/2) centered_image[i][j]=image[i+n/2][j-n/2]; //third
      if(i>=n/2 && j>=n/2) centered_image[i][j]=image[i-n/2][j-n/2];//fourth
    }
  }

  for(int i=0; i<n; i++) for(int j=0; j<n; j++) image[i][j]=centered_image[i][j];
  for(int i=0; i<n; i++)  delete  centered_image[i];
  delete  centered_image;

  return image;


}


//=====================================================================
//=====================================================================


double **modulus(int n, double **re, double **im){

  double **mod=allocate_matrix(n, n, mod);
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      mod[i][j]=re[i][j]*re[i][j]+im[i][j]*im[i][j];
      mod[i][j]=sqrt(mod[i][j]);
    }
  }
  return mod;
  for(int i=0; i<n; i++){
    delete mod[i];
  }
  delete mod;

}

//=====================================================================
//=====================================================================


void print_file(int n, int offset, double range,  double **image, char image_file[]){
  //This function recieves the length of one side of the image matrix 'n' (the total size is n*n),
  //it also recieves  the image matrix itself 'image', and the file name.
  //Then it simply prints to a file whose name is contained in 'image_file[]'
 
  double step=double(range/N);

  int counter=0;
  FILE *img;
  img=fopen(image_file,"w");	//This image file can be either the reconstructed image or its fourier transform
  
  fprintf(stderr, "\nCreating %s\n", image_file);
  fprintf(stderr, "step = %0.3f\nmaximum range = %0.3f \n", step, step*n/2);
  
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){	
      double x=i*step-n/2*step;
      double y=j*step-n/2*step;
      if(i+offset<n && j+offset<n){
	fprintf(img, "%lf  %lf  %lf\n", x, y, image[i+offset][j+offset]);
	if(image[i+offset][j+offset]!=0) counter++;	    
      }
      else fprintf(img, "%lf  %lf  %lf\n", x, y, 0.0);

    }
    fprintf(img, "\n");//This extra space is only for gnuplot purposes. 
  }	
  fclose(img);
  if(counter==0) fprintf(stderr, "Empty image!\n");
  if(counter!=0) fprintf(stderr, "Non empty image\n");
  
}

//=====================================================================
//=====================================================================

double **fftw_matrix(int n, fftw_complex *a, double **matrix, int flag){

  for (int i = 0; i < n; i++ ){
    for (int j = 0; j < n; j++ ){
      if(flag==0) matrix[i][j]=double(a[i*n+j][0]);//Real part
      if(flag==1) matrix[i][j]=double(a[i*n+j][1]);//Imaginary part
      if(isnan(double(a[i*n+j][0])) || isnan(double(a[i*n+j][1]))){
	fprintf(stderr, "NAN IMAGE (fftw_matrix)\n");
	exit(1);
      }
    }
  }
  return matrix;
}

//====================================================================
//====================================================================


fftw_complex *fftw_format(int n, double **real, double **imaginary, fftw_complex *a){

  //For a 2D complex NX by NY array used by FFTW, we need to access elements   as follows:
  // a[i*ny+j][0] is the real      part of A(I,J).
  // a[i*ny+j][1] is the imaginary part of A(I,J)..

  a = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * n * n );

  for (int i = 0; i < n; i++ ){
    for (int j = 0; j < n; j++ ){
      a[i*n+j][0] = double(real[i][j]);
      a[i*n+j][1] = double(imaginary[i][j]);
    }
  }

  return a;

}

//====================================================================
//====================================================================
/*
double **allocate_matrix(int i_max, int j_max, double **f ){
  //Allocates memory for a matrix
  //double **f=new double*[i_max];
  f=new double*[i_max];
  for(int j=0; j<j_max; j++){
    f[j]=new double[j_max];
  }
  return f;
}
*/
//=====================================================================
//=====================================================================

double **truncate_fft(int n1, int n2, double **image1, double **image2){
  //This function takes an image1 (fft) centered at 0,0 of dimension n1xn1 and returns an image2 of dimension n2xn2

  for(int i=0; i<n2/2; i++){
    for(int j=0; j<n2/2; j++){
      image2[i][j]=image1[i][j]; //First quadrant
    }   
  }  
  for(int i=0; i<n2/2; i++){
    for(int j=0; j<n2/2; j++){
      image2[i+n2/2][j]=image1[n1-n2/2+i][j]; //Second quadrant
    }
  }
  for(int i=0; i<n2/2; i++){
    for(int j=0; j<n2/2; j++){
      image2[i][j+n2/2]=image1[i][n1-n2/2+j]; //Third quadrant
    }
  }
  for(int i=0; i<n2/2; i++){
    for(int j=0; j<n2/2; j++){
      image2[i+n2/2][j+n2/2]=image1[n1-n2/2+i][n1-n2/2+j]; //fourth quadrant
    }
  }


  return image2;

}

//=====================================================================
//=====================================================================

double **read_pgm(int n, double fov, char *file_name){

  fprintf(stderr, "#Reading file: %s\n", file_name);
    
    int res=N;
    double width=fov; //
    
    FILE *pgm;
    
    int i, j;
    
    //char *crap1[10]; //FIX THIS
    char c1;
    char c2;
    int npix=0;
    int crap3=0;
    int crap4=0;
    int flag=0;
    pgm=fopen(file_name,"r");   
    //fscanf(pgm, "%s\n%d %d\n%d\n", &crap1, &npix, &crap3, &crap4);
    fscanf(pgm, "%c%c\n%d %d\n%d\n", &c1, &c2, &npix, &crap3, &crap4);
    
    
    fprintf(stderr, "Pix number %d\n", npix );
    double **image=allocate_matrix(npix, npix, image);
    
    for(j=0; j<npix; j++){
      for(i=0; i<npix; i++){
	fscanf(pgm, "%lf", &image[i][j]);
	if (image[i][j]<5) image[i][j]=0; //TODO: REMOVE THIS!
	if(image[i][j]!=0 && flag==0){ fprintf(stderr, "Non zero image\n"); flag=1;}
      }
    }
    
    fclose(pgm);
    
    return image;
}


//=====================================================================
//=====================================================================

void read_pupil(struct ppl *ppl){
  //This function reads a pupil parameter file (.par) and saves aperture positions and diameters.
  fprintf(stderr, "The pupil file name is %s\n", ppl->file_name);
  FILE *pupil_file=fopen(ppl->file_name, "r");
  int i=0;
  //fscanf(pupil_file, "%d\n%lf\n%lf\n", &ppl->n_sub, &ppl->size, &ppl->step);
  fscanf(pupil_file, "%d\n%lf\n", &ppl->n_sub, &ppl->size);
  
  fprintf(stderr, "Number of sub-apertures is: %d\n", ppl->n_sub);
  //fprintf(stderr, "Physical pupil size is: %0.3f m\n", ppl->size); 
  ppl->pos_x=new int[ppl->n_sub];
  ppl->pos_y=new int[ppl->n_sub];
  ppl->diameter=new int[ppl->n_sub];
  ppl->diameter_c=new int[ppl->n_sub];

  while(fscanf(pupil_file, "%d %d %d\n", &ppl->pos_x[i], &ppl->pos_y[i], &ppl->diameter[i])!=EOF ) {
    ppl->pos_x[i]=(int)floor(ppl->pos_x[i]*ppl->pscale+0.5);
    ppl->pos_y[i]=(int)floor(ppl->pos_y[i]*ppl->pscale+0.5);
    ppl->diameter[i]=(int)floor(ppl->diameter[i]*ppl->pscale+0.5);
    if(ppl->pos_x[i]>N || ppl->pos_x[i]>N){
      fprintf(stderr, "Pupil scaling error!\n");
      exit(1);
    }
    ppl->diameter_c[i]=2; 
    i++;
 }
  fclose(pupil_file);
  fprintf(stderr, "Done reading pupil\n");
}


//=====================================================================
//=====================================================================

void main_spot_find(int n, double threshold, struct ppl *ppl, double **image){
  //This function receives a pgm file, finds the sub-aperture positions and saves them in the ppl structure
  fprintf(stderr, "Finding sub-pupil positions\n");

  struct spot_params params; //defined in spot_find.c
  struct image img; //defined in spot_find.c
  
  params.brightness=0;

  img.pixel=new float[n*n];
  img.x=new float[n*n];
  img.y=new float[n*n];

  ppl->pos_x=new int[n*n];
  ppl->pos_y=new int[n*n];
  ppl->diameter=new int[n*n];

  int j=0;  
  for(int l=0; l<n; l++){
    for(int k=0; k<n; k++){
      img.pixel[j]=image[k][l];
      img.x[j]=k;
      img.y[j]=l;
      j++;
    }
  }

  int i=1;  
  while(params.brightness>=0){ //Brightness will be set to -1 when it stops finding spots    
    img.pixel=find_spot(img.pixel, n*n, threshold, -(double)i, &params, &img); //defined in spot_find.c
    if(params.radius>TINY){
      ppl->pos_x[i-1]=params.xcm;
      ppl->pos_y[i-1]=params.ycm;
      ppl->diameter[i-1]=(int)float(params.radius*2+0.5);
      //ppl->diameter_c[i]=2;
      i++;
    }
  }
  ppl->n_sub=i-1;
  fprintf(stderr, "Done finding sub-apertures\n%d sub-aperutes\n", ppl->n_sub);
}

//=====================================================================
//=====================================================================


int readpilot(struct plt *plt) {
  //TODO: check for errors, missing parameters or conflicts
  //Reads the pilot file
  fprintf(stderr, "The pilot file name is %s\n",plt->pilotfile);
  
  FILE *fpilot;           //Pilot file unit  
  //Opens the pilot file
  if((fpilot=fopen(plt->pilotfile,"r"))==NULL){
    fprintf(stderr, "Failed opening the pilot file\n");  
    exit(1);
  }

  int ieof=0;
  char flag[FLAG_LENGTH];
  while(!ieof){  
    /*we look for a flag*/
    look_for_flag(&ieof,flag,fpilot);
    if (ieof) break; 
        
    //Pupil .par file name
    if (strstr(flag,"PPLPAR")!=NULL) {
      get_word(fpilot, plt->pupil_par_name, &ieof);
    }

    //Error in sub-pupil location
    if (strstr(flag,"PPLERR")!=NULL) {
      plt->ppl_uncertain=1;
      fscanf(fpilot,"%lf",&(plt->ppl_err));
    }

    //Piston error in each sub-pupil
    if (strstr(flag,"PPLPHS")!=NULL) {
      plt->ppl_phase=1;
      fscanf(fpilot,"%lf",&(plt->ppl_phase_error));
    }

    //Pupil .par file generated from fits or pgm file of pupil
    if (strstr(flag,"GENPAR")!=NULL) {
      plt->gen_par_flag=1;
      get_word(fpilot, plt->gen_par_name, &ieof);
    }

    //Pupil fits file name
    if (strstr(flag,"PUPILN")!=NULL) {
      get_word(fpilot, plt->pupil_fits_name, &ieof);
    }

    //Pristine fits image name
    if (strstr(flag,"PRIFTS")!=NULL) {
      plt->fits_flag=1;
      get_word(fpilot, plt->pristine_name, &ieof);
      fscanf(fpilot, "%d", &plt->fits_npix);
    }

    //Pristine pgm image name
    if (strstr(flag,"PRIPGM")!=NULL) {
      plt->pgm_flag=1;
      get_word(fpilot, plt->pristine_name, &ieof);
      fscanf(fpilot, "%d", &plt->pgm_npix);
    }


    //Resampled pristine image name
    if (strstr(flag,"RPRIST")!=NULL) {
      get_word(fpilot, plt->res_pristine_name, &ieof);
    }

    //fft magnitude fits name
    if (strstr(flag,"FFTMAG")!=NULL) {
      get_word(fpilot, plt->mod_fits_name, &ieof);
    }

    //Pupil pgm file name
    if (strstr(flag,"PPLPGM")!=NULL) {
      get_word(fpilot, plt->pupil_pgm_name, &ieof);
      plt->pupil_pgm_flag=1;
    }

    //Pupil fits file name
    if (strstr(flag,"PPLFTS")!=NULL) {
      get_word(fpilot, plt->pupil_fits_read, &ieof);
      plt->pupil_fits_flag=1;
    }


    //Fizeau image file name
    if (strstr(flag,"FIZEAU")!=NULL) {
      get_word(fpilot, plt->fizeau_fits_name, &ieof);
    }

    //FFT Fizeau image file name
    if (strstr(flag,"FFTFIZ")!=NULL) {
      get_word(fpilot, plt->fft_fizeau_name, &ieof);
    }

    //FFT Fizeau image file name
    if (strstr(flag,"FFTDEN")!=NULL) {
      get_word(fpilot, plt->fft_densified_name, &ieof);
    }

    //mtf image file name
    if (strstr(flag,"MODTRF")!=NULL) {
      get_word(fpilot, plt->mtf_fits_name, &ieof);
    }

    //Densified image file name
    if (strstr(flag,"MICHEL")!=NULL) {
      get_word(fpilot, plt->densified_fits_name, &ieof);
    }
    
    //Pristine image physical size in mas
    if (strstr(flag,"IMGSIZ")!=NULL) {
      fscanf(fpilot,"%lf",&(plt->image_size));
    }

    //Pupil size in meters
    if (strstr(flag,"PPLSIZ")!=NULL) {
      fscanf(fpilot,"%lf",&(plt->pupil_size));
    }

    //Wavelength in nm
    if (strstr(flag,"LAMBDA")!=NULL) {
      fscanf(fpilot,"%lf",&(plt->lambda));
    }

    //Densification factor
    if (strstr(flag,"DFACTR")!=NULL) {
      fscanf(fpilot,"%lf",&(plt->gamma));
    }

    //Densification factor correction
    if (strstr(flag,"DFCORR")!=NULL) {
      plt->dfcorr=1;
    }  

    //Average number of photons and thermal photons per pixel
    if (strstr(flag,"PHOTON")!=NULL) {
      fscanf(fpilot,"%lf %lf",&(plt->nphoton), &(plt->tphoton));
      plt->photon=1;
    }

    //Average number of photons per pixel
    if (strstr(flag,"PSFOUT")!=NULL) {
      plt->psf_flag=1;
    }
    
    
  }
  
  //Closes the pilot file
  fclose (fpilot);
  
  //TODO: POSSIBLE ERROR MESSAGES use flags
  fprintf(stderr, "Done reading the pilot file\n");
  return 0;
}

//============================================
//============================================


void initialize_pilot(struct plt *plt){
  //TODO: create error messages in case some options are not provided or conflicting options
  plt->lambda=550E-9;
  plt->photon=0;
  plt->tphoton=0;
  plt->pupil_pgm_flag=0;
  plt->pupil_fits_flag=0;
  plt->gen_par_flag=0;
  plt->fits_flag=0;
  plt->pgm_flag=0;
  plt->pgm_npix=0;
  plt->ppl_uncertain=0;
  plt->ppl_phase=0;
  plt->ppl_phase_error=0;
  plt->ppl_err=0;
  plt->psf_flag=0;
  plt->dfcorr=0;
  plt->pupil_size=-1;
  plt->pilotfile=new  char[STRING_LENGTH];  
  plt->pristine_name=new char[STRING_LENGTH];  
  plt->res_pristine_name=new char[STRING_LENGTH];  
  plt->mod_fits_name=new char[STRING_LENGTH];  
  plt->mtf_fits_name=new char[STRING_LENGTH];  
  plt->pupil_par_name=new char[STRING_LENGTH];
  plt->pupil_pgm_name=new char[STRING_LENGTH];
  plt->gen_par_name=new char[STRING_LENGTH];
  plt->pupil_fits_name=new char[STRING_LENGTH];
  plt->pupil_fits_read=new char[STRING_LENGTH];
  plt->fizeau_fits_name=new char[STRING_LENGTH];
  plt->fft_fizeau_name=new char[STRING_LENGTH];
  plt->fft_densified_name=new char[STRING_LENGTH];
  plt->densified_fits_name=new char[STRING_LENGTH];

}

//============================================
//============================================

double **initialize_matrix(int n, double **image, double init){
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      image[i][j]=init;
    }
  }
  return image;
}

//=====================================================================
//=====================================================================

double **random_phase(int n, double sigma, double **phase, struct ppl *ppl){
  fprintf(stderr, "Including phase errors\n");
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  long seed;
  seed = time (NULL);    
  gsl_rng_set (r, seed);       

  phase=initialize_matrix(n, phase, 0);
  for(int i=0; i<ppl->n_sub; i++){
    double random=gsl_ran_gaussian (r, sigma);
    for(int l=ppl->pos_y[i]-int(ppl->diameter[i]/2); l<=ppl->pos_y[i]+int(ppl->diameter[i]/2); l++){
      for(int k=ppl->pos_x[i]-int(ppl->diameter[i]/2); k<=ppl->pos_x[i]+int(ppl->diameter[i]/2); k++){
	if((k-ppl->pos_x[i])*(k-ppl->pos_x[i])+(l-ppl->pos_y[i])*(l-ppl->pos_y[i])<=int(ppl->diameter[i]/2)*int(ppl->diameter[i]/2)) phase[k][l]=random; //Phased within sub-pupils
      }
    }
  }

  /*for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(ppl->pupil[i][j]>0){
	phase[i][j]=gsl_ran_gaussian (r, sigma);
	//printf("%f  ", phase[i][j]);
      }
      }*/


  return phase;

}

//=====================================================================
//=====================================================================


void ppl_uncertainty(double sigma, struct ppl *ppl){

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  for(int i=0; i<ppl->n_sub; i++){
    ppl->pos_x[i]=ppl->pos_x[i]+gsl_ran_gaussian (r, sigma);
    ppl->pos_y[i]=ppl->pos_y[i]+gsl_ran_gaussian (r, sigma);
  }
 

}

//=====================================================================
//=====================================================================
