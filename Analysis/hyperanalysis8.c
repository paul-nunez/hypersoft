
// Compile as: g++ hyperanalysis3.c -o hyperanalysis -lm -lcfitsio -lgsl -lgslcblas -lfftw3

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

# include <fftw3.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#define N 1024
#define STRING_LENGTH 100
#define FLAG_LENGTH 6
#define PI 3.141592
#define TINY 1e-30
#define SMALL 1e-5
#define RAD_TO_MAS 206264849.15927

#include "../Tools/fits.c" //Includes functions to generate fits files
#include "../Tools/misc.c"
#include "../Tools/spot_find.c"

double **generate_mask(int n, int mask_rad, double **mask);
void checknan(int n, double **image);
//double **allocate_matrix(int i_max, int j_max, double **image); //In ../Tools/misc.c
double **initialize_matrix(int n, double **image, double init);
double **scale_image(int n, double desired_size, double input_size, double **prist, double **scaled);
double **read_pgm(int n, double fov, char *file_name);
fftw_complex *fftw_format(int n, double **real, double **imaginary, fftw_complex *a);
double **fftw_matrix(int n, fftw_complex *a, double **matrix, int flag);
fftw_complex *normalize_fftw(int n, fftw_complex *a);
double **transpose(int n, double **image);
//double **change_to_2d(int n, double *image, double **data, int q); //In ../Tools/misc.c
void print_file(int n, int offset, double step,  double **image, char image_file[]);
double **modulus(int n, double **re, double **im);
double **shift(int n, double **image);
double **multiply(int n, double **image1, double **image2);
fftw_complex *autocorrelate_fftw(int n, double **image, double **phase);
void display(char *name, int n, int offset, double range, double **image);
double **normalize(int n, double **image);
double **remove_offset(int n, double **image, double offset);
fftw_complex *conv_fftw(int n, fftw_complex *img, fftw_complex *mtf, fftw_complex *conv);
double **generate_pupil(int n, struct ppl *ppl);
void read_pupil(struct ppl *ppl);
void main_spot_find(int n, double threshold, struct ppl *ppl, double **image);
fftw_complex *rich_lucy(int n, int n_iter, fftw_complex *image, fftw_complex *mtf, fftw_complex *mask, fftw_complex *prist, char *error_name);
double **densify(int n, double gamma, struct mtfunc *mtfunc, double **image); //TODO: STILL ISSUES HERE?
double **dilute(int n, double gamma, struct mtfunc *mtfunc, double **image);
void redefine_mtf(int n, int flag, double gamma, struct mtfunc *mtfunc);
double **poisson(int n, double **image);
double **include_photons(int n, double avgpp, double **image);
double find_ppl_spacing(int n, struct ppl *ppl);
void find_mtf_centers(int n, struct ppl *ppl, struct mtfunc *mtf);
double find_error(int n, fftw_complex *conv_img, fftw_complex *conv_guess);
double *correlation(int n, fftw_complex *image1, fftw_complex *image2);
double find_error_matrix(int n, double **A, double **B); //TODO: NOT USED?
int readpilot(struct plt *plt); //in Tools/pilot_read.c
void initialize_pilot(struct plt *plt);
double **read_file(int n, int offset, double range,  double **image, char image_file[]);

struct plt{
  int prist_flag;
  int psf_flag;
  char *pilotfile;
  double lambda;
  char *pristine_name;
  char *pupil_pgm_name;
  char *pgm_par_name;
  char *gen_par_name;
  int gen_par_flag;
  double pupil_size;
  int pupil_pgm_flag;
  char *pupil_fits_read;
  int pupil_fits_flag;
  char *pupil_par_name;
  char *pupil_fits_name;
  char *mod_fits_name;
  char *mtf_fits_name;
  double tphoton;
  double nphoton;
  int dfcorr;
  double gamma;
  char *fft_fizeau_name;
  char *fizeau_fits_name;
  char *fft_densified_name;
  char *densified_fits_name;
  char *data_fits_name;
  char *de_densified_fits_name;
  char *rl_fits_name;
  char *error_name;
  double image_size;
  int nrl;
  int mask_flag;
  double mask;
};

struct ppl{
  char *file_name;
  int n_sub;
  int *pos_x;
  int *pos_y;
  double step;
  double size;
  double pscale;
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
  int n_sub;
};

int main(int n, char **arg){

  struct plt plt;
  struct ppl ppl;
  struct mtfunc mtfunc;

  if(arg[1]==NULL || strcmp(arg[1], "-h")==0 || n<1){
        fprintf(stderr, "\nPlease use a pilot file with the -p flag\n\n");    
	fprintf(stderr, "Example:\n ./hyperanalysis -p analysis.pilot\n");
	exit(1);
  }  


  initialize_pilot(&plt); //Initialize pilot variables

  if (strcmp(arg[1], "-p")==0){
    strcpy(plt.pilotfile, arg[2]);
    fprintf(stderr, "Pilot file name: %s\n", plt.pilotfile);
    readpilot(&plt);
  }
  else{
    fprintf(stderr, "Pilot file not given?\n");
    exit(1);
  }

  if(plt.pupil_pgm_flag==1) fprintf(stderr, "Pupil .pgm file provided\n");
  if(plt.pupil_pgm_flag==0) fprintf(stderr, "Pupil .par file provided\n");

  //Maximal size of pupil plane can be defined by the pristine image pixel size and lambda.
  //The pupil size is defined by the user. Sub-aperture pixel positons must be scaled 
  //by the ratio of the maximal pupil size and the user-specified pupil size.
  
  //Allocate memory for pupil
  ppl.pupil=allocate_matrix(N, N, ppl.pupil);
  ppl.pupil_centers=allocate_matrix(N, N, ppl.pupil_centers);

  //Set physical pupil size if it is supplied in the pilot file
  if(plt.pupil_size!=-1) ppl.size=plt.pupil_size;
  ppl.step=ppl.size/N;  //Pupil pixel size
  fprintf(stderr, "Image size %lf mas\nPupil size %lf \nPupil pixel size %lf\n", plt.image_size, plt.pupil_size, ppl.step);
  plt.lambda=plt.lambda*1E-9; //nanometers

  //CURRENTLY NOT SCALING PUPIL. INTRODUCES PROBLEMS IN DENSIFICATION
  double maximal_pupil_size=(plt.lambda*(N/plt.image_size)*RAD_TO_MAS);
  ppl.pscale=1; 
  fprintf(stderr, "Maximal pupil size is : %lf\n", maximal_pupil_size);

  //Read pupil FITS file if supplied in pilot file
  if(plt.pupil_fits_flag==1){
    double *ppl1=new double[2*N*N];
    ppl1=readfits(plt.pupil_fits_read, N);
    fprintf(stderr, "Fits read\n");
    ppl.pupil=change_to_2d(N, ppl1, ppl.pupil, 0);
    ppl.pupil=normalize(N, ppl.pupil);
  }
  
  //Read pupil PGM file if supplied in pilot file
  if(plt.pupil_pgm_flag==1){
    ppl.pupil=read_pgm(N, ppl.size, plt.pupil_pgm_name);
    ppl.pupil=normalize(N, ppl.pupil);   
    //Transpose pupil
    ppl.pupil=transpose(N, ppl.pupil); //Only necessary to be consistent with David's images
  }
  
  //Read pupil parameter file unless a pupil pgm/fits file is given
  ppl.file_name=new char[STRING_LENGTH];
  strcpy(ppl.file_name, plt.pupil_par_name); 
  if(plt.pupil_pgm_flag==0 && plt.pupil_fits_flag==0){
    read_pupil(&ppl); //This can also scale pupil. TODO: when scaling pupil, would be nice to center it as well
    ppl.pupil=transpose(N, ppl.pupil); //Only necessary to be consistent with David's images
  }
  if(plt.pupil_size!=-1){
    ppl.size=plt.pupil_size;
    ppl.step=ppl.size/N; 
  }

  //Generate pupil from parameter .par file unless a pgm file is given
  if(plt.pupil_pgm_flag==0 && plt.pupil_fits_flag==0) ppl.pupil=generate_pupil(N, &ppl);
  
  //Create pupil fits file
  fitsout(N, ppl.step, ppl.pupil,  plt.pupil_fits_name); 

  double image_step=plt.image_size/N;

  //find mtf centers with pupil information
  find_mtf_centers(N, &ppl, &mtfunc);

  //Create OTF from input pupil
  double *otf1=new double[2*N*N];
  double **re_otf=allocate_matrix(N, N, re_otf);
  double **im_otf=allocate_matrix(N, N, im_otf);  
  double **pupil_phase=allocate_matrix(N, N, pupil_phase);
  pupil_phase=initialize_matrix(N, pupil_phase, 0);  
  fprintf(stderr, "Autocorrelating pupil\n");
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

  //For proper densification: find corrected gamma, we must find minumum mtf spacing
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

  fprintf(stdout, "%0.3lf  ", gamma);

  //Up to now, we should have read pupil, scaled it according to image size, found mtf positions, found otf.
  //Now we need to read fits data, re-dilute (starting from re-defined mtf positions), and deconvolve

  double **new_re_otf=allocate_matrix(N, N, new_re_otf);
  double **new_im_otf=allocate_matrix(N, N, new_im_otf);
  new_re_otf=densify(N, gamma, &mtfunc, new_re_otf); 
  new_re_otf=densify(N, gamma, &mtfunc, new_im_otf); 
  redefine_mtf(N, 1, gamma, &mtfunc); //Redifine positions of spots in mtf

  //Read data from fits file
  char *input_name=new char[50];
  strcpy(input_name, plt.data_fits_name);
  double *sim_data_read=new double[2*N*N];
  sim_data_read=readfits(input_name, N); 
  double **sim_data=allocate_matrix(N, N, sim_data);
  sim_data=change_to_2d(N, sim_data_read, sim_data, 0);
  sim_data=transpose(N, sim_data); //Only necessary to be consistent with David's images
  checknan(N, sim_data);

  //*******DATA ANALYSIS STARTS HERE********
  //Substract constant (thermal) background if user desires
  if(plt.tphoton>0) sim_data=remove_offset(N, sim_data, plt.tphoton);  

  //Create data fits file
  char *data_fits_check=new char[STRING_LENGTH];
  strcpy(data_fits_check, "./Images/Secondary/input_data.fits"); //READ FROM PILOT
  fitsout(N, 1, sim_data,  data_fits_check); 

  //Following part is to re-dilute
  sim_data=shift(N, sim_data);
  //Take fftwn
  double **im_data=allocate_matrix(N, N, im_data);
  im_data=initialize_matrix(N, im_data, 0);
  fftw_complex *data=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  data=fftw_format(N, sim_data, im_data, data);//Creates input for fftw
  fftw_complex *data_fftw = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  fftw_plan plan_forward = fftw_plan_dft_2d ( N, N, data, data_fftw, FFTW_FORWARD,  FFTW_ESTIMATE );
  fftw_execute ( plan_forward );  
  double **re_sim_data=allocate_matrix(N, N, re_sim_data);
  re_sim_data=fftw_matrix(N, data_fftw, re_sim_data, 0);
  double **im_sim_data=allocate_matrix(N, N, im_sim_data);
  im_sim_data=fftw_matrix(N, data_fftw, im_sim_data, 1);
  re_sim_data=shift(N, re_sim_data);
  im_sim_data=shift(N, im_sim_data);
  double **re_sim_data2=allocate_matrix(N, N, re_sim_data2);
  double **im_sim_data2=allocate_matrix(N, N, im_sim_data2);
  double inv_gamma=1/gamma;
  
  //re_sim_data2=densify(N, inv_gamma, &mtfunc, re_sim_data); 
  re_sim_data2=dilute(N, inv_gamma, &mtfunc, re_sim_data); 
  //im_sim_data2=densify(N, inv_gamma, &mtfunc, im_sim_data); 
  im_sim_data2=dilute(N, inv_gamma, &mtfunc, im_sim_data); 
  re_sim_data2=shift(N, re_sim_data2);
  im_sim_data2=shift(N, im_sim_data2);

  //Inverse fftw
  fftw_complex *fftw_data=fftw_format(N, re_sim_data2, im_sim_data2, fftw_data);//Creates input for fftw
  fftw_complex *diluted_data = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  fftw_plan plan_backward = fftw_plan_dft_2d ( N, N, fftw_data, diluted_data, FFTW_BACKWARD,  FFTW_ESTIMATE );
  fftw_execute ( plan_backward );  
  sim_data=fftw_matrix(N, diluted_data, sim_data, 0);
  sim_data=shift(N, sim_data);
  sim_data=normalize(N, sim_data); //This is the rediluted data
  
  fitsout(N, 1, sim_data,  plt.de_densified_fits_name);  //TODO :units  

  sim_data=shift(N, sim_data); 
  mtf_diluted=shift(N, mtf_diluted);  
  
  //Deconvolve with Richardson-Lucy
  double **zero=allocate_matrix(N, N, zero);//Imaginary part of data is zero
  zero=initialize_matrix(N, zero, 0);
  fftw_complex *sim_data_fftw=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  fftw_complex *mtf_diluted_fftw=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  mtf_diluted_fftw=fftw_format(N, mtf_diluted, mtf_diluted, mtf_diluted_fftw);
  sim_data_fftw=fftw_format(N, sim_data, zero, sim_data_fftw);

  //In case pristine image is provided for comparison with deconvolved image
  double **prist=allocate_matrix(N, N, prist);
  fftw_complex *prist_fftw=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  if(plt.prist_flag==1){
    double *prist1=new double[2*N*N];
    prist1=readfits(plt.pristine_name, N);
    prist=change_to_2d(N, prist1, prist, 0);
    prist=normalize(N, prist);
    prist=transpose(N, prist);
    prist=shift(N, prist);
    //fitsout(N, 1, prist, "prist.fits");
    prist_fftw=fftw_format(N, prist, zero, prist_fftw);
  }
  if(plt.prist_flag==0) prist_fftw=fftw_format(N, zero, zero, prist_fftw);

  int mask_rad=plt.mask*N/plt.image_size;  
  double **mask=allocate_matrix(N, N, mask);
  if(plt.mask_flag==1) mask=generate_mask(N, mask_rad, mask);
  else mask=initialize_matrix(N, mask, 1);
  mask=shift(N, mask);
  fftw_complex *mask_fftw=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * N * N );
  mask_fftw=fftw_format(N, mask, zero, mask_fftw);
  //Here I call Richardson-Lucy alg.
  sim_data_fftw=rich_lucy(N, plt.nrl, sim_data_fftw, mtf_diluted_fftw, mask_fftw, prist_fftw, plt.error_name); 
  sim_data=fftw_matrix(N, sim_data_fftw, sim_data, 0);
  sim_data=shift(N, sim_data);

  //sim_data=transpose(N, sim_data);  
  fitsout(N, image_step, sim_data,  plt.rl_fits_name); 
  
  return 0;
  
}

//=====================================================================
//=====================================================================

double find_error_matrix(int n, double **A, double **B){

  double sum=0;
  double sum1=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      sum=sum+(A[i][j]-B[i][j])*(A[i][j]-B[i][j]);
      sum1=sum1+(A[i][j]*B[i][j]);      
    }
  }
  double error=0;
  error=sum/sum1;
  if(sum1==0){fprintf(stderr, "Division by 0 in find_data_error\n"); exit(1);}

  return error;
}

//=====================================================================
//=====================================================================

double **generate_mask(int n, int mask_rad, double **mask){

  mask=initialize_matrix(n, mask, 1);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(i>mask_rad/2 && i<n-mask_rad/2) mask[i][j]=0; //Apply mask for faster convergence	
      if(j>mask_rad/2 && j<n-mask_rad/2) mask[i][j]=0; //Apply mask for faster convergence	      
      //if(pow(i-n/2,2)+pow(j-n/2,2)<pow(mask_rad,2)) mask[i][j]=1;
      //else mask[i][j]=0;
    }
  }
  mask=shift(n, mask);
  return mask;

}

//=====================================================================
//=====================================================================

double find_error(int n, fftw_complex *ref, fftw_complex *conv_guess){

  double sum=0;
  double sum1=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      sum=sum+(ref[i*n+j][0]-conv_guess[i*n+j][0])*(ref[i*n+j][0]-conv_guess[i*n+j][0]);
      sum1=sum1+(ref[i*n+j][0]*ref[i*n+j][0]);      
      //sum=sum+(ref[i*n+j][0]-conv_guess[j*n+i][0])*(ref[i*n+j][0]-conv_guess[j*n+i][0]);
      //sum1=sum1+(ref[i*n+j][0]*ref[i*n+j][0]);      
    }
  }
  double error=0;
  error=sum/sum1;
  if(sum1==0){fprintf(stderr, "Division by 0 in find_data_error\n"); exit(1);}

  return error;
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
  fftw_complex *mtfc=(fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * n * n );
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

double **include_photons(int n, double avgpp, double **image){ //TODO: ERASE THIS?
  fprintf(stderr, "Including photons per pixel");
  //First find average of image
  double sum=0;
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      sum=sum+image[i][j];
    }
  }
  double average=sum/(n*n);

  //Now include number of photons per pixel without noise
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      image[i][j]=image[i][j]*avgpp/average;
    }
  }

  //Now add Poisson noise
  image=poisson(N, image); 

  return image;

}

//=====================================================================
//=====================================================================


double **poisson(int n, double **image){
  //Taken from http://www.gnu.org/software/gsl/manual/html_node/Random-Number-Distribution-Examples.html

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
     
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

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
  fprintf(stderr, "Scaling image by %0.3f\n", s);
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

double **densify(int n, double gamma, struct mtfunc *mtfunc, double **image){
  //TODO:  must deal with overlapping frequencies
  fprintf(stderr, "Densifying\n");

  double **image_d=allocate_matrix(n, n, image_d);
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      image_d[i][j]=0; //Initialize
    }
  }
  
  mtfunc->dpos_x=new int[mtfunc->n_sub];
  mtfunc->dpos_y=new int[mtfunc->n_sub];

  //gamma=find_next_gamma(n, gamma, mtfunc);

  for(int k=0; k<mtfunc->n_sub; k++){
    //double a=1.1;
    double a=1.25;
    double d=sqrt((mtfunc->pos_x[k]-n/2)*(mtfunc->pos_x[k]-n/2)+(mtfunc->pos_y[k]-n/2)*(mtfunc->pos_y[k]-n/2));
    int dx=int(floor(mtfunc->pos_x[k]-n/2+0.5));
    int dy=int(floor(mtfunc->pos_y[k]-n/2+0.5));

    //mtfunc->dpos_x[k]=int(floor((mtfunc->pos_x[k]-n/2)-(gamma*(mtfunc->pos_x[k]-n/2)+0.5)))+1; //TODO:UNDERSTAND THIS EXTRA PIXEL
    mtfunc->dpos_x[k]=int(floor((mtfunc->pos_x[k]-n/2)-(gamma*(mtfunc->pos_x[k]-n/2))+0.5));
    //mtfunc->dpos_y[k]=int(floor((mtfunc->pos_y[k]-n/2)-(gamma*(mtfunc->pos_y[k]-n/2)+0.5)))+1;
    mtfunc->dpos_y[k]=int(floor((mtfunc->pos_y[k]-n/2)-(gamma*(mtfunc->pos_y[k]-n/2))+0.5));

    int start_y=int(floor(mtfunc->pos_y[k]-a*mtfunc->radius[k]+0.5));
    int start_x=int(floor(mtfunc->pos_x[k]-a*mtfunc->radius[k]+0.5));
    int end_y=int(floor(mtfunc->pos_y[k]+a*mtfunc->radius[k]+0.5));
    int end_x=int(floor(mtfunc->pos_x[k]+a*mtfunc->radius[k]+0.5));
    for(int j=start_y; j<end_y; j++){
      for(int i=start_x; i<end_x; i++){
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


double **dilute(int n, double gamma, struct mtfunc *mtfunc, double **image){
  //This function should be used after knowing the frequency displacements, which are found with the densify function.
  fprintf(stderr, "Diluting\n");

  double **image_d=allocate_matrix(n, n, image_d);
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      image_d[i][j]=0; //Initialize
    }
  }

  double a=1.25;
  for(int k=0; k<mtfunc->n_sub; k++){
    int start_y=int(floor(mtfunc->pos_y[k]-a*mtfunc->radius[k]+0.5));
    int start_x=int(floor(mtfunc->pos_x[k]-a*mtfunc->radius[k]+0.5));
    int end_y=int(floor(mtfunc->pos_y[k]+a*mtfunc->radius[k]+0.5));
    int end_x=int(floor(mtfunc->pos_x[k]+a*mtfunc->radius[k]+0.5));
    for(int j=start_y; j<end_y; j++){
      for(int i=start_x; i<end_x; i++){
	if(i>0 && i<n && j>0 &&j<n){
	  if(i+mtfunc->dpos_x[k]>=0 && i+mtfunc->dpos_x[k]<n && j+mtfunc->dpos_y[k]>=0 && j+mtfunc->dpos_y[k]<n){
	    image_d[i+mtfunc->dpos_x[k]][j+mtfunc->dpos_y[k]]=image[i][j]; 
	  }
	}
      }
    }

  }
  return image_d;
  
}

//=====================================================================
//=====================================================================


double **densify_old(int n, double gamma, struct mtfunc *mtfunc, double **image){
  //TODO:  must deal with overlapping frequencies
  fprintf(stderr, "Densifying\n");
  
  double **image_d=allocate_matrix(n, n, image_d);
  double **image_old=allocate_matrix(n, n, image_d);
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      image_d[i][j]=0; //Initialize
      image_old[i][j]=0; //Initialize
    }
  }
  
  mtfunc->dpos_x=new int[mtfunc->n_sub];
  mtfunc->dpos_y=new int[mtfunc->n_sub];

    //gamma=gamma-1;//Usual definition of densification factor
  for(int k=0; k<mtfunc->n_sub; k++){
    //double a=1.0;
    double a=1.2;
    double d=sqrt((mtfunc->pos_x[k]-n/2)*(mtfunc->pos_x[k]-n/2)+(mtfunc->pos_y[k]-n/2)*(mtfunc->pos_y[k]-n/2));
    int dx=int(floor(mtfunc->pos_x[k]-n/2+0.5));
    int dy=int(floor(mtfunc->pos_y[k]-n/2+0.5));
    double alpha=(mtfunc->pos_y[k]-n/2)/(mtfunc->pos_x[k]-n/2);//Conserved quantity
    if(abs(dx)>0){ 
      mtfunc->dpos_x[k]=int(floor(((mtfunc->pos_x[k]-n/2)/abs(mtfunc->pos_x[k]-n/2))*(gamma*d/sqrt(1+alpha*alpha))+0.5));
      mtfunc->dpos_y[k]=int(floor(((mtfunc->pos_y[k]-n/2)/abs(mtfunc->pos_y[k]-n/2))*(gamma*sqrt(alpha*alpha)*d/sqrt(1+alpha*alpha))+0.5));
    }
    if(dy==0){
      mtfunc->dpos_x[k]=int(floor(((mtfunc->pos_x[k]-n/2)/abs(mtfunc->pos_x[k]-n/2))*(gamma*d/sqrt(1+alpha*alpha))+0.5));
      mtfunc->dpos_y[k]=int(floor(gamma*sqrt(alpha*alpha)*d/sqrt(1+alpha*alpha)+0.5));          
    }
    if(dx==0){ 
      mtfunc->dpos_x[k]=0;
      mtfunc->dpos_y[k]=int(floor(gamma*d+0.5));
    }
    if(dx==0 && dy<0){ 
      mtfunc->dpos_x[k]=0;
      mtfunc->dpos_y[k]=-int(floor(gamma*d+0.5));
    }
    int start_y=int(floor(mtfunc->pos_y[k]-a*mtfunc->radius[k]+0.5));
    int start_x=int(floor(mtfunc->pos_x[k]-a*mtfunc->radius[k]+0.5));
    int end_y=int(floor(mtfunc->pos_y[k]+a*mtfunc->radius[k]+0.5));
    int end_x=int(floor(mtfunc->pos_x[k]+a*mtfunc->radius[k]+0.5));
    for(int j=start_y; j<end_y; j++){
      for(int i=start_x; i<end_x; i++){
	if(i>0 && i<n && j>0 &&j<n){
	  if(i-mtfunc->dpos_x[k]>0 && i-mtfunc->dpos_x[k]<n && j-mtfunc->dpos_y[k]>0 && j-mtfunc->dpos_y[k]<n){
	    image_d[i-mtfunc->dpos_x[k]][j-mtfunc->dpos_y[k]]=image[i][j]; 
	    //Deals with frequency overlaps?
	    //image_d[i-mtfunc->dpos_x[k]][j-mtfunc->dpos_y[k]]=image[i][j]+image[i-mtfunc->dpos_x[k]][j-mtfunc->dpos_y[k]]; 	    
	  }
	}
      }
    }

  }

  return image_d;
  
}

//=====================================================================
//=====================================================================

fftw_complex *conv_fftw(int n, fftw_complex *img, fftw_complex *mtf, fftw_complex *conv){
  
  fftw_complex *fft_img=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n*n);

  fftw_plan plan_forward = fftw_plan_dft_2d ( n, n, img, fft_img, FFTW_FORWARD,  FFTW_ESTIMATE );//Plan optimal fft algthm  
  fftw_execute(plan_forward); //Take fft

  //Multiply in Fourier space
  for (int i = 0; i < n; i++ ){
    for (int j = 0; j < n; j++ ){
      if(isnan(fft_img[i*n+j][0]) || isnan(fft_img[i*n+j][1])){fprintf(stderr, "NAN IMAGE (conv)\n"); exit(1);}
      if(isnan(mtf[i*n+j][0]) || isnan(mtf[i*n+j][1])){fprintf(stderr, "NAN MTF (conv)\n"); exit(1);}
      img[i*n+j][0]=fft_img[i*n+j][0]*mtf[i*n+j][0];//Real part
      img[i*n+j][1]=fft_img[i*n+j][1]*mtf[i*n+j][1];//Imaginary part
    }
  }

  fftw_plan plan_backward = fftw_plan_dft_2d ( N, N, img, conv, FFTW_BACKWARD,  FFTW_ESTIMATE );//Plan optimal fft algthm
  fftw_execute(plan_backward); //Take inverse fftw
  
  fftw_destroy_plan ( plan_forward );
  fftw_destroy_plan ( plan_backward );

  fftw_free(fft_img);

  return conv;

}


//=====================================================================
//=====================================================================

fftw_complex *rich_lucy(int n, int n_iter, fftw_complex *image, fftw_complex *mtf, fftw_complex *mask, fftw_complex *prist, char *error_name){

  fprintf(stderr, "Perfoming %d Richardson-Lucy iterations\n", n_iter);
  image=normalize_fftw(n, image);

  fftw_complex *image_0=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n*n);

  //Initialize starting image
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      //image_0[i*n+j][0]=image[i*n+j][0]*mask[i*n+j][0]; //This is the convolved data
      image_0[i*n+j][0]=image[i*n+j][0];
      image[i*n+j][1]=0; 
      image_0[i*n+j][1]=0; //Imaginary part is set to 0
      if(isnan(image[i*n+j][0]) || isnan(image[i*n+j][1])){
	fprintf(stderr, "NAN IMAGE\n");
	exit(1);
      }
    }
  }
  image_0=normalize_fftw(n, image_0); //OJO

  fftw_complex *object=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n*n);
  fftw_complex *object_raw=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n*n); //OJO
  fftw_complex *hk=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n*n);

  //Initialize best guess of object to be reconstructed
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      //object[i*n+j][0]=1;
      object[i*n+j][0]=image_0[i*n+j][0];
      object[i*n+j][1]=0;
      hk[i*n+j][0]=object[i*n+j][0]; 
      hk[i*n+j][1]=object[i*n+j][1];
    }
  }

  hk=conv_fftw(n, hk, mtf, hk); 

  double error1, error2, emse, ecorr;

  fprintf(stderr, "RL Iteration: \n");
  //*************
  //FILE *error_file=fopen(error_name, "a");
  FILE *error_file=fopen(error_name, "w");
  for(int iter=0; iter<n_iter; iter++){ //RL main loop
    fprintf(stderr, "%d  ", iter+1);
    
    //Set image=image/hk (hk is the conv of psf and object estimate)
    for(int m=0; m<2; m++){
      for(int i=0; i<n; i++){
	for(int j=0; j<n; j++){
	  if(hk[i*n+j][m]==0) image[i*n+j][m]=0.0;	  
	  else image[i*n+j][m]=image_0[i*n+j][m]/hk[i*n+j][m]; //This is the important line in this loop
	}
      }
    }
    
    //Now convolve (image/conv(guess & psf)) and psf
    image=conv_fftw(n, image, mtf, image);

    //Now multiply ((image/conv(guess & psf)) with the guess
    for(int i=0; i<n; i++){
	for(int j=0; j<n; j++){
	  object[i*n+j][0]=image[i*n+j][0]*object[i*n+j][0];
	  object[i*n+j][1]=0; //Image must be real
	  object[i*n+j][0]=object[i*n+j][0]*mask[i*n+j][0]; 	  
	  hk[i*n+j][0]=object[i*n+j][0];
	  hk[i*n+j][1]=0; //Image must be real
	  if(iter==n_iter-1){ 
	    object_raw[i*n+j][0]=object[i*n+j][0];
	    object_raw[i*n+j][1]=0;
	  }
	}
    }
    object=normalize_fftw(n, object); //OJO

    hk=conv_fftw(n, hk, mtf, hk);

    hk=normalize_fftw(n, hk); //OJO

    //Print error
    double error_y=find_error(n, image_0, hk);
    double error_x=find_error(n, prist, object);
    double *corr_mse=new double[2];
    corr_mse=correlation(n, prist, object);    
    fprintf(stderr, "Error %e  %e  %e  %e\n", error_y, error_x, corr_mse[0], corr_mse[1]);
    fprintf(error_file, "%e  %e  %e  %e\n", error_y, error_x, corr_mse[0], corr_mse[1]);
    error1=error_y;
    error2=error_x;
    ecorr=corr_mse[0];
    emse=corr_mse[1];
  }//end main loop

  fclose(error_file);
  fprintf(stderr, "\n");

  fprintf(stdout, "%e  %e  %e  %e\n", error1, error2, ecorr, emse);

  //return object;
  return object_raw;

}


//=====================================================================
//=====================================================================


double **generate_pupil(int n, struct ppl *ppl){

  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){
      ppl->pupil[i][j]=0;
      ppl->pupil_centers[i][j]=0;
    }
  }

  for(int i=0; i<ppl->n_sub; i++){
    //if(ppl->diameter[i]<2) ppl->diameter[i]=2; //TODO: WHY CAN'T THIS BE 1?
    for(int l=ppl->pos_y[i]-int(ppl->diameter[i]/2); l<=ppl->pos_y[i]+int(ppl->diameter[i]/2); l++){
      for(int k=ppl->pos_x[i]-int(ppl->diameter[i]/2); k<=ppl->pos_x[i]+int(ppl->diameter[i]/2); k++){
	if((k-ppl->pos_x[i])*(k-ppl->pos_x[i])+(l-ppl->pos_y[i])*(l-ppl->pos_y[i])<=int(ppl->diameter[i]/2)*int(ppl->diameter[i]/2)) ppl->pupil[k][l]=1;
      }
    }
  }

  return ppl->pupil;

}

//=====================================================================
//=====================================================================

fftw_complex *normalize_fftw(int n, fftw_complex *a){
  double max_re=0;
  double max_im=0;
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(max_re<a[i*n+j][0]) max_re=a[i*n+j][0];
      if(max_im<a[i*n+j][1]) max_im=a[i*n+j][1];
    }
  }
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      a[i*n+j][0]=a[i*n+j][0]/max_re;
      a[i*n+j][1]=a[i*n+j][1]/max_re;
    }
  }
  
  return a;
}

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

//========================================
//========================================

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
  return centered_image;

  for(int i=0; i<n; i++){
    delete  centered_image[i];
  }
  delete  centered_image;

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

double *change_to_1d(int n, double **image, int q){
  double *data;
  int datalen=n*n;
  data = new double [2*datalen];
  int k=0;
  data[0]=0.0;
  for(int j=0; j<N; j++){
    for(int i=0; i<N; i++){
      if(q==0){
	data[k]=image[i][j];      
	if((k+1)<2*datalen) data[k+1]=0; //Imaginary part of image
      k=k+2;
      }
      if(q==1){
	data[k]=0;      	
	if((k+1)<2*datalen) data[k+1]=image[i][j]; //Imaginary part of image
	k=k+2;
      }    
    }    
  }
  return data;
  delete data;
}


//====================================================================
//====================================================================
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
  fprintf(stderr, "Physical pupil size is: %0.3f m\n", ppl->size);
  ppl->pos_x=new int[ppl->n_sub];
  ppl->pos_y=new int[ppl->n_sub];
  ppl->diameter=new int[ppl->n_sub];
  ppl->diameter_c=new int[ppl->n_sub];


  while(fscanf(pupil_file, "%d %d %d\n", &ppl->pos_x[i], &ppl->pos_y[i], &ppl->diameter[i])!=EOF ) {
    ppl->pos_x[i]=(int)floor(ppl->pos_x[i]*ppl->pscale+0.5);
    ppl->pos_y[i]=(int)floor(ppl->pos_y[i]*ppl->pscale+0.5);
    ppl->diameter[i]=(int)floor(ppl->diameter[i]*ppl->pscale+0.5);
    //ppl->diameter[i]=5;
    if(ppl->pos_x[i]>N || ppl->pos_x[i]>N){
      fprintf(stderr, "Pupil scaling error!\n");
      exit(1);
    }
    ppl->diameter_c[i]=2; 
    i++;    
 }

  fclose(pupil_file);

}

//======================================================
//======================================================

int readpilot(struct plt *plt) {
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

    //Pristine fits image name
    if (strstr(flag,"PRIST")!=NULL) {
      plt->prist_flag=1;
      get_word(fpilot, plt->pristine_name, &ieof);
    }
        
    //Pupil .par file name
    if (strstr(flag,"PPLPAR")!=NULL) {
      get_word(fpilot, plt->pupil_par_name, &ieof);
    }
    //Pupil pgm .par file name
    if (strstr(flag,"PGMPAR")!=NULL) {
      get_word(fpilot, plt->pgm_par_name, &ieof);
    }

    //Pupil fits file name
    if (strstr(flag,"PUPILN")!=NULL) {
      get_word(fpilot, plt->pupil_fits_name, &ieof);
    }

    //Pristine image name
    if (strstr(flag,"PRISTI")!=NULL) {
      get_word(fpilot, plt->pristine_name, &ieof);
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

    //Average number of photons per pixel
    if (strstr(flag,"PHOTON")!=NULL) {
      fscanf(fpilot,"%lf",&(plt->nphoton));
    }

    //Average number of photons per pixel
    if (strstr(flag,"TPHOTO")!=NULL) {
      fscanf(fpilot,"%lf",&(plt->tphoton));
    }

    //Average number of photons per pixel
    if (strstr(flag,"PSFOUT")!=NULL) {
      plt->psf_flag=1;
    }

    //Deconvolved image file name
    if (strstr(flag,"ANAFIZ")!=NULL) {
      get_word(fpilot, plt->de_densified_fits_name, &ieof);
    }

    //Data image file name
    if (strstr(flag,"DATAIN")!=NULL) {
      get_word(fpilot, plt->data_fits_name, &ieof);
    }

    //Number of iterations in RL deconcolution
    if (strstr(flag,"RLITER")!=NULL) {
      fscanf(fpilot,"%d",&(plt->nrl));
    }

    //Mask for deconvolution
    if (strstr(flag,"RLMASK")!=NULL) {
      plt->mask_flag=1;
      fscanf(fpilot,"%lf",&(plt->mask));
    }
    
    //Deconvolved image file name
    if (strstr(flag,"ERRORF")!=NULL) {
      get_word(fpilot, plt->error_name, &ieof);
    }

    //Deconvolved image file name
    if (strstr(flag,"RESULT")!=NULL) {
      get_word(fpilot, plt->rl_fits_name, &ieof);
    }
    
  }
  
  //Closes the pilot file
  fclose (fpilot);
  
  //TODO: POSSIBLE ERROR MESSAGES use flags
  fprintf(stderr, "Done reading the pilot file\n");
  return 0;
}

//==================================================
//==================================================

void initialize_pilot(struct plt *plt){

  plt->lambda=550E-9;
  plt->pupil_pgm_flag=0;
  plt->pupil_fits_flag=0;
  plt->psf_flag=0;
  plt->prist_flag=0;
  plt->dfcorr=0;
  plt->tphoton=0;
  plt->pilotfile=new char[STRING_LENGTH];  
  plt->data_fits_name=new char[STRING_LENGTH];  
  plt->pristine_name=new char[STRING_LENGTH];  
  plt->mod_fits_name=new char[STRING_LENGTH];  
  plt->mtf_fits_name=new char[STRING_LENGTH];  
  plt->pupil_par_name=new char[STRING_LENGTH];
  plt->pupil_pgm_name=new char[STRING_LENGTH];
  plt->pgm_par_name=new char[STRING_LENGTH];
  plt->pupil_fits_name=new char[STRING_LENGTH];
  plt->fizeau_fits_name=new char[STRING_LENGTH];
  plt->pupil_fits_read=new char[STRING_LENGTH];
  plt->de_densified_fits_name=new char[STRING_LENGTH];
  plt->fft_fizeau_name=new char[STRING_LENGTH];
  plt->fft_densified_name=new char[STRING_LENGTH];
  plt->densified_fits_name=new char[STRING_LENGTH];
  plt->error_name=new char[STRING_LENGTH];
  plt->rl_fits_name=new char[STRING_LENGTH];
  plt->mask_flag=0;

}

//===============================================
//===============================================

double **initialize_matrix(int n, double **image, double init){
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      image[i][j]=init;
    }
  }
  return image;
}

//===============================================
//===============================================

void checknan(int n, double **image){
  
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if(isnan(image[i][j])){
	fprintf(stderr, "NAN IMAGE\n");
	exit(1);
      }
    }
  }

}

//===============================================
//===============================================

double **read_file(int n, int offset, double range,  double **image, char image_file[]){
  //This function recieves the length of one side of the image matrix 'n' (the total size is n*n),
  //it also recieves  the image matrix itself 'image', and the file name.
  //Then it simply prints to a file whose name is contained in 'image_file[]'
 
  double step=double(range/N);

  int counter=0;
  FILE *img;
  img=fopen(image_file,"r");	//This image file can be either the reconstructed image or its fourier transform
  
  fprintf(stderr, "\nReading %s\n", image_file);
  fprintf(stderr, "step = %0.3f\nmaximum range = %0.3f \n", step, step*n/2);
  
  for(int j=0; j<n; j++){
    for(int i=0; i<n; i++){	
      double x=i*step-n/2*step;
      double y=j*step-n/2*step;
      if(i+offset<n && j+offset<n){
	fscanf(img, "%lf  %lf  %lf\n", &x, &y, &image[i+offset][j+offset]);
	if(image[i+offset][j+offset]!=0) counter++;	    
      }

    }
    fscanf(img, "\n");//This extra space is only for gnuplot purposes. 
  }	
  fclose(img);
  if(counter==0) fprintf(stderr, "Empty image!\n");
  if(counter!=0) fprintf(stderr, "Non empty image\n");

  return image;

}

//=====================================================================
//=====================================================================

void main_spot_find(int n, double threshold, struct ppl *ppl, double **image){
  //This function receives a pgm file, finds the sub-aperture positions and saves them in the ppl structure
  fprintf(stderr, "Finding spots in pupil\n");

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

double **remove_offset(int n, double **image, double offset){
  //Useful to remove constant background
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      image[i][j]=image[i][j]-offset;
      if(image[i][j]<0) image[i][j]=0; //Positivity constraint in data?
    }
  }

  return image;
}

//=====================================================================
//=====================================================================

double *correlation(int n, fftw_complex *image1, fftw_complex *image2){
  //This function finds the mean square error
  //The dimension of the correlation is size_corr*size*corr 

  double int1=0;
  double average1=0;
  for(int j=0; j<n; j++){ 
    for(int i=0; i<n; i++){ 
      average1=image1[i*n+j][0]+average1;
    }
  }
  int1=average1;
  average1=average1/(n*n);
  
  
  double int2=0;
  double average2=0;
  for(int j=0; j<n; j++){ 
    for(int i=0; i<n; i++){ 
      average2=image2[i*n+j][0]+average2;
    }
  }
  int2=average2;
  average2=average2/(n*n);
  
  //Find squared standard deviations
  double epsilon=0.1;
  double sum1=0, sum2=0;
  for(int l=0; l<n; l++){ 
    for(int k=0; k<n; k++){ 
      sum1=(image1[k*n+l][0]-average1)*(image1[k*n+l][0]-average1)+sum1;
    }
  }
  
  for(int l=0; l<n; l++){ 
    for(int k=0; k<n; k++){ 
      sum2=(image2[k*n+l][0]-average2)*(image2[k*n+l][0]-average2)+sum2;	
    }
  }
  sum1=sum1/(n*n); //Will divide correlation by sqrt of these values
  sum2=sum2/(n*n);
  
  
  double diff;
  time_t  t1, t2; //To calculate computation time
  //fprintf(stderr, "Calculating correlation of size %d\n", n);
  double sum=0;
  double max=0;
  int step=1; //step size in correlation computation
  
  double mse=0; //Mean square error
  double r=0; //Correlation
  
  sum=0;
  for(int j=0; j<n; j=j+step){//i and j are the summation indices
    for(int i=0; i<n; i=i+step){
	sum=(image1[i*n+j][0]-average1)*(image2[i*n+j][0]-average2)+sum;
	mse=(image1[i*n+j][0]/int1-image2[i*n+j][0]/int2)*(image1[i*n+j][0]/int1-image2[i*n+j][0]/int2)+mse;
    }
  }
  r=sum/(n*n)*1/sqrt(sum1*sum2);
  mse=mse/(n*n);
  
  //fprintf(stderr, "\ncorrelation is %f\n", r);
  //fprintf(stderr, "\nmse is %e\n\n", mse);
  
  double *rtn=new double[2];
  rtn[0]=r;
  rtn[1]=mse;
  
  return rtn;
  
}

//=====================================================================
//=====================================================================
