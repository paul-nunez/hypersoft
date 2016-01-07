#include <fitsio.h>
#ifndef MACFITS
//#include "../cfitsio/include/fitsio.h"
#endif
#ifdef MACFITS
#include "/cfitsio/include/fitsio.h"
#endif

void printerror( int status);

void fitsout(int mmax, float step,  double **image, char image_file[]);

//==================================================
//==================================================
void fitsout(int mmax, float step,  double **image1, char image_file[]){
  //Function that recieves the square root of the  dimension of the imag, the image, and the output name
  //It generates fits output

  fprintf(stderr, "\nCreating '%s' output\n", image_file);

  float xstep=step;
  float ystep=step;
  int maxx=mmax;
  int maxy=mmax;
  long naxes[2] = { maxx, maxy }; // image is mmaxx pixels wide by mmaxy 
  remove(image_file);   // Delete old file if it already exists 
  
  fitsfile *fptr;  // pointer to the FITS file, defined in fitsio.h 
  int status=0;   

  double **image;
  image=new double*[mmax];
  //double *image[mmax];
  
  // allocate memory for the whole image 
  image[0] = (double*)malloc( naxes[0] * naxes[1]
                                * sizeof(double) );


  // initialize pointers to the start of each row of the image 
  for( int i=1; i<naxes[1]; i++ )
    image[i] = image[i-1] + naxes[0];
  

  // initialize the values in the image with values of image1
  //For reasons that I don't understand I need to do this because fits_write_img does not understand **
  for (int j = 0; j < naxes[1]; j++)
    {   for (int i = 0; i < naxes[0]; i++)
        {
	  image[j][i] = image1[i][j];
        }
    }
  

  /*
  //For some reason this declaration does not work. fits_write_img does not seem to work with **
  //image=new double *[mmax];
  for(int i=0; i<mmax; i++){
    image[i]=new double[mmax];
    }

  for(int i=0;i<mmax;i++){
    for(int j=0;j<mmax;j++){
      //image[i][j]=3.5*cos(sqrt((i-0.5*mmax)*(i-0.5*mmax)+(j-0.5*mmax)*(j-0.5*mmax))*0.1);
      image[i][j]=image1[i][j];
    }
    }
  */

  if (fits_create_file(&fptr, image_file, &status)) // create new FITS file 
    printerror( status );           // call printerror if error occurs 
  
  
  int bitpix   =  DOUBLE_IMG; // 16-bit unsigned short pixel values //
  long naxis    =   2;  // 2-dimensional image
  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    printerror( status );          
  
  // write the array of unsigned integers to the FITS file
  long fpixel=1; //first pixel to write
  long nelements=naxes[0] * naxes[1]; //number of pixels to write
  //if ( fits_write_img(fptr, TDOUBLE, fpixel, nelements, image, &status) )
  if ( fits_write_img(fptr, TDOUBLE, fpixel, nelements, image[0], &status) )
    printerror( status );
  
  free( image[0] );  /* free previously allocated memory */

  //Header information
  /*char radesys[3];
  strcpy(radesys,"LIN");
  char paulsfix[15];
  strcpy(paulsfix,"RADESYS");
  char paulsfix2[25];
  strcpy(paulsfix2, "Linear");
  if ( fits_update_key(fptr, TSTRING, paulsfix, radesys,
    paulsfix2, &status) )
    printerror( status );
  printf("hello\n");
  char ctype1[1];
  strcpy(ctype1,"");
  strcpy(paulsfix, "CTYPE1");
  strcpy(paulsfix2, "?");
  
  if ( fits_update_key(fptr, TSTRING, paulsfix, ctype1,
		       paulsfix2, &status) )
    printerror( status );           

  char ctype2[1];
  strcpy(ctype2,"");
  strcpy(paulsfix, "CTYPE2");
  if ( fits_update_key(fptr, TSTRING, paulsfix, ctype2,
		       paulsfix2, &status) )
    printerror( status );           
  
  
  float xxref;  //CRVAL1
  if(maxx%2==0) xxref=-0.5*xstep;
  else xxref=0.0;
  strcpy(paulsfix, "CRVAL1");
  strcpy(paulsfix2, "Total Exposure Time");
  if ( fits_update_key(fptr, TFLOAT, paulsfix, &xxref,
		       paulsfix2, &status) )
    printerror( status );           
  float yyref=0.0;  //CRVAL2
  if(maxy%2==0) yyref=-0.5*ystep;
  else yyref=0.0;
  strcpy(paulsfix, "CRVAL2");
  if ( fits_update_key(fptr, TFLOAT, paulsfix, &yyref,
		       paulsfix2, &status) )
    printerror( status );           
  
  long xiref; //CRPIX1
  if(maxx%2==0) xiref=maxx/2;
  else xiref=(maxx-1)/2+1;
  strcpy(paulsfix, "CRPIX1");
  strcpy(paulsfix2, "Central pixel in x");
  if ( fits_update_key(fptr, TLONG, paulsfix, &xiref,
		       paulsfix2, &status) )
    printerror( status );
  
  long yiref=1; //CRPIX2
  if(maxy%2==0) yiref=maxy/2;
  else yiref=(maxy-1)/2+1;
  strcpy(paulsfix, "CRPIX2");
  strcpy(paulsfix2, "Central pixel in y");
  if ( fits_update_key(fptr, TLONG, paulsfix, &yiref,
		       paulsfix2, &status) )
    printerror( status );           
  
  float xscale=xstep;
  strcpy(paulsfix, "CDELT1");
  strcpy(paulsfix2, "X scale");
  if ( fits_update_key(fptr, TFLOAT, paulsfix, &xscale,
		       paulsfix2, &status) )
    printerror( status );           
  float yscale=ystep;
  strcpy(paulsfix, "CDELT2");
  strcpy(paulsfix2, "Y scale");
  if ( fits_update_key(fptr, TFLOAT, paulsfix, &yscale,
		       paulsfix2, &status) )
    printerror( status );           
  */
  if ( fits_close_file(fptr, &status) ) // close the file 
    printerror( status );           
  
  return;
}

//==================================================
//==================================================

double *readfits(char *fits_name, int nmax)

    /************************************************************************/
    /* Read a FITS image and determine the minimum and maximum pixel values */
    /************************************************************************/
{
  double *output=new double[2*nmax*nmax];// PAUL
  fprintf(stderr, "\nReading fits data\n");

    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status,  nfound, anynull;
    long naxes[2], fpixel, nbuffer, npixels, ii;

    #define buffsize 1000
    //#define buffsize 1000000
    float datamin, datamax, nullval, buffer[buffsize];
    //char filename[]  = "atestfil.fit";     /* name of existing FITS file   */
    char *filename=new char[50];
    strcpy(filename, fits_name);

    status = 0;

    if ( fits_open_file(&fptr, filename, READONLY, &status) )
         printerror( status );

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
         printerror( status );

    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1;
    nullval  = 0;                /* don't check for null values in the image */
    datamin  = 1.0E30;
    datamax  = -1.0E30;

    int i=0;
    while (npixels > 0)
    {
      nbuffer = npixels;
      if (npixels > buffsize)
        nbuffer = buffsize;     /* read as many pixels as will fit in buffer */

      /* Note that even though the FITS images contains unsigned integer */
      /* pixel values (or more accurately, signed integer pixels with    */
      /* a bias of 32768),  this routine is reading the values into a    */
      /* float array.   Cfitsio automatically performs the datatype      */
      /* conversion in cases like this.                                  */

      if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
                  buffer, &anynull, &status) )
           printerror( status );

      for (ii = 0; ii < nbuffer; ii++)  {
	//if(isnan(double(buffer[ii]))) buffer[ii]=0;
	output[i]=double(buffer[ii]); //PAUL copy image
	output[i+1]=0; //Imaginary part (assumes image is real)
	if(isnan(double(output[i]))){fprintf(stderr, "Nan error in readfits\n"); exit(1);}
	i=i+2;

        if ( buffer[ii] < datamin )
            datamin = buffer[ii];

        if ( buffer[ii] > datamax )
            datamax = buffer[ii];
      }
      npixels -= nbuffer;    /* increment remaining number of pixels */
      fpixel  += nbuffer;    /* next pixel to be read in image */

    }

    //printf("\nMin and max image pixels =  %.0f, %.0f\nN pix %d\n", datamin, datamax, i);
    //fprintf(stderr, "N pixels = %d\n", i);

    if ( fits_close_file(fptr, &status) )
         printerror( status );


    return output;
}


//==================================================
//==================================================
void printerror( int status)
{
  // Print out cfitsio error messages and exit program 
  if (status){
    fits_report_error(stderr, status); // print error report  
    exit( status ); // terminate the program, returning error status 
  }
  return;
}
//==================================================
//==================================================
