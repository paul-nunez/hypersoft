#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define L 50
#define N 1024

main(int n, char **arg){
  
  if(arg[1]==NULL || strcmp(arg[1], "-h")==0 || n<6){
        fprintf(stderr, "\nPlease use kewords followed by the argument\n\n");    
	fprintf(stderr, "-o Output parameter file path (.par file)\n");
	fprintf(stderr, "-n Number of sub-pupils\n");
	fprintf(stderr, "-d Diameter of sub-pupils\n");
	fprintf(stderr, "-s Physical size in meters of pupil\n");
	fprintf(stderr, "-l Length (m) of mirror array \n");
	fprintf(stderr, "-w Width (m) of mirror array \n");
	fprintf(stderr, "-h Prints this help message\n\n");
	fprintf(stderr, "Example:\n ./generate_periodic_pupil -o pupil.par -n 100 -s 1000 -l 1000 -w 1000 -d 1\n");
	exit(1);
  }  

  char *output=new char[L];
  int number=-1;
  float len=-1;
  float wid=-1;
  float diam=-1;
  float size=-1;
    
  for(int i=0; i<n; i++){
      if (strcmp(arg[i], "-o")==0) strcpy(output, arg[i+1]);
      if (strcmp(arg[i], "-n")==0) number=atoi(arg[i+1]);
      if (strcmp(arg[i], "-l")==0) len=atof(arg[i+1]);      
      if (strcmp(arg[i],"-w")==0) wid=atof(arg[i+1]);      
      if (strcmp(arg[i],"-d")==0) diam=atof(arg[i+1]);      
      if (strcmp(arg[i],"-s")==0) size=atof(arg[i+1]);      
  }

  fprintf(stderr, "Parameter file: %s\n", output);
  fprintf(stderr, "Number of sub-pupils: %d\n", number);
  fprintf(stderr, "Diameter: %f\n", diam);
  fprintf(stderr, "Length: %f\n", len);
  fprintf(stderr, "Width: %f\n", wid);  

  if(output==NULL || number==-1 || len==-1 || wid==-1 || diam==-1 || size==-1){
    fprintf(stderr, "Execution error\n");
    exit(1);
  }

  FILE *pupil_par=fopen(output, "w");

  int length=int(floor(len*N/size+0.5));
  int width=int(floor(wid*N/size+0.5));
  int diameter=int(floor(diam*N/size+0.5));

  int nx=int(floor(length/sqrt(number)+0.5));
  int ny=int(floor(width/sqrt(number)+0.5));

  fprintf(stderr, "length (pixels): %d\nwidth (pixels): %d\nstep x: %d\nstep y: %d\n", length, width, nx, ny);

  float step=size/N;
  float dtheta=2*3.141592/number;


  fprintf(pupil_par, "%d\n%f\n%f\n", number, size, step);

  for(int k=0; k<number; k++){
    int i=N/2+int(floor(length/2*cos(k*dtheta)+0.5));
    int j=N/2+int(floor(length/2*sin(k*dtheta)+0.5));
    fprintf(pupil_par, "%d  %d  %d\n", i, j, diameter);
  }
  fclose(pupil_par);

  return 0;

}
