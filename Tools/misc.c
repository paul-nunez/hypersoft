double **allocate_matrix(int i_max, int j_max, double **image); 
double **change_to_2d(int n, double *image, double **data, int q);
int get_word(FILE *funit, char name[STRING_LENGTH], int *ieof);

//====================================================================
//====================================================================

double **allocate_matrix(int i_max, int j_max, double **f ){
  //Allocates memory for a matrix
  //double **f=new double*[i_max];
  f=new double*[i_max];
  for(int j=0; j<j_max; j++){
    f[j]=new double[j_max];
  }
  return f;
}

//=====================================================================
//=====================================================================

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

//=====================================================================
//=====================================================================
int get_word(FILE *funit, char name[STRING_LENGTH], int *ieof)
//TODO: MOVE TO TOOLS
//Read a word defined a a string with no blanks. The first encountered 
//blank will indicate the end of the string. 
{
  char last_char;          /*Last character read on the line*/ 
  int fnamlen=0;           /*File name length*/
  
  while(!*ieof){
    /*we test for end of file*/
    *ieof=feof(funit);
    if(*ieof) return 0;  
    /*reads one character*/
    fscanf(funit,"%c",&last_char);
    
    /*skip any blank before the name*/
    if (last_char==' ' && fnamlen==0) continue;
    /*any blank encountered when reading the name indicates 
      the end of the name*/
    if (last_char==' ' || last_char=='\n') break;
    
    name[fnamlen]=last_char;
    fnamlen=fnamlen+1;
    if (fnamlen>STRING_LENGTH) 
      return fprintf(stderr, "Word too long");
  }
  name[fnamlen]='\0';
  return 0;   
}

//================================================================
//================================================================



int look_for_flag(int *ieof, char flag[FLAG_LENGTH], FILE *unit)
{
  //This function looks for a flag which is a * followed by a string of length 
  //FLAGLENGTH. There maybe a number of blanks between the * and the string. 
  //They will be ignored. The returned flag is not string not including the *
  char last_char;
  int line_begins=1;
  int i;
  
  //As long as end of file has not been reached
  while(!*ieof){
    //we test for end of file
    *ieof=feof(unit);
    if(*ieof) return 0;  
    //reads one character
    fscanf(unit,"%c",&last_char);
    //skip any blank
    if (last_char==' ') continue;
    //At this point we found a character that is not a '' 
    if (line_begins){
      line_begins=0;
      if (last_char=='*') {
	last_char=' ';
	//skips all space following the '*' 
	while(last_char==' ') {
	  *ieof=feof(unit);
	  if(*ieof) return 0;  
	  fscanf(unit,"%c",&last_char);
	}
	//We found the first character of a flag
	flag[0]=last_char;
	//and we read the others
	for (i=1;i<FLAG_LENGTH;i++) {
	  *ieof=feof(unit);
	  if(*ieof) return 0;  
	  fscanf(unit,"%c",&flag[i]);
	}
	/*we found a flag and now return*/
	return 0;
      }
    } 
    /*if we found a carriage return we are going to start a new line*/
    if (last_char=='\n') line_begins=1;
  }
  return 0;
}

//================================================================
//================================================================
