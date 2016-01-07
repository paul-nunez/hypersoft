 
struct NODE {
  int number;
  int j;
  int head;
  int tail;
  struct NODE *next;
};


struct llist {
  int number;
  int j; //pixel reference
  int head;
  int tail;
  struct NODE *next;
};


struct spot_params{
    float brightness;
    float xcm;
    float ycm;
    float radius;
    int counter;
};

struct image{
  float *pixel;
  float *x;
  float *y;
};



//*********************


float *find_spot(float *pixel, int n, float threshold, float replace, struct spot_params *params, struct image *img);

void append_node(struct NODE *llist, int num, int i);

void delete_node(struct NODE *llist, int num);

int search_value(struct NODE *llist, int num);

void add_to_params(int j, float pixel, struct spot_params *params, struct image *img);



//**********************

float *find_spot(float *pixel, int n, float threshold, float replace, struct spot_params *params, struct image *img){

  int flag=0;
  int num = 0;
  int input = 1;
  int retval = 0;
  struct NODE *llist;
  
  llist = (struct NODE *)malloc(sizeof(struct NODE));
  llist->number = 0;
  llist->j=0;
  llist->next = NULL;
  
  //Find the first pixel in spot
  params->brightness=0;
  params->xcm=0;
  params->ycm=0;
  params->counter=0;
  int i=(int)sqrt(n)+2;  
  //int i=(int)sqrt(n);  
  //printf("N %d\n", i);
  while(flag==0){
    if(pixel[i]>=threshold){
      //fprintf(stderr, "Pixel value = %f\n", pixel[i]);
      add_to_params(i, pixel[i], params, img); //Aids in brightness and cm calculation
      pixel[i]=replace;
      num=1;
      llist->head=num;
      append_node(llist, num, i);   
      flag=1;
      //fprintf(stderr, "%d %lf %lf", i, pixel[i], pixel[i+1]);
      //frintf(stderr, "\nFound spot at j=%d (replacement = %f)\n", i-1, replace);
    }  
    if(i==n-1){
      fprintf(stderr, "\nNo spot found\n");
      params->brightness=-1.0;
      flag=1;
    }

    i++;    

  }


  //Now start adding nodes as long as neighbors meet threshold
  flag=0;
  int j=i;

  if(i==n){
    flag=1; //No spot found
  }

  while(flag==0){ //Flag will be 1 when list is done
    j=search_value(llist, llist->head); //Search for head of list
    //j=num; //CHANGE
    
    //Now search all of it's neighbors to see if we need to add elements to the list
    if(pixel[j]>=threshold){
      add_to_params(j, pixel[j], params, img); //Aids in brightness and cm calculation
      num++; 
      llist->head=num;
      //num++; 
      pixel[j]=replace; 
      flag=1;
    }    
    if(pixel[j+1]>=threshold){
      add_to_params(j+1, pixel[j+1], params, img);
      num++; 
      llist->tail=num;
      append_node(llist, num, j+1); 
      pixel[j+1]=replace; 
      flag=1;
    }
    if(pixel[j-1]>=threshold){
      add_to_params(j-1, pixel[j-1], params, img);
      num++; 
      llist->tail=num;
      append_node(llist, num, j-1); 
      pixel[j-1]=replace; 
      flag=1;
    }
    if(pixel[j+(int)sqrt(n)-1]>=threshold){ 
      add_to_params(j+(int)sqrt(n)-1, pixel[j+(int)sqrt(n)-1], params, img);
      num++; 
      llist->tail=num;
      append_node(llist, num, j+(int)sqrt(n)-1); 
      pixel[j+(int)sqrt(n)-1]=replace; 
      flag=1;
    }
    if(pixel[j+(int)sqrt(n)]>=threshold){ 
      add_to_params(j+(int)sqrt(n), pixel[j+(int)sqrt(n)], params, img);
      num++; 
      llist->tail=num;
      append_node(llist, num, j+(int)sqrt(n)); 
      pixel[j+(int)sqrt(n)]=replace; 
      flag=1;
    }
    if(pixel[j+(int)sqrt(n)+1]>=threshold){ 
      add_to_params(j+(int)sqrt(n)+1, pixel[j+(int)sqrt(n)+1], params, img);
      num++; 
      llist->tail=num;
      append_node(llist, num, j+(int)sqrt(n)+1); 
      pixel[j+(int)sqrt(n)+1]=replace; 
      flag=1;
    }
    if(pixel[j-(int)sqrt(n)-1]>=threshold && j-(int)sqrt(n)-1>=0){ //CHANGE oct 31 2012
      add_to_params(j-(int)sqrt(n)-1, pixel[j-(int)sqrt(n)-1], params, img);
      num++; 
      llist->tail=num;
      append_node(llist, num, j-(int)sqrt(n)-1); 
      pixel[j-(int)sqrt(n)-1]=replace; 
      flag=1;
    }
    if(pixel[j-(int)sqrt(n)]>=threshold){ 
      add_to_params(j-(int)sqrt(n), pixel[j-(int)sqrt(n)], params, img);
      num++; 
      llist->tail=num;
      append_node(llist, num, j-(int)sqrt(n)); 
      pixel[j-(int)sqrt(n)]=replace; 
      flag=1;
    }
    if(pixel[j-(int)sqrt(n)+1]>=threshold){ 
      add_to_params(j-(int)sqrt(n)+1, pixel[j-(int)sqrt(n)+1], params, img);
      num++; 
      llist->tail=num;
      append_node(llist, num, j-(int)sqrt(n)+1); 
      pixel[j-(int)sqrt(n)+1]=replace; 
      flag=1;
    }
    if(num==1){ //If there is only one pixel that meets threshold, the loop should stop
      llist->head=llist->tail;
      flag=1;
    }
    if(flag==1){
      delete_node(llist, llist->head); //Delete node whose neighbors have already been checked
    }

    flag=0;

    if(llist->head==llist->tail) flag=1;        //Done with list
    if(llist->head!=llist->tail) llist->head++; //go to next item on list
    
    //fprintf(stderr, "(%d, %d), ", llist->head, llist->tail);

  }

  /*
  FILE *mod;
  mod=fopen("modified_image.txt", "w");
  for(int k=0; k<n; k++){    
    int l=k%(int)sqrt(n);
    int m=(int)k/sqrt(n);
    fprintf(mod, "%d  %d  %f\n", l, m, pixel[k]);
    if(l==sqrt(n)-1) fprintf(mod, "\n");
  }
  fclose(mod);
  */
  float step=img->x[1]-img->x[0];
  params->xcm=params->xcm/params->brightness;
  params->ycm=params->ycm/params->brightness;
  params->radius=sqrt(params->counter*step*step/PI);
  //fprintf(stderr, "brightness = %f\nxcm = %f\nycm = %f\nradius = %f\n", params->brightness, params->xcm, params->ycm, params->radius);

  return pixel;

}


//************************

void append_node(struct NODE *llist, int num, int i) {
 while(llist->next != NULL)
  llist = llist->next;

 llist->next = (struct NODE *)malloc(sizeof(struct NODE));
 llist->next->number = num;
 llist->next->j=i;
 llist->next->next = NULL;
}

//************************


void delete_node(struct NODE *llist, int num) {
 struct NODE *temp;
 temp = (struct NODE *)malloc(sizeof(struct NODE));

 if(llist->number == num) {
  /* remove the node */
  temp = llist->next;
  free(llist);
  llist = temp;
 } else {
  while(llist->next->number != num)
   llist = llist->next;

  temp = llist->next->next;
  free(llist->next);
  llist->next = temp;
 }   
}


//************************

int search_value(struct NODE *llist, int num) {
  //printf(".", num);
  
  int retval = -1;
  int i = 1;
  
  while(llist->next != NULL) {
    if(llist->next->number == num)
      return llist->next->j;
    else
      i++;
    
    llist = llist->next;
  }
  
  return retval;
}


//*************************

void add_to_params(int j, float pixel, struct spot_params *params, struct image *img){
  //Aids in brightness and cm calculation
  params->brightness=params->brightness+pixel; //Add to brightness
  params->xcm=pixel*img->x[j]+params->xcm;
  params->ycm=pixel*img->y[j]+params->ycm;
  params->counter++;
}

//*************************
