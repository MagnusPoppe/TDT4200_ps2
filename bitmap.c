#include "bitmap.h"

/* save 24-bits bmp file, buffer must be in bmp format: upside-down */
void savebmp(char *name,uchar *buffer,int x,int y) {
	FILE *f=fopen(name,"wb");
	if(!f) {
		printf("Error writing image to disk.\n");
		return;
	}
	unsigned int size=x*y*3+54;
	uchar header[54]={'B','M',size&255,(size>>8)&255,(size>>16)&255,size>>24,0,
                    0,0,0,54,0,0,0,40,0,0,0,x&255,x>>8,0,0,y&255,y>>8,0,0,1,0,24,0,0,0,0,0,0,
                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	int ierr = fwrite(header,1,54,f);
    if (ierr != 54) printf("ERROR WRITING HEADER.");
	ierr = fwrite(buffer,1,IMG_X*IMG_Y*3,f);
    if (ierr != IMG_X*IMG_Y*3) printf("ERROR WRITING IMAGE.");
	fclose(f);
}


/* given iteration number, set a colour */
void colorize(uchar* p, cell my_cell) {
	if(my_cell.color == WHITE){
    p[0] = 255;
    p[1] = 255;
    p[2] = 255;
  }
	else if(my_cell.color == ROCK) {
    p[0] = 0;
    p[1] = 0;
    p[2] = 255;
  }
	else if(my_cell.color == SCISSOR) {
    p[0] = 0;
    p[1] = 255;
    p[2] = 0;
  }
	else if(my_cell.color == PAPER) {
    p[0] = 255;
    p[1] = 0;
    p[2] = 0;
  }
}

void make_bmp(cell** image, int index){

  // create nice image from iteration counts. take care to create it upside
  // down (bmp format)

  // used to be calloc with 1 as arg 2, thus no sizeof(??)
  unsigned char* buffer = malloc(IMG_X*IMG_Y*3);
  for(int i=0;i<IMG_X;i++) {
    for(int j=0;j<IMG_Y;j++) {
      int p=((IMG_Y-j-1)*IMG_X+i)*3;
      colorize(buffer+p, image[i][j]);
    }
  }

  char filename [50];
  sprintf(filename, "data/CA-%d.bmp", index);
  /* write image to disk */
  savebmp(filename,buffer,IMG_X,IMG_Y);
  free(buffer);
}
