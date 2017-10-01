#include "RPS.h"
#include <time.h>

void process_chunk(cell* buffer, int buf_size, int frames);

int main(int argc, char** argv){

  printf("WARNING: THIS CODE WILL GENERATE A TON OF DATA BEFORE MAKING THE VIDEO\n");
  printf("THE FRAMES WILL BE STORED IN THE data/ DIRECTORY. YOU MUST DELETE THIS AFTER RUNNING THIS PROGRAM\n");
  printf("Press any key to continue\n");
  getchar();

  srand(time(NULL));
  cell* img_buffer = malloc(BUF_SIZE);
  if(img_buffer == NULL){
    printf("Uh oh, seems like your OS doesn't want to allocate the video buffer\n");
    printf("Try to change the BUF_SIZE define in bitmap.h (divide it by 2?)");
  }

  int frame_size = sizeof(cell)*IMG_X*IMG_Y;
  int frames_per_buffer = BUF_SIZE/frame_size;

  printf("frame size: %d kb\n", frame_size/1024);
  printf("frames per buffer: %d\n", frames_per_buffer);


  // Yep, that's a bona fide triple pointer
  // This makes sense when you remember we have many images, a la images[ImageNo][x][y]
  cell*** images = calloc(1024, sizeof(cell**));
  for(int ii = 0; ii < 1024; ii++){
    images[ii] = malloc(img_buffer, ii);
  }


  // Seed some CAs
  for(int ii = 0; ii < 100; ii++){
    int rx = rand() % (IMG_X - 1);
    int ry = rand() % (IMG_Y - 1);
    int rt = rand() % 4;

    images[0][rx][ry].color = rt;
    images[0][rx][ry].strength = 1;
    printf("NEW CELL (%d, %d) with color %d and strength %d\n", rx, ry, rt, 1);
  }


  int frames_total = 5000;
  // In the for loop further down you can see why we subtract 1.
  // we would run out of bounds by 1 if we used all frames in the buffer
  int temp = frames_total % (frames_per_buffer - 1);
  int chunks = frames_total/(frames_per_buffer - 1);
  if(temp != 0){ chunks++; }


  for(int ii = 0; ii < chunks; ii++){

    printf("processing chunk %d out of %d\n", ii, chunks);
    int frames_to_render;
    if(((ii + 1)*(frames_per_buffer - 1)) < frames_total){
      frames_to_render = (frames_per_buffer - 1);
    }
    else{
      frames_to_render = frames_total - (ii*(frames_per_buffer - 1));
    }


    printf("rendering frames\n");
    for(int jj = 0; jj < frames_to_render; jj++){
      iterate_image(images[jj], images[jj+1]);
    }

    printf("saving frames to disk\n");
    for(int jj = 0; jj < frames_to_render; jj++){
      make_bmp(images[jj], (jj + ii*(frames_per_buffer - 1)));
    }

    for(int jj = 0; jj < IMG_Y*IMG_X; jj++){
      img_buffer[jj] = img_buffer[(IMG_X*IMG_Y*(frames_per_buffer - 1)) + jj];
    }

  }

  free(img_buffer);
}
