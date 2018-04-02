#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#define pi 3.1415926535897932384

FILE *dat;

int main(){
    int i,j,k;
    dat = fopen("toroide.txt","w");
    fprintf(dat," \"x\",\"y\",\"value\"\n");
    for(i=-10*13.5;i<=10*13.5;i++){
        for(j=-10*13.5;j<=10*13.5;j++){
            if((10-pow(pow(i/10.0,2)+pow(j/10.0,2),0.5),2>=0)&&(12.25-pow(10-pow(pow(i/10.0,2)+pow(j/10.0,2),0.5),2)>=0)){
                    fprintf(dat,"%f,%f,%f\n",i/10.0,j/10.0,pow(12.25-pow(10-pow(pow(i/10.0,2)+pow(j/10.0,2),0.5),2), 0.5));
            }
            else{
                    fprintf(dat,"%f,%f,%f\n",i/10.0,j/10.0,0.0);
            }
        }
    }
    fclose(dat);
    return(0);
}
