#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#define pi 3.1415926535897932384
float r, R;

FILE *dat;

int main(){
    int i,j,k;
    printf("Introduzca el radio mayor");
    scanf("%f",R);
    printf("Introduzca el radio menor");
    scanf("%f",r);
    dat = fopen("toroide.txt","w");
    fprintf(dat," \"x\",\"y\",\"value\"\n");
    for(i=-10*(r+R);i<=10*(r+R);i++){
        for(j=-10*(r+R);j<=10*(r+R);j++){
            if((R-pow(pow(i/10.0,2)+pow(j/10.0,2),0.5),2>=0)&&(r*r-pow(R-pow(pow(i/10.0,2)+pow(j/10.0,2),0.5),2)>=0)){
                    fprintf(dat,"%f,%f,%f\n",i/10.0,j/10.0,pow(r*r-pow(R-pow(pow(i/10.0,2)+pow(j/10.0,2),0.5),2), 0.5));
            }
            else{
                    fprintf(dat,"%f,%f,%f\n",i/10.0,j/10.0,0.0);
            }
        }
    }
    fclose(dat);
    return(0);
}
