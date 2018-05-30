#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <tgmath.h>

#define MAXPART 1000
#define pi 3.1415926535897932384
#define epce 8.8541878176e-12
#define epsi 1//78.5
#define qe 1.60219e-19
#define kb 1.3806488e-23
#define na 6.02217e23
const int maxl = 20;
const int maxw = 20;

float alea(void);
void leer_datos_iniciales(void);
void imprimir_datos_iniciales(void);
void imprimir_celda_plasma(int a);
void condiciones_iniciales(void);
void arreglo_inicial(void);
void crear_matriz_plasma(void);
void calc_carga(int a);

float distanciacelda(int a, int b);
int signo(float a);
void mover_particulas(void);
void metropolis_plasma(int a);

float de_plasma(void);
float calcular_de_mov(void);
void gdr_plasma(void);
void actu_salida(void);
float energia(void);

//////////////////////////////Constantes
int pasos, terma, actu;
float esc, densidad, volumen, nelectronesr, tempe, tempee, tempei;
long long int nh20, nhp, nh2p, nelectrones;
float R;
float m1, m2, m3;
int v1, v2;
///////////////////////////////////////////////////////////////////////Variables globales
int p=0, tipo;
int rechazo;
int nceldas, n1, n2, n1i, n2i;                                //numero de celda Nueva/Vieja Fnicial/Final

///////////////////////////////////////////////////////////////////////Contadores
int c_mov=1,c_mova=1;
int rechazo_neg=0, rechazo_met_mov=0;

struct smatriz_plasma{
	float rho,phi,x,y;
	int carga,h20,hp,h2p,electrones;
    float k1x, k1y, k1z, k2x, k2y, k2z, volumen;
}matriz_plasma[61*50+1];
////////////////////////////////////////////////////////////////////////Variables para salida
FILE *dat,*dat2;
char salidac[10];

main(){
    int i,j,k;
    srand((unsigned)time(NULL));
    system("mkdir datos");

	condiciones_iniciales();
	crear_matriz_plasma();
	arreglo_inicial();

    c_mov=c_mova=1;
    rechazo_neg=rechazo_met_mov=0;

    /*for(i=1;i<=npart;i++){
        printf("\ni:%i,x:%f,y:%f,z:%f,q:%f",i,part[i].x,part[i].y,part[i].z,part[i].carga);
    }
    getchar();*/
    for(p=1;p<=pasos;p++){
        rechazo = 0;
        c_mov++;
        mover_particulas();

        if(rechazo == 0){
            metropolis_plasma(tipo);
        }
        else{
            rechazo_neg++;
            matriz_plasma[n1] = matriz_plasma[n1i];
            matriz_plasma[n2] = matriz_plasma[n2i];
        }
        if(p%actu==0){
            actu_salida();
        }
        if(p%1000==0){
            if(p%10000==0){
                printf("\r   Paso: %i Acept. Mov: %1.5f Rechazo negativas: %1.5f Rechazo met: %1.5f Energia: %e",p,c_mova/(c_mov*1.0),rechazo_neg/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0),energia());
            }
            else{
                printf("\r   Paso: %i Acept. Mov: %1.5f Rechazo negativas: %1.5f Rechazo met: %1.5f",p,c_mova/(c_mov*1.0),rechazo_neg/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0));
            }
        }
    }
}
///////////////////////////////////////////////////////////////////////////////////FIN DE MAIN
////////////////////////////////////////////////////////////////////////////////////
float alea(void){
return((float)rand()/RAND_MAX);
}
////////////////////////////////////////////////////////////////////////////////////
void leer_datos_iniciales(){
    dat=fopen("entrada.txt","r");
    fscanf(dat,"Numero de pasos: %i\n", &pasos);
    fscanf(dat,"Actualizacion: %i\n", &actu);
    fscanf(dat,"Termalizacion: %i\n", &terma);
    fscanf(dat,"Valencias: %i, %i\n",&v1,&v2);
    fscanf(dat,"Masas (kg): %f, %f, %f\n",&m1,&m2,&m3);
    fscanf(dat,"R (cm): %f\n",&R);
    fscanf(dat,"Temperatura (K): %f\n",&tempe);
    fscanf(dat,"Temperatura electrones (eV): %f\n",&tempee);
    fscanf(dat,"Temperatura iones (eV): %f\n",&tempei);
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_datos_iniciales(){
    printf("Numero de pasos: %i\n", pasos);
    printf("Actualizacion: %i\n", actu);
    printf("Termalizacion: %i\n", terma);
    printf("Valencias: %i, %i\n",v1,v2);
    printf("Masas (kg): %e, %e, %e\n",m1,m2,m3);
    printf("R (cm): %f\n",R);
    printf("Temperatura (K): %f\n",tempe);
    printf("Temperatura electrones (eV): %f\n",tempee);
    printf("Temperatura iones (eV): %f\n",tempei);
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_iniciales(){
    leer_datos_iniciales();
    esc = 1e-2;//diam*1e-9;

    densidad = 2e19;
    volumen = pi*R*1e-2*R*1e-2*0.125*1e-6;
    nelectronesr = (densidad*volumen);
    nelectrones = nelectronesr;
    nh20 = 9*nelectrones;
    nhp = 0.8*nelectrones;
    nh2p = 0.2*nelectrones+1;
    printf("\ndensidad: %e volumen: %e nnr: %e nelectrones: %I64d",densidad,volumen,nelectronesr,nelectrones);
    printf("\nnh20: %I64d nhp: %I64d nh2p: %I64d cargatotal: %I64d\n",nh20,nhp,nh2p,nh2p+nhp-nelectrones);
    //getchar();

    n1i = 1408;
    n2i = 1409;

    tempei = tempei*qe/kb;
    tempee = tempee*qe/kb;

    printf("\ntempei: %f tempee: %f", tempei, tempee);
    //getchar();

    imprimir_datos_iniciales();
}

////////////////////////////////////////////////////////////////////////////////////
void crear_matriz_plasma(){
    int i,j,contadorm=0,contadorvec=0;
    float RR;
    densidad = 2e19;
    volumen = pi*R*1e-2*R*1e-2*0.125*1e-6;
    RR = 60000000;                                                                   //Radio en cm, todo lo demas en nm
    contadorm++;
    matriz_plasma[1].rho = 1/sqrt(2);
    matriz_plasma[1].phi = pi/8;
    matriz_plasma[1].x = matriz_plasma[1].rho*cos(matriz_plasma[1].phi);
    matriz_plasma[1].y = matriz_plasma[1].rho*sin(matriz_plasma[1].phi);
    matriz_plasma[1].volumen = pi/8;
    for(i=2;i<=(int)(RR/(1000000));i++){
        for(j=1;j<=(int)((pi*i)*0.25);j++){
            contadorm++;
            matriz_plasma[contadorm].rho = sqrt(pow(i,2)+pow(i-1,2))/sqrt(2);
            matriz_plasma[contadorm].phi = pi/(8*((int)((pi*i)*0.25)))+(j-1)*(pi/(4*((int)((pi*i)*0.25))));
            matriz_plasma[contadorm].x = matriz_plasma[contadorm].rho*cos(matriz_plasma[contadorm].phi);
            matriz_plasma[contadorm].y = matriz_plasma[contadorm].rho*sin(matriz_plasma[contadorm].phi);
            matriz_plasma[contadorm].volumen =  (pi*(i*i-(i-1)*(i-1))*0.125*1e-2*1e-2*1e-6)/(int)((pi*i)*0.25);
        }
    }
    nceldas=contadorm;
}
////////////////////////////////////////////////////////////////////////////////////
void calc_carga(int a){
    matriz_plasma[a].carga = matriz_plasma[a].h2p + matriz_plasma[a].hp - matriz_plasma[a].electrones;
}
////////////////////////////////////////////////////////////////////////////////////
void arreglo_inicial(){
    int i,ni,dummy,res=1000;

    dat = fopen("datos/posiciones.dat","r");
    for(i=1;i<=1407;i++){
        fscanf(dat,"%i\t%i\t%i\t%i\t%i\t%i\n",&dummy,&matriz_plasma[i].electrones,&matriz_plasma[i].h20,&matriz_plasma[i].hp,&matriz_plasma[i].h2p,&matriz_plasma[i].carga);
    }
    fclose(dat);
    if(matriz_plasma[1].electrones==0){
        for(i=1;i<=(int)(nelectrones/res);i++){
            ni = (int)((alea())*1407)+1;
            if(ni==1408)ni=1407;
            matriz_plasma[ni].electrones+=res;
        }
        for(i=1;i<=(int)(nh20/res);i++){
            ni = (int)((alea())*1407)+1;
            if(ni==1408)ni=1407;
            matriz_plasma[ni].h20+=res;
        }
        for(i=1;i<=(int)(nhp/res);i++){
            ni = (int)((alea())*1407)+1;
            if(ni==1408)ni=1407;
            matriz_plasma[ni].hp+=res;
        }
        for(i=1;i<=(int)(nh2p/res);i++){
            ni = (int)((alea())*1407)+1;
            if(ni==1408)ni=1407;
            matriz_plasma[ni].h2p+=res;
        }
        ni = (int)((alea())*1407)+1;
        if(ni==1408)ni=1407;
        matriz_plasma[ni].electrones+=nelectrones%res;
        ni = (int)((alea())*1407)+1;
        if(ni==1408)ni=1407;
        matriz_plasma[ni].h20+=nh20%res;
        ni = (int)((alea())*1407)+1;
        if(ni==1408)ni=1407;
        matriz_plasma[ni].hp+=nhp%res;
        ni = (int)((alea())*1407)+1;
        if(ni==1408)ni=1407;
        matriz_plasma[ni].h2p+=nh2p%res;

        for(i=1;i<=1407;i++){
            calc_carga(i);
        }
    }
    actu_salida();
}
////////////////////////////////////////////////////////////////////////////////////
void mover_particulas(void){
    int nale;
    float tale;
    n1 = int(nceldas*alea())+1;
    if(n1==nceldas+1)n1=nceldas;
    n2 = int(nceldas*alea())+1;
    if(n2==nceldas+1)n2=nceldas;
    while(n2==n1){
        n2 = int(nceldas*alea())+1;
        if(n2==nceldas+1)n2=nceldas;
    }
    //if(p<1){
        tale = alea();
        nale = int(alea()*(4000000-0))+0;
        matriz_plasma[n1i] = matriz_plasma[n1];
        matriz_plasma[n2i] = matriz_plasma[n2];

        /*printf("\nn1: %i n2: %i nale: %i",n1,n2,nale);
        imprimir_celda_plasma(n1);
        imprimir_celda_plasma(n2);*/

        if(tale<0.25){
            tipo = 1;
            matriz_plasma[n1].electrones +=  nale;
            matriz_plasma[n2].electrones -=  nale;
            if(matriz_plasma[n2].electrones<0) rechazo=1;
        }
        else if(tale<0.5){
            tipo = 3;
            matriz_plasma[n1].h20 +=  nale;
            matriz_plasma[n2].h20 -=  nale;
            if(matriz_plasma[n2].h20<0) rechazo=1;
        }
        else if(tale<0.75){
            tipo = 2;
            matriz_plasma[n1].h2p +=  nale;
            matriz_plasma[n2].h2p -=  nale;
            if(matriz_plasma[n2].h2p<0) rechazo=1;
        }
        else{
            tipo = 2;
            matriz_plasma[n1].hp +=  nale;
            matriz_plasma[n2].hp -=  nale;
            if(matriz_plasma[n2].hp<0) rechazo=1;
        }
        calc_carga(n1);
        calc_carga(n2);
    //}

        /*imprimir_celda_plasma(n1);
        imprimir_celda_plasma(n2);
        getchar();*/

    /*imprimir_celda_plasma(n1i);
    imprimir_celda_plasma(n2i);
    imprimir_celda_plasma(n1);
    imprimir_celda_plasma(n2);
    getchar();*/
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda_plasma(int a){
    printf("\nCelda: %i electrones: %i h20: %i h2p: %i hp: %i carga: %i",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,matriz_plasma[a].h2p,matriz_plasma[a].hp,matriz_plasma[a].carga);
}
////////////////////////////////////////////////////////////////////////////////////
float distanciacelda(int a,int b){
    float dist,dist2,xx,yy,pp;

    dist = sqrt(pow(matriz_plasma[a].x-matriz_plasma[b].x,2)+pow(matriz_plasma[a].y-matriz_plasma[b].y,2));

    if(fabs(matriz_plasma[a].phi-matriz_plasma[b].phi)>pi*0.125){
        pp = matriz_plasma[b].phi + signo(matriz_plasma[a].phi-matriz_plasma[b].phi)*pi*0.25;
        xx = matriz_plasma[b].rho*cos(pp);
        yy = matriz_plasma[b].rho*sin(pp);
    }
    else{
        return(dist);
    }
    dist2 = sqrt(pow(matriz_plasma[a].x-xx,2)+pow(matriz_plasma[a].y-yy,2));
    if(dist<dist2){
        return(dist);
    }
    else{
        return(dist2);
    }
}
////////////////////////////////////////////////////////////////////////////////////
int signo(float a){
    int sign;
    if(a>=0){
        sign = 1;
    }
    else{
        sign = -1;
    }
    return(sign);
}
////////////////////////////////////////////////////////////////////////////////////
float de_plasma(void){
    int i,j,k,k1;
    int contadorvec_i=1, contadorvec_f=1;
    float eicoul = 0, efcoul = 0, eisig = 0, efsig = 0, elai = 0, elaf = 0;
    float d, z, r1, r2;
    float dem;

    for(i=1;i<=nceldas;i++){
        if(matriz_plasma[i].carga!=0){
            if(i!=n1){
                d = distanciacelda(n1i,i);
                eicoul += (qe*qe*matriz_plasma[n1i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                d = distanciacelda(n1,i);
                efcoul += (qe*qe*matriz_plasma[n1].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                /*printf("\n n1 eic: %e efc: %e",(qe*qe*matriz_plasma[n1i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc),
                (qe*qe*matriz_plasma[n1].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc));*/
            }
            if(i!=n2){
                d = distanciacelda(n2i,i);
                eicoul += (qe*qe*matriz_plasma[n2i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                d = distanciacelda(n2,i);
                efcoul += (qe*qe*matriz_plasma[n2].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                /*printf("\n n2 eic: %e efc: %e",(qe*qe*matriz_plasma[n2i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc),
                (qe*qe*matriz_plasma[n2].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc));*/
            }
            //getchar();
        }
    }
    dem = efcoul - eicoul;
    //printf("\n plasma eic: %e efc: %e de: %e",eicoul,efcoul,dem);
    //getchar();
    return(dem);
}
////////////////////////////////////////////////////////////////////////////////////
void metropolis_plasma(int a){
    float zeta, argexp, emet;
    int i;
    zeta = alea();
    if(a==1){
        argexp = -de_plasma()/(kb*tempee);
    }
    else if(a==2){
        argexp = -de_plasma()/(kb*tempei);
    }
    else{
        argexp = -de_plasma()/(kb*tempee);
    }

    if(argexp>=100)
    {
        emet = 2.0;
    }
    else
    {
        emet = exp(argexp);
    }

    if(emet>=zeta){
        c_mova++;
    }
    else{
        matriz_plasma[n1] = matriz_plasma[n1i];
        matriz_plasma[n2] = matriz_plasma[n2i];
        rechazo_met_mov++;
    }
}
////////////////////////////////////////////////////////////////////////////////////
void actu_salida(void){
	int i;

    sprintf(salidac,"datos/posiciones.dat");
    dat=fopen(salidac,"w");
    for(i=1;i<=1407;i++){
        if(matriz_plasma[i].carga == 0){
            fprintf(dat,"%4.0i\t%9.0i\t%9.0i\t%9.0i\t%9.0i\t        0\n",i,matriz_plasma[i].electrones,matriz_plasma[i].h20,matriz_plasma[i].hp,matriz_plasma[i].h2p);
        }
        else{
            fprintf(dat,"%4.0i\t%9.0i\t%9.0i\t%9.0i\t%9.0i\t%9.0i\n",i,matriz_plasma[i].electrones,matriz_plasma[i].h20,matriz_plasma[i].hp,matriz_plasma[i].h2p,matriz_plasma[i].carga);
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/electron.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=1407;i++){
        //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
        fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
    }
    fclose(dat);

    sprintf(salidac,"datos/h2+.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h2+\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=1407;i++){
        //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].h2p);
        fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].h2p);
    }
    fclose(dat);

    sprintf(salidac,"datos/h20.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h20\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=1407;i++){
        //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].h20);
        fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].h20);
    }
    fclose(dat);

    sprintf(salidac,"datos/hp.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"hp\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=1407;i++){
        //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].hp);
        fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].hp);
    }
    fclose(dat);

    sprintf(salidac,"datos/todas.dat");
    dat=fopen(salidac,"w");
    fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
    for(i=1;i<=1407;i++){
        fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones+matriz_plasma[i].h2p+matriz_plasma[i].h20+matriz_plasma[i].hp);
    }
    fclose(dat);

    sprintf(salidac,"datos/carga.dat");
    dat=fopen(salidac,"w");
    fprintf(dat, "\"x\", \"y\", \"carga\"\n");
    for(i=1;i<=1407;i++){
        fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].carga);
    }
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
float energia(void){
    int i,j;
    float energia_total=0,d;
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=nceldas;j++){
            if((matriz_plasma[i].carga!=0)&&(matriz_plasma[j].carga!=0)){
                if(i!=j){
                    d = distanciacelda(i,j);
                    energia_total += (qe*qe*matriz_plasma[i].carga*matriz_plasma[j].carga)/(4*pi*epce*epsi*d*esc);
                }
            }
        }
    }
    return(energia_total/2);
}
////////////////////////////////////////////////////////////////////////////////////
