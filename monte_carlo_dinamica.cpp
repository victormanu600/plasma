#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <tgmath.h>
#include <windows.h>
#include <vector>

using namespace std;


#define MAXPART 1000
#define pi 3.1415926535897932384
#define epce 8.8541878176e-12
#define epsi 1//78.5
#define muce 1.25663706144e-6
#define qe 1.60219e-19
#define kb 1.3806488e-23
#define na 6.02217e23
const int maxl = 20;
const int maxw = 20;

float alea(void);
float alea_f(float a, float b);
int alea_i(int a, int b);
void leer_datos_iniciales(void);
void leer_campo_magnetico(void);
void imprimir_datos_iniciales(void);
void imprimir_celda_plasma(int a);
void imprimir_celda_plasma_rec(int a);
void condiciones_iniciales(void);
void arreglo_inicial(void);
void crear_matriz_plasma(void);
void crear_matriz_plasma_rec(void);
void calc_carga(int a);

float distanciacelda(int a, int b);
float distancianormal(int a, int b);
float dist_ima(int a, int b, int c);
int signo(float a);
void mover_particulas(void);
void metropolis_plasma(int a);
void dinamica(void);

float de_plasma(void);
float calcular_de_mov(void);
void gdr_plasma(void);
void actu_salida(void);
void salida_prom(void);
float energia(void);
float ecoulomb(int a, int b, int c);
void beta(void);

void prueba_dinamica(void);
float aaadist_normal(float mu, float sigma);
void aaaprueba_dist(void);
void aaprueba(double &a);

vector<float> vect(float a[3]);
float punto(vector<float> a, vector<float> b);
vector<float> cruz(vector<float> a, vector<float> b);
void print_vector(vector<float> a);
int celda_dir(int a);
int celda_dir_rec(int a);
int celda_dir_rec_error(int a);
int celda_dir_rec_hp(int a);
int celda_dir_rec_error_hp(int a);
int celda_dir_rec_neg(int a);

void printposiciones(int a);

//////////////////////////////Constantes
int pasos, actu, terma, reso;
float esc, densidad, volumen, nelectronesr, tempee, tempei;
long long int nh20, nhp, nh2p, nelectrones, npart[4];
float R;
const int R_reso = 60*1*2+1;                   //Primer numero es R, segundo es reso, 2 es por que es de -R a R y +1 para comenzar los arreglos desde 1.
float me, mh20, mh2p, mhp, masa[4], carga[4];
///////////////////////////////////////////////////////////////////////Variables globales
int p=0, tipo, ienergia;
int rechazo;
int nceldas, n1, n2, n1i, n2i;                                //numero de celda Nueva/Vieja Inicial/Final
char cero[] = "         0";
///////////////////////////////////////////////////////////////////////Contadores
int c_mov=0,c_mova=1;
int c_uno=0,c_unoa=0,c_dos=0,c_dosa=0,c_tres=0,c_tresa=0;
int rechazo_neg=0, rechazo_met_mov=0;

long long int gelectron[11260], gh20[11260], gh2p[11260], ghp[11260], nprom[4][11260];

struct smatriz_plasma{
	float rho,phi,x,y;
    float vxe, vye, vxh2p, vyh2p, vxhp, vyhp, vxh20, vyh20;
    float vx[4], vy[4];
    float ve, vhp, vh2p, v[4];
	int carga;
	float t[3];
	long long int h20,hp,h2p,electrones,npart[4];
//}matriz_plasma[8*61*50+1];
//}matriz_plasma[7*61*60+1];
//}matriz_plasma[12346+10000];
}matriz_plasma[R_reso*R_reso+1+2];
//}matriz_plasma[3000*3000+1];

struct smatriz_plasma2{
	float rho,phi,x,y;
    float vxe, vye, vxh2p, vyh2p, vxhp, vyhp, vxh20, vyh20;
    float vx[4], vy[4];
    float ve, vhp, vh2p, v[4];
	int carga;
	float t[3];
	long long int h20,hp,h2p,electrones,npart[4];
}matriz_plasma2[R_reso*R_reso+1+2];
//}matriz_plasma2[2000*2000+1];

float ndt[2000*2000+1]={0};

/*struct smatriz_plasma2{
	float rho,phi,x,y;
    float vxe, vye, vxh2p, vyh2p, vxhp, vyhp, vxh20, vyh20;
    float ve, vhp, vh2p;
	int carga;
	float t[3];
	long long int h20,hp,h2p,electrones;
    float k1x, k1y, k1z, k2x, k2y, k2z;
//}matriz_plasma[8*61*50+1];
}matriz_plasma2[7*61*61+1];*/

//}matriz_plasma[8*61*50+1];

struct sarreglo{                        //pongo el +1 porque me gusta empezar los ciclos desde i,j = 1
    float z;
    float r[36+1+24];                      //36 valores posibles para r [0.0,3.5] de 0.1 en 0.1
    float br[36+1+24];
    float bz[36+1+24];
    float absb[36+1+24];
    float A[36+1];
    float dbz[36+1];
    float dbr[36+1];
    float fieldindex[36+1];
}arreglo[371+1];                        //371 valores posibles para z [0.0,37.0] de 0.1 en 0.1

struct struct_pruebaa{
    float x,y;
}struct_prueba[10];

struct scampo_magnetico{
    float br, bz, absb;
}campo_magnetico[60+1];

//int nmatriz_plasma[60+1][7*60];
int nmatriz_plasma[R_reso+1][R_reso+1];
//int nmatriz_plasma[2000+1][2000+1]={0};
////////////////////////////////////////////////////////////////////////Variables para salida
FILE *dat, *dat2;
char salidac[10];

main(){
    int i,j,k;
    srand((unsigned)time(NULL));
    /*system("mkdir datos");
    system("mkdir datos/datitos");*/
    CreateDirectory ("datos", NULL);
    CreateDirectory ("datos/promedios", NULL);
    printf("R_reso*R_reso+1: %i\n",R_reso*R_reso+1);
    int histograaaama[11]={0};
    for(i=1;i<=1000000000;i++){
        if(i%100000==0)printf("\ri: %i\talea_i: %i",i,alea_i(1,10));
        histograaaama[ alea_i(1,10) ]++;
    }
    dat = fopen("histograma_aleai.txt","w");
    for(i=1;i<=10;i++){
        fprintf(dat,"%i\t%i\n",i,histograaaama[i]);
    }
    fclose(dat);
    getchar();
    //getchar();
    /*float numero = 2.5;
    int floor_x;
    floor_x = floor(numero);
    printf("flotante: %f floor(flotante): %i",numero,floor_x);
    getchar();*/
    /*aaaprueba_dist();
    printf("Ya quedo la prueba v:");
    double fx[2]={2,3};
    printf("ANTES fx_0:%f fx_0:%f",fx[0],fx[1]);
    aaprueba(fx[0]);
    aaprueba(fx[1]);
    printf("\nDESPUES fx_0:%f fx_0:%f",fx[0],fx[1]);
    getchar();
    getchar();*/
    /*bool aaa=false;
    if(aaa){
        printf("Que shooooow aaa es true");
    }else{
        printf("Que shooooow aaa es false");
    }*/
    /*nmatriz_plasma[1][1]=1;nmatriz_plasma[1][2]=2;nmatriz_plasma[2][1]=3;nmatriz_plasma[2][2]=4;
	matriz_plasma[1].x = 1;	matriz_plasma[1].y = 0;
    matriz_plasma[1].rho = sqrt(matriz_plasma[1].x*matriz_plasma[1].x + matriz_plasma[1].y*matriz_plasma[1].y);
    matriz_plasma[1].phi = atan2(matriz_plasma[1].y,matriz_plasma[1].x);
	matriz_plasma[1].vxe = 0; matriz_plasma[1].vye = 1;
    printf("\ncelda_dir 1,1: %i", celda_dir(1) );
    getchar();*/
    /*printf("INICIANDO PRUEBA DINAMICA\n");
    prueba_dinamica();
    printf("TERMINANDO PRUEBA DINAMICA\n");
    getchar();*/

    /*vector<float> A, B, AXB;
    float a[3] = {1.2,1.2,1.2}, b[3] = {0,0,1.2};
    A = vect(a);
    B = vect(b);
    AXB = cruz(A,B);
    //float A[3], B[3];
    //A[0] = 1; A[1] = 1; A[2] = 0;
    //B[0] = 3; B[1] = 1; B[2] = 0;
    //cout << "A:\t" << A[0] << "\t" << A[1] << "\t" << A[2] << endl;
    //cout << "B:\t" << B[0] << "\t" << B[1] << "\t" << B[2] << endl;
    //cout << "Producto punto entre A y B: " << punto(A,B);
    cout << "Size of A: " << A.size() << endl;
    cout << "A:\t"; print_vector(A);
    cout << "B:\t"; print_vector(B);
    cout << "A.B:\t" << punto(A,B) << endl;
    cout << "AxB:\t"; print_vector(AXB);

    return(0);*/

    //PONER NCELDAS AQUI ABAJO
    for(i=1;i<=R_reso*R_reso+1;i++){
        //matriz_plasma[i].electrones=matriz_plasma[i].h20=matriz_plasma[i].h2p=matriz_plasma[i].hp=matriz_plasma[i].carga = 0;
        matriz_plasma[i].carga=matriz_plasma[i].electrones=matriz_plasma[i].h20=matriz_plasma[i].h2p=matriz_plasma[i].hp=0;
        matriz_plasma[i].phi=matriz_plasma[i].rho=matriz_plasma[i].ve=matriz_plasma[i].vh2p=matriz_plasma[i].vhp=matriz_plasma[i].vxe=matriz_plasma[i].vxh20=matriz_plasma[i].vxh2p=matriz_plasma[i].vxhp=matriz_plasma[i].vye=matriz_plasma[i].vyh20=matriz_plasma[i].vyh2p=matriz_plasma[i].vyhp=matriz_plasma[i].x=matriz_plasma[i].y=0;
        for(j=1;j<=4;j++){
            matriz_plasma[i].t[j]=0;
            matriz_plasma[i].v[j]=0;
            matriz_plasma[i].vy[j]=0;
            matriz_plasma[i].vx[j]=0;
        }
    }
	condiciones_iniciales();
    printf("int(R*reso): %i\n",int(R*reso));
	//crear_matriz_plasma();
	crear_matriz_plasma_rec();
	printf("Chequeo de celdas con rho>60mm\n");
	//getchar();
	//getchar();
	//getchar();
	for(i=1;i<=nceldas;i++){
        if(i%100==0)printf("\rCelda actual: %i",i);
        if(matriz_plasma[i].rho>60){
            imprimir_celda_plasma_rec(i);
            getchar();
        }
	}
	arreglo_inicial();
	leer_campo_magnetico();
	/*for(i=1;i<=nceldas;i++){
        printf("\nCELDA ORIGEN:\n");
        imprimir_celda_plasma_rec(i);
        printf("\nCELDA DESTINO:\n");
        if(celda_dir_rec(i)==0){
            imprimir_celda_plasma_rec(celda_dir_rec(i));printf("Que show");
            getchar();
        }
        else{
            imprimir_celda_plasma_rec(celda_dir_rec(i));
        }
	}
	getchar();*/

	/*for(i=1;i<=nceldas;i++){
        printf("\ncelda_dir: %i", celda_dir(i) );
	}
    getchar();*/
    for(i=1;i<=nceldas;i++){
        //matriz_plasma[nmatriz_plasma[ 62 ][ 62 ]].vxe = 0;
        //matriz_plasma[nmatriz_plasma[ 62 ][ 62 ]].vye = 1;
        matriz_plasma[i].vxe = 0;
        matriz_plasma[i].vye = 3;
        matriz_plasma[i].vxhp = 0;
        matriz_plasma[i].vyhp = 3;
        if(matriz_plasma[i].vye!=3){
            printf("\nvx: %f vy: %f\n",matriz_plasma[i].vxe,matriz_plasma[i].vye);
            getchar();
        }
    }
    /*matriz_plasma[nmatriz_plasma[ 801 ][ 901 ]].vxe = 1;
    matriz_plasma[nmatriz_plasma[ 801 ][ 901 ]].vye = 0;
    matriz_plasma[nmatriz_plasma[ 811 ][ 901 ]].vxe = -2;
    matriz_plasma[nmatriz_plasma[ 811 ][ 901 ]].vye = 0;*/
    /*printf("INICIANDO DINAMICA\n");
    sprintf(salidac,"datos/electron.dat",i);
    dat=fopen(salidac,"w");
    fprintf(dat, "#X\tY\tZ\n");
    sprintf(salidac,"datos/hp.dat",i);
    dat2=fopen(salidac,"w");
    fprintf(dat2, "#X\tY\tZ\n");
    fclose(dat);
    fclose(dat2);
    dat=fopen("datos/electron.dat","a");
    dat2=fopen("datos/hp.dat","a");
    for(i=1;i<=1000;i++){
        printf("\rdinamica_i: %i",i);
        dinamica();
        if(i%10==0){
            //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
            for(j=1;j<=nceldas;j++){
                if(matriz_plasma[j].electrones!=0){
                    fprintf(dat,"%f\t%f\t%I64d%i\n", matriz_plasma[j].x, matriz_plasma[j].y, matriz_plasma[j].electrones, i);
                }
                if(matriz_plasma[j].hp!=0){
                    fprintf(dat2,"%f\t%f\t%I64d%i\n", matriz_plasma[j].x, matriz_plasma[j].y, matriz_plasma[j].hp, i);
                }
            }
        }
        //getchar();
    }
    fclose(dat);
    fclose(dat2);
    printf("\nDESPUES DINAMICA\n");
    getchar();getchar();getchar();getchar();getchar();*/
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
        if(p>terma){
            gdr_plasma();
            if(p%actu==0){
                actu_salida();
                salida_prom();
                beta();
            }
        }
        if(p%1000==0||p<=10){
            if(p%actu==0){
                //printf("\rPaso: %i Acept. Mov: %1.5f Rechazo negativas: %1.5f Rechazo met: %1.5f Energia: %e",p,c_mova/(c_mov*1.0),rechazo_neg/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0),energia());
                printf("\rPaso: %i A. tipo1: %1.5f A. tipo2: %1.5f A. tipo3: %1.5f uno: %1.5f Energia: %e",p,c_unoa/(c_uno*1.0),c_dosa/(c_dos*1.0),c_tresa/(c_tres*1.0),(c_unoa+c_dosa+c_tresa+rechazo_met_mov+rechazo_neg)/(c_mov*1.0),energia());
            }
            else{
                //printf("\rPaso: %i Acept. Mov: %1.5f Rechazo negativas: %1.5f Rechazo met: %1.5f",p,c_mova/(c_mov*1.0),rechazo_neg/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0));
                printf("\rPaso: %i A. tipo1: %1.5f A. tipo2: %1.5f A. tipo3: %1.5f uno: %1.5f ",p,c_unoa/(c_uno*1.0),c_dosa/(c_dos*1.0),c_tresa/(c_tres*1.0),(c_unoa+c_dosa+c_tresa+rechazo_met_mov+rechazo_neg)/(c_mov*1.0));
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
float alea_f(float a, float b){
    return( ( a + ( (float)rand()/RAND_MAX )*( b - a ) ) );
}
////////////////////////////////////////////////////////////////////////////////////
int alea_i(int a, int b){
    return( ( a + int( ( (float)rand()/(RAND_MAX+1) )*( 1 + b - a ) ) ) );
}
////////////////////////////////////////////////////////////////////////////////////
void leer_datos_iniciales(){
    float dummy;
    dat=fopen("entrada.txt","r");
    fscanf(dat,"Numero de pasos: %i\n", &pasos);
    fscanf(dat,"Actualizacion: %i\n", &actu);
    fscanf(dat,"Termalizacion: %i\n", &terma);
    fscanf(dat,"Resolucion: %i\n", &reso);
    fscanf(dat,"Masas (kg): %f, %f, %f, %f\n",&me,&mh20,&mh2p,&mhp);
    fscanf(dat,"R (mm): %f\n",&R);
    fscanf(dat,"Temperatura electrones (eV): %f\n",&tempee);
    fscanf(dat,"Temperatura iones (eV): %f\n",&tempei);
    fclose(dat);
    dat = fopen("datos/energia.dat","r");
    while(fscanf(dat,"\n%i\t%f",&ienergia,&dummy)!=EOF);
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_datos_iniciales(){
    printf("Numero de pasos: %i\n", pasos);
    printf("Actualizacion: %i\n", actu);
    printf("Termalizacion: %i\n", terma);
    printf("Resolucion: %i\n", reso);
    printf("Masas (kg): %e, %e, %e, %e\n",me,mh20,mh2p,mhp);
    printf("R (mm): %f\n",R);
    printf("Temperatura electrones (K): %f\n",tempee);
    printf("Temperatura iones (K): %f\n\n",tempei);
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_iniciales(){
    printf("Condiciones iniciales\n");
    leer_datos_iniciales();
    esc = 1e-3;//diam*1e-9;

    densidad = 2e19;
    //volumen = pi*R*esc*R*esc*0.125*1e-6;
    volumen = pi*R*esc*R*esc*1e-6;
    nelectronesr = (densidad*volumen);
    nelectrones = nelectronesr;
    nh20 = 9*nelectrones;
    nhp = 0.8*nelectrones;
    nh2p = 0.2*nelectrones+1;
    printf("\ndensidad: %e volumen: %e nnr: %e nelectrones: %I64d",densidad,volumen,nelectronesr,nelectrones);
    printf("\nnh20: %I64d nhp: %I64d nh2p: %I64d cargatotal: %I64d\n",nh20,nhp,nh2p,nh2p+nhp-nelectrones);
    //getchar();

    tempei = tempei*qe/kb;
    tempee = tempee*qe/kb;

    printf("te: %f ti: %f\n",tempee,tempei);

    //printf("\ntempei: %f tempee: %f\n", tempei, tempee);
    //getchar();

    imprimir_datos_iniciales();
}

////////////////////////////////////////////////////////////////////////////////////
void crear_matriz_plasma(){
    int i,j,k,contadorm,contadorvec=0;
    float vi, fvi;
    float fraccion = 0.25;

    contadorm=0;
    for(k=1;k<=2/fraccion;k++){
    //for(k=1;k<=1;k++){
        contadorm++;
        nmatriz_plasma[1][k] = contadorm;
        //matriz_plasma[1].rho = 1/sqrt(2);
        matriz_plasma[contadorm].rho = 1/sqrt(2);
        matriz_plasma[contadorm].phi = pi/8+(k-1)*pi*fraccion;
        printf("contador: %d phi: %f\n",contadorm,matriz_plasma[contadorm].phi);
        matriz_plasma[contadorm].x = matriz_plasma[contadorm].rho*cos(matriz_plasma[contadorm].phi);
        matriz_plasma[contadorm].y = matriz_plasma[contadorm].rho*sin(matriz_plasma[contadorm].phi);

        /*for(i=2;i<=R;i++){
            for(j=1;j<=(int)((pi*i)*fraccion);j++){
                contadorm++;
                nmatriz_plasma[i][j+(k-1)*(int)((pi*i)*fraccion)]=contadorm;
                matriz_plasma[contadorm].rho = sqrt(pow(i,2)+pow(i-1,2))/sqrt(2);
                //matriz_plasma[contadorm].rho = i-0.5;
                matriz_plasma[contadorm].phi = pi/(8*((int)((pi*i)*fraccion)))+(j-1)*(pi/(4*((int)((pi*i)*fraccion)))) + (k-1)*pi*fraccion;
                matriz_plasma[contadorm].x = matriz_plasma[contadorm].rho*cos(matriz_plasma[contadorm].phi);
                matriz_plasma[contadorm].y = matriz_plasma[contadorm].rho*sin(matriz_plasma[contadorm].phi);
            }
        }*/
        for(i=2;i<=R;i++){
            for(j=1;j<=(int)((pi*i)*fraccion);j++){
                contadorm++;
                nmatriz_plasma[i][j+(k-1)*(int)((pi*i)*fraccion)]=contadorm;
                matriz_plasma[contadorm].rho = sqrt(pow(i,2)+pow(i-1,2))/sqrt(2);
                //matriz_plasma[contadorm].rho = i-0.5;
                matriz_plasma[contadorm].phi = 0.5*fraccion*pi/((int)((pi*i)*fraccion))+(j-1)*(fraccion*pi/(((int)((pi*i)*fraccion)))) + (k-1)*pi*fraccion;
                matriz_plasma[contadorm].x = matriz_plasma[contadorm].rho*cos(matriz_plasma[contadorm].phi);
                matriz_plasma[contadorm].y = matriz_plasma[contadorm].rho*sin(matriz_plasma[contadorm].phi);
                /*do{
                    vi=alea_f(-4*sqrt((kb*tempee)/me), 4*sqrt((kb*tempee)/me) );
                    fvi=exp(-vi*vi*me/(2*kb*tempee))*pow(me/(2*pi*kb*tempee), 1.0/2.0);
                }while(alea()*pow(me/(2*pi*kb*tempee), 1.0/2.0)>fvi);
                matriz_plasma[contadorm].vxe = vi;*/
            }
        }
    }
    nceldas=contadorm;
    n1i = nceldas+1;
    n2i = nceldas+2;
    for(i=1;i<=nceldas;i++){
        do{
            vi=alea_f(-4*sqrt((kb*tempee)/me), 4*sqrt((kb*tempee)/me) );
            fvi=exp(-vi*vi*me/(2*kb*tempee))*pow(me/(2*pi*kb*tempee), 1.0/2.0);
        }while(alea()*pow(me/(2*pi*kb*tempee), 1.0/2.0)>fvi);
        matriz_plasma[i].vxe = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempee)/me), 4*sqrt((kb*tempee)/me) );
            fvi=exp(-vi*vi*me/(2*kb*tempee))*pow(me/(2*pi*kb*tempee), 1.0/2.0);
        }while(alea()*pow(me/(2*pi*kb*tempee), 1.0/2.0)>fvi);
        matriz_plasma[i].vye = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/mhp), 4*sqrt((kb*tempei)/mhp) );
            fvi=exp(-vi*vi*mhp/(2*kb*tempei))*pow(mhp/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(mhp/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vxhp = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/mhp), 4*sqrt((kb*tempei)/mhp) );
            fvi=exp(-vi*vi*mhp/(2*kb*tempei))*pow(mhp/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(mhp/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vyhp = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/mh2p), 4*sqrt((kb*tempei)/mh2p) );
            fvi=exp(-vi*vi*mh2p/(2*kb*tempei))*pow(mh2p/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(mh2p/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vxh2p = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/mh2p), 4*sqrt((kb*tempei)/mh2p) );
            fvi=exp(-vi*vi*mh2p/(2*kb*tempei))*pow(mh2p/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(mh2p/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vyh2p = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/mh20), 4*sqrt((kb*tempei)/mh20) );
            fvi=exp(-vi*vi*mh20/(2*kb*tempei))*pow(mh20/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(mh20/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vxh20 = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/mh20), 4*sqrt((kb*tempei)/mh20) );
            fvi=exp(-vi*vi*mh20/(2*kb*tempei))*pow(mh20/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(mh20/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vyh20 = vi;

        matriz_plasma[i].ve = sqrt( matriz_plasma[i].vxe*matriz_plasma[i].vxe + matriz_plasma[i].vye*matriz_plasma[i].vye );
        matriz_plasma[i].vhp = sqrt( matriz_plasma[i].vxhp*matriz_plasma[i].vxhp + matriz_plasma[i].vyhp*matriz_plasma[i].vyhp );
        matriz_plasma[i].vh2p = sqrt( matriz_plasma[i].vxh2p*matriz_plasma[i].vxh2p + matriz_plasma[i].vyh2p*matriz_plasma[i].vyh2p );
    }

    printf("nceldas: %i\tn1i: %i\tn2i: %i",nceldas,n1i,n2i);
    dat = fopen("celdas.dat","w");
    fprintf(dat,"#X Y Z\n");
    for(i=1;i<=nceldas;i++){
        fprintf(dat,"%e\t%e\t%e\n",matriz_plasma[i].x,matriz_plasma[i].y, sqrt( matriz_plasma[i].x*matriz_plasma[i].x + matriz_plasma[i].y*matriz_plasma[i].y ) );
    }
    fclose(dat);
    dat = fopen("celdas2.dat","w");
    for(k=1;k<=2/fraccion;k++){
        fprintf(dat,"%d\t%d\t%d\n",1,k,nmatriz_plasma[1][k]);
        for(i=1;i<=R;i++){
            for(j=1;j<=(int)((pi*i)*0.25);j++){
                fprintf(dat,"%d\t%d\t%d\n",i,j+(k-1)*(int)((pi*i)*0.25),nmatriz_plasma[i][j+(k-1)*(int)((pi*i)*0.25)]);
            }
        }
    }
    fclose(dat);

    long int clases[201]={0};
    for(i=1;i<=nceldas;i++){
        clases[ int( 200*(matriz_plasma[i].vxe+4*sqrt((kb*tempee)/me))/(8*sqrt((kb*tempee)/me)) )+1]++;
    }
    printf("Que show");
    //getchar();
    dat = fopen("histograma.dat","w");
    fprintf(dat,"#X\tY\n");
    for(i=1;i<=200;i++){
        fprintf(dat,"%e\t%i\n", -4*sqrt((kb*tempee)/me)+ (i-1)*(8*sqrt((kb*tempee)/me))/(200) , clases[i] );
    }
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void crear_matriz_plasma_rec(){
    int i,j,k,contadorm=0,contadorm2=11291,contadorvec=0;
    float vi, fvi;
    printf("Crear matriz plasma rec\n");
    printf("Asignando posiciones\n");
    /*printf("sizeof(plasma[1]): %i\tsizeof(plasma[2]): %i",sizeof(matriz_plasma[1]),sizeof(matriz_plasma[2]));
    printf("sizeof(plasma): %i\tplasma/plasma[1]: %i\n",sizeof(matriz_plasma),sizeof(matriz_plasma)/sizeof(matriz_plasma[1]));
    printf("sizeof(matriz[1][1]): %i\tsizeof(matriz): %i",sizeof(nmatriz_plasma[1][1]),sizeof(nmatriz_plasma));*/
    //getchar();
    for(i=-int(R*reso);i<=int(R*reso);i++){
        for(j=-int(R*reso);j<=int(R*reso);j++){
            if(i*i+j*j<=int(R*reso)*int(R*reso)){
                contadorm++;
                nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1]=contadorm;
                matriz_plasma[contadorm].x = (1.0*(i+0.1))/reso;
                matriz_plasma[contadorm].y = (1.0*(j+0.1))/reso;
                matriz_plasma[contadorm].rho = sqrt(i*i+j*j)/reso;
                //matriz_plasma[contadorm].rho = sqrt((matriz_plasma[contadorm].x*reso-0.1)*(matriz_plasma[contadorm].x*reso-0.1)+(matriz_plasma[contadorm].y*reso-0.1)*(matriz_plasma[contadorm].y*reso-0.1));
                matriz_plasma[contadorm].phi = (atan2(matriz_plasma[contadorm].y,matriz_plasma[contadorm].x)+((atan2(matriz_plasma[contadorm].y,matriz_plasma[contadorm].x))<0?2*pi:0.0));
            }
            else{
                contadorm2++;
                nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1]=contadorm2;
                matriz_plasma[contadorm2].x = (1.0*(i+0.1))/reso;
                matriz_plasma[contadorm2].y = (1.0*(j+0.1))/reso;
                matriz_plasma[contadorm2].rho = sqrt(i*i+j*j)/reso;
                //matriz_plasma[contadorm2].rho = sqrt(matriz_plasma[contadorm2].x*matriz_plasma[contadorm2].x+matriz_plasma[contadorm2].y*matriz_plasma[contadorm2].y);
                matriz_plasma[contadorm2].phi = (atan2(matriz_plasma[contadorm2].y,matriz_plasma[contadorm2].x)+((atan2(matriz_plasma[contadorm2].y,matriz_plasma[contadorm2].x))<0?2*pi:0.0));
            }
        }
    }
    /*int int_maxx = 600;
    for(i=-int_maxx;i<=int_maxx;i++){
        for(j=-int_maxx;j<=int_maxx;j++){
            if(i+int(R*reso)+1>1201||j+int(R*reso)+1>1201){
                printf("i: %i\tj: %i",i+int(R*reso)+1,j+int(R*reso)+1);
                getchar();
            }
            if(i*i+j*j<=int_maxx*int_maxx){
                contadorm++;
                nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1]=contadorm;
                matriz_plasma[contadorm].x = 0.1*(i+0.1);
                matriz_plasma[contadorm].y = 0.1*(j+0.1);
                matriz_plasma[contadorm].rho = sqrt((matriz_plasma[contadorm].x-0.01)*(matriz_plasma[contadorm].x-0.01)+(matriz_plasma[contadorm].y-0.01)*(matriz_plasma[contadorm].y-0.01));
                matriz_plasma[contadorm].phi = (atan2(matriz_plasma[contadorm].y,matriz_plasma[contadorm].x)+((atan2(matriz_plasma[contadorm].y,matriz_plasma[contadorm].x))<0?2*pi:0.0));
            }
            else{
                contadorm2++;
                nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1]=contadorm2;
                matriz_plasma[contadorm2].x = 0.1*(i+0.1);
                matriz_plasma[contadorm2].y = 0.1*(j+0.1);
                matriz_plasma[contadorm2].rho = sqrt((matriz_plasma[contadorm2].x-0.01)*(matriz_plasma[contadorm2].x-0.01)+(matriz_plasma[contadorm2].y-0.01)*(matriz_plasma[contadorm2].y-0.01));
                matriz_plasma[contadorm2].phi = (atan2(matriz_plasma[contadorm2].y,matriz_plasma[contadorm2].x)+((atan2(matriz_plasma[contadorm2].y,matriz_plasma[contadorm2].x))<0?2*pi:0.0));
            }
        }
    }*/
    nceldas=contadorm;
    n1i = nceldas+1;
    n2i = nceldas+2;
    printf("\nnceldas: %i\tn1i: %i\tn2i: %i\n",nceldas,n1i,n2i);
    printf("Asignando velocidades\n");
    for(i=1;i<=nceldas;i++){
        if(i%1000==0)printf("\rCelda actual: %i",i);
        matriz_plasma[i].vxe = aaadist_normal(0.0,sqrt(kb*tempee/me));
        matriz_plasma[i].vye = aaadist_normal(0.0,sqrt(kb*tempee/me));
        matriz_plasma[i].vxhp = aaadist_normal(0.0,sqrt(kb*tempei/mhp));
        matriz_plasma[i].vyhp = aaadist_normal(0.0,sqrt(kb*tempei/mhp));
        matriz_plasma[i].vxh2p = aaadist_normal(0.0,sqrt(kb*tempei/mh2p));
        matriz_plasma[i].vyh2p = aaadist_normal(0.0,sqrt(kb*tempei/mh2p));
        matriz_plasma[i].vxh20 = aaadist_normal(0.0,sqrt(kb*tempei/mh20));
        matriz_plasma[i].vyh20 = aaadist_normal(0.0,sqrt(kb*tempei/mh20));

        matriz_plasma[i].ve = sqrt( matriz_plasma[i].vxe*matriz_plasma[i].vxe + matriz_plasma[i].vye*matriz_plasma[i].vye );
        matriz_plasma[i].vhp = sqrt( matriz_plasma[i].vxhp*matriz_plasma[i].vxhp + matriz_plasma[i].vyhp*matriz_plasma[i].vyhp );
        matriz_plasma[i].vh2p = sqrt( matriz_plasma[i].vxh2p*matriz_plasma[i].vxh2p + matriz_plasma[i].vyh2p*matriz_plasma[i].vyh2p );
    }
    /*for(i=1;i<=nceldas;i++){
        if(i%1000==0)printf("\rCelda actual: %i",i);
        matriz_plasma[i].vxe = aaadist_normal(0.0,1);
        matriz_plasma[i].vye = aaadist_normal(0.0,1);
        matriz_plasma[i].vxhp = aaadist_normal(0.0,1);
        matriz_plasma[i].vyhp = aaadist_normal(0.0,1);
        matriz_plasma[i].vxh2p = aaadist_normal(0.0,1);
        matriz_plasma[i].vyh2p = aaadist_normal(0.0,1);
        matriz_plasma[i].vxh20 = aaadist_normal(0.0,1);
        matriz_plasma[i].vyh20 = aaadist_normal(0.0,1);

        matriz_plasma[i].ve = sqrt( matriz_plasma[i].vxe*matriz_plasma[i].vxe + matriz_plasma[i].vye*matriz_plasma[i].vye );
        matriz_plasma[i].vhp = sqrt( matriz_plasma[i].vxhp*matriz_plasma[i].vxhp + matriz_plasma[i].vyhp*matriz_plasma[i].vyhp );
        matriz_plasma[i].vh2p = sqrt( matriz_plasma[i].vxh2p*matriz_plasma[i].vxh2p + matriz_plasma[i].vyh2p*matriz_plasma[i].vyh2p );
    }*/
    printf("\rImprimiendo archivos de celda\n");
    dat = fopen("celdas.dat","w");
    fprintf(dat,"#X Y Z\n");
    for(i=1;i<=nceldas;i++){
        fprintf(dat,"%e\t%e\t%e\n",matriz_plasma[i].x,matriz_plasma[i].y, matriz_plasma[i].rho );
    }
    fclose(dat);
    dat = fopen("celdas2.dat","w");
    fprintf(dat,"i\tj\tx\ty\tncelda\n");
    for(i=1;i<=2*R*reso+1;i++){
        for(j=1;j<=2*R*reso+1;j++){
            fprintf(dat,"%d\t%d\t%f\t%f\t%d\n",i,j,matriz_plasma[nmatriz_plasma[i][j]].x,matriz_plasma[nmatriz_plasma[i][j]].y, nmatriz_plasma[i][j]);
        }
    }
    fclose(dat);

    printf("Haciendo histograma de velocidades\n");

    long int clases[201]={0};
    for(i=1;i<=nceldas;i++){
        clases[ int( 200*(matriz_plasma[i].vxe+5*sqrt((kb*tempee)/me))/(10*sqrt((kb*tempee)/me)) )+1]++;
    }
    //getchar();
    dat = fopen("histograma.dat","w");
    fprintf(dat,"#X\tY\n");
    for(i=1;i<=200;i++){
        fprintf(dat,"%e\t%i\n", -5*sqrt((kb*tempee)/me)+ (i-1)*(10*sqrt((kb*tempee)/me))/(200) , clases[i] );
    }
    fclose(dat);
    printf("Termino de hacer histograma de velocidades\n");
}
////////////////////////////////////////////////////////////////////////////////////
void calc_carga(int a){
    matriz_plasma[a].carga = matriz_plasma[a].h2p + matriz_plasma[a].hp - matriz_plasma[a].electrones;
}
////////////////////////////////////////////////////////////////////////////////////
void arreglo_inicial(){
    int i,j,ni,dummy,res=10000;
    long long int cargatotal=0,ee=0, h200=0, h2pp=0, hpp=0;

    printf("\rAsignando arreglo inicial\n");
    if(ienergia!=0){
        printf("Continuando posicion final de corrida anterior\n");
        dat = fopen("datos/posiciones.dat","r");
        for(i=1;i<=nceldas;i++){
            fscanf(dat,"\n%i\t%i\t%i\t%i\t%i\t%i",&dummy,&matriz_plasma[i].electrones,&matriz_plasma[i].h20,&matriz_plasma[i].hp,&matriz_plasma[i].h2p,&matriz_plasma[i].carga);
            if((matriz_plasma[i].electrones==0)||(matriz_plasma[i].h20==0)||(matriz_plasma[i].hp==0)||(matriz_plasma[i].h2p==0)||(matriz_plasma[i].carga==0)){
                //imprimir_celda_plasma(i);
            }
        }
        fclose(dat);
    }else{
        printf("Arreglo nuevo\n");
        for(i=1;i<=nceldas;i++){
            matriz_plasma[i].electrones = nelectrones/nceldas;
            matriz_plasma[i].hp = nhp/nceldas;
            matriz_plasma[i].h2p = nh2p/nceldas;
            matriz_plasma[i].h20 = nh20/nceldas;
        }
        for(i=1;i<=nelectrones%nceldas;i++){
            matriz_plasma[i].electrones += 1;
        }
        for(i=1;i<=nhp%nceldas;i++){
            matriz_plasma[i].hp += 1;
        }
        for(i=1;i<=nh2p%nceldas;i++){
            matriz_plasma[i].h2p += 1;
        }
        for(i=1;i<=nh20%nceldas;i++){
            matriz_plasma[i].h20 += 1;
        }
        printf("Calculando carga\n");
        for(i=1;i<=nceldas;i++){
            calc_carga(i);
        }
    }

    /*if(matriz_plasma[1].electrones==0){
        printf("\nAsignando electrones");
        for(i=1;i<=(int)(nelectrones/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].electrones+=res;
        }
        printf("\nAsignando moleculas de hidrogeno");
        for(i=1;i<=(int)(nh20/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].h20+=res;
        }
        printf("\nAsignando hidrones");
        for(i=1;i<=(int)(nhp/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].hp+=res;
        }
        printf("\nAsignando deuterones");
        for(i=1;i<=(int)(nh2p/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].h2p+=res;
        }
        ni = (int)((alea())*nceldas)+1;
        if(ni==nceldas+1)ni=nceldas;
        matriz_plasma[ni].electrones+=nelectrones%res;
        ni = (int)((alea())*nceldas)+1;
        if(ni==nceldas+1)ni=nceldas;
        matriz_plasma[ni].h20+=nh20%res;
        ni = (int)((alea())*nceldas)+1;
        if(ni==nceldas+1)ni=nceldas;
        matriz_plasma[ni].hp+=nhp%res;
        ni = (int)((alea())*nceldas)+1;
        if(ni==nceldas+1)ni=nceldas;
        matriz_plasma[ni].h2p+=nh2p%res;

        for(i=1;i<=nceldas;i++){
            calc_carga(i);
        }
    }*/
    /*for(i=1;i<=nceldas;i++){
        if(i==nmatriz_plasma[ int((R+20)*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].electrones=nelectrones/2.0;
            matriz_plasma[i].h2p=0;
            matriz_plasma[i].hp=0;
            matriz_plasma[i].h20=0;
        }
        else if(i==nmatriz_plasma[ int((R+15)*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].electrones=nelectrones/2.0;
            matriz_plasma[i].h2p=0;
            matriz_plasma[i].hp=0;
            matriz_plasma[i].h20=0;
        }else if(i==nmatriz_plasma[ int((R-30)*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].electrones=0;
            matriz_plasma[i].h2p=nh2p/2.0;
            matriz_plasma[i].hp=nhp/2.0;
            matriz_plasma[i].h20=nh20/2.0;
        }else if(i==nmatriz_plasma[ int((R-25)*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].electrones=0;
            matriz_plasma[i].h2p=nh2p/2.0;
            matriz_plasma[i].hp=nhp/2.0;
            matriz_plasma[i].h20=nh20/2.0;
        }else{
            matriz_plasma[i].electrones=0;
            matriz_plasma[i].h2p=0;
            matriz_plasma[i].hp=0;
            matriz_plasma[i].h20=0;
        }
        calc_carga(i);
    }*/
    /*if((matriz_plasma[1].electrones==0)&&(matriz_plasma[1].h20==0)&&(matriz_plasma[1].hp==0)&&(matriz_plasma[1].h2p==0)&&(matriz_plasma[1].carga==0)){
        for(i=1;i<=5;i++){
            for(j=1;j<=2;j++){
                matriz_plasma[ nmatriz_plasma[30+i][15+j] ].electrones = nelectrones/10;
                matriz_plasma[ nmatriz_plasma[30+i][15+j] ].h2p = nh2p/10;
                matriz_plasma[ nmatriz_plasma[30+i][15+j] ].hp = nhp/10;
                matriz_plasma[ nmatriz_plasma[30+i][15+j] ].h20 = nh20/10;
                cargatotal += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].h2p + matriz_plasma[ nmatriz_plasma[30+i][15+j] ].hp - matriz_plasma[ nmatriz_plasma[30+i][15+j] ].electrones;
                printf("\ni: %i j: %i electrones: %I64d h20: %I64d hp: %I64d h2p: %I64d carga: %I64d",i,j,nelectrones/10,nh20/10,nhp/10,nh2p/10,-nelectrones/10+nh2p/10+nhp/10);
                ee += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].electrones;
                h200 += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].h20;
                h2pp += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].h2p;
                hpp += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].hp;
                calc_carga( nmatriz_plasma[30+i][15+j] );
            }
        }
        matriz_plasma[500].electrones = 4;
        matriz_plasma[500].h20 = 6;
        matriz_plasma[500].hp = 5;
        matriz_plasma[500].h2p = 9;
        calc_carga(500);
        cargatotal += matriz_plasma[500].hp + matriz_plasma[500].h2p - matriz_plasma[500].electrones;
        printf("\nCargatotal: %I64d ee: %I64d h200: %I64d hpp: %I64d h2pp: %I64d",cargatotal,ee-nelectrones,h200-nh20,hpp-nhp,h2pp-nh2p);
    }*/
    sprintf(salidac,"datos/posiciones2.dat");
    dat=fopen(salidac,"w");
    for(i=1;i<=nceldas;i++){
        printposiciones(i);
    }
    fclose(dat);
    //printf("\nYa termino posiciones2");
    //getchar();
    p = terma;
    actu_salida();
    p = 0;
}
////////////////////////////////////////////////////////////////////////////////////
void mover_particulas(void){
    long long int nale;
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
        matriz_plasma[n1i] = matriz_plasma[n1];
        matriz_plasma[n2i] = matriz_plasma[n2];

        /*printf("\nn1: %i n2: %i nale: %i",n1,n2,nale);
        imprimir_celda_plasma(n1);
        imprimir_celda_plasma(n2);*/

        if(tale<0.25){
            tipo = 1;
            c_uno++;
            nale = alea()*matriz_plasma[n2].electrones;
            nale += 1;
            if(nale==matriz_plasma[n2].electrones+1)nale=matriz_plasma[n2].electrones;
            matriz_plasma[n1].electrones +=  nale;
            matriz_plasma[n2].electrones -=  nale;
            if((matriz_plasma[n2].electrones<0)||(matriz_plasma[n1].electrones<0)) rechazo=1;
        }
        else if(tale<0.5){
            tipo = 2;
            c_dos++;
            nale = alea()*matriz_plasma[n2].hp;
            nale += 1;
            if(nale==matriz_plasma[n2].hp+1)nale=matriz_plasma[n2].hp;
            matriz_plasma[n1].hp +=  nale;
            matriz_plasma[n2].hp -=  nale;
            if((matriz_plasma[n2].hp<0)||(matriz_plasma[n1].hp<0)) rechazo=1;
        }
        else if(tale<0.75){
            tipo = 2;
            c_dos++;
            nale = alea()*matriz_plasma[n2].h2p;
            nale += 1;
            if(nale==matriz_plasma[n2].h2p+1)nale=matriz_plasma[n2].h2p;
            matriz_plasma[n1].h2p +=  nale;
            matriz_plasma[n2].h2p -=  nale;
            if((matriz_plasma[n2].h2p<0)||(matriz_plasma[n1].h2p<0)) rechazo=1;
        }
        else{
            tipo = 3;
            c_tres++;
            nale = alea()*matriz_plasma[n2].h20;
            nale += 1;
            if(nale==matriz_plasma[n2].h20+1)nale=matriz_plasma[n2].h20;
            matriz_plasma[n1].h20 +=  nale;
            matriz_plasma[n2].h20 -=  nale;
            if((matriz_plasma[n2].h20<0)||(matriz_plasma[n1].h20<0)) rechazo=1;
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
    printf("\nCelda: %i rho: %f phi: %f\nelectrones: %I64d h20: %I64d hp: %I64d h2p: %I64d carga: %i",a,matriz_plasma[a].rho,matriz_plasma[a].phi,matriz_plasma[a].electrones,matriz_plasma[a].h20,matriz_plasma[a].hp,matriz_plasma[a].h2p,matriz_plasma[a].carga);
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda_plasma_rec(int a){
    printf("\nCelda: %i x: %f y: %f\nelectrones: %I64d h20: %I64d hp: %I64d h2p: %I64d carga: %i",a,matriz_plasma[a].x,matriz_plasma[a].y,matriz_plasma[a].electrones,matriz_plasma[a].h20,matriz_plasma[a].hp,matriz_plasma[a].h2p,matriz_plasma[a].carga);
}
////////////////////////////////////////////////////////////////////////////////////
float distanciacelda(int a, int b){
    float dist,dist2,xx,yy,pp;

    dist = sqrt(pow(matriz_plasma[a].x-matriz_plasma[b].x,2)+pow(matriz_plasma[a].y-matriz_plasma[b].y,2));
    dist2 = sqrt(pow(matriz_plasma[a].x-matriz_plasma[b].x,2)+pow(matriz_plasma[a].y-matriz_plasma[b].y,2));

    if(fabs(matriz_plasma[a].phi-matriz_plasma[b].phi)>pi*0.125){
        pp = matriz_plasma[b].phi + signo(matriz_plasma[a].phi-matriz_plasma[b].phi)*pi*0.25;
        xx = matriz_plasma[b].rho*cos(pp);
        yy = matriz_plasma[b].rho*sin(pp);
        dist2 = sqrt(pow(matriz_plasma[a].x-xx,2)+pow(matriz_plasma[a].y-yy,2));
    }
    else{
        return(dist);
    }
    if(dist<dist2){
        return(dist);
    }
    else{
        return(dist2);
    }
}
////////////////////////////////////////////////////////////////////////////////////
float distancianormal(int a, int b){
    float dist;
    dist = sqrt(pow(matriz_plasma[a].x-matriz_plasma[b].x,2)+pow(matriz_plasma[a].y-matriz_plasma[b].y,2));
    return(dist);
}
////////////////////////////////////////////////////////////////////////////////////
float dist_ima(int a, int b, int c){
    float dist,xx,yy,pp;

    pp = matriz_plasma[b].phi +(c-1)*pi*0.25;
    xx = matriz_plasma[b].rho*cos(pp);
    yy = matriz_plasma[b].rho*sin(pp);

    dist = sqrt(pow(matriz_plasma[a].x-xx,2)+pow(matriz_plasma[a].y-yy,2));
    return(dist);
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
    int i;
    float ei = 0, ef = 0, d, dem;
    for(i=1;i<=nceldas;i++){
        if(matriz_plasma[i].carga!=0){
            if((i!=n1)&&(i!=n2)){
                //d = distanciacelda(n1i,i);
                d = distancianormal(n1i,i);
                ei += (qe*qe*matriz_plasma[n1i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                //d = distanciacelda(n1,i);
                d = distancianormal(n1,i);
                ef += (qe*qe*matriz_plasma[n1].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                /*ei += 8*ecoulomb(n1i,i,1);
                ef += 8*ecoulomb(n1,i,1);*/
                //d = distanciacelda(n2i,i);
                d = distancianormal(n2i,i);
                ei += (qe*qe*matriz_plasma[n2i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                //d = distanciacelda(n2,i);
                d = distancianormal(n2,i);
                ef += (qe*qe*matriz_plasma[n2].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                /*ei += 8*ecoulomb(n2i,i,1);
                ef += 8*ecoulomb(n2,i,1);*/
            }
            //getchar();
        }
    }
    //d = distanciacelda(n1i,n2i);
    d = distancianormal(n1i,n2i);
    ei += (qe*qe*matriz_plasma[n1i].carga*matriz_plasma[n2i].carga)/(4*pi*epce*epsi*d*esc);
    //ei += 8*ecoulomb(n1i,n2i,1);
    //d = distanciacelda(n1,n2);
    d = distancianormal(n1,n2);
    ef += (qe*qe*matriz_plasma[n1].carga*matriz_plasma[n2].carga)/(4*pi*epce*epsi*d*esc);
    // += 8*ecoulomb(n1,n2,1);
    /*for(i=1;i<=8;i++){
        ei += ecoulomb(n1i,n1i,i+1) + ecoulomb(n2i,n2i,i+1);
        ef += ecoulomb(n1,n1,i+1) + ecoulomb(n2,n2,i+1);
    }*/
    dem = ef - ei;
    //printf("\n plasma eic: %e efc: %e de: %e",ei,ef,dem);
    //getchar();
    return(dem);
}
////////////////////////////////////////////////////////////////////////////////////
float ecoulomb(int a, int b, int c){
    int i, j;
    float d, ec=0;
    for(i=c;i<=8;i++){
        d = dist_ima(a,b,i);
        if(d==0){
            printf("d = 0, a: %i b: %i c: %i",a,b,c);
            getchar();
        }
        ec += (qe*qe*matriz_plasma[a].carga*matriz_plasma[b].carga)/(4*pi*epce*epsi*d*esc);
    }
    return(ec);
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
        if(a==1){
            c_unoa++;
        }
        else if(a==2){
            c_dosa++;
        }
        else{
            c_tresa++;
        }
    }
    else{
        matriz_plasma[n1] = matriz_plasma[n1i];
        matriz_plasma[n2] = matriz_plasma[n2i];
        rechazo_met_mov++;
    }
}
////////////////////////////////////////////////////////////////////////////////////
void dinamica(void){
    int i,j;
    int celdaobj,nmov;
    float vxx, vyy, dd, acel[2]={0,0}, B = 0.1*mhp/qe;
    for(i=1;i<=nceldas;i++){
        matriz_plasma2[i].electrones = 0;
        matriz_plasma2[i].hp = 0;
        matriz_plasma2[i].vxe = 0;
        matriz_plasma2[i].vye = 0;
        matriz_plasma2[i].vxhp = 0;
        matriz_plasma2[i].vyhp = 0;
    }
    for(i=1;i<=nceldas;i++){//mover particulas debido a la dinamica(8
        if(matriz_plasma[i].electrones!=0){
            //celdaobj = celda_dir(i);
            celdaobj = celda_dir_rec(i);
            //celdaobj = celda_dir_rec_neg(i);
            dd = distancianormal(i,celdaobj);
            ndt[i] = dd/sqrt(matriz_plasma[i].vxe*matriz_plasma[i].vxe+matriz_plasma[i].vye*matriz_plasma[i].vye);
            //printf("\nELECTRON\tdistancia: %f\tdt: %f",dd,ndt[i]);
            if(dd==0){
                printf("\n\nDD = 0\n");
                celdaobj = celda_dir_rec_error(i);
                getchar();
            }
            /*acel[0] = -qe*matriz_plasma[i].vye*B/me;
            acel[1] = qe*matriz_plasma[i].vxe*B/me;*/
            vxx = matriz_plasma[i].vxe;
            vyy = matriz_plasma[i].vye;
            //matriz_plasma[i].vxe = vxx + ndt[i]*acel[0];
            //matriz_plasma[i].vye = vyy + ndt[i]*acel[1];
            matriz_plasma[i].vxe = vxx*cos(ndt[i]*qe*B/me) + vyy*sin(-ndt[i]*qe*B/me);
            matriz_plasma[i].vye = -vxx*sin(-ndt[i]*qe*B/me) + vyy*cos(ndt[i]*qe*B/me) ;
            matriz_plasma2[celdaobj].vxe += matriz_plasma[i].electrones*matriz_plasma[i].vxe;
            matriz_plasma2[celdaobj].vye += matriz_plasma[i].electrones*matriz_plasma[i].vye;
            matriz_plasma2[celdaobj].electrones += matriz_plasma[i].electrones;
            /*if(matriz_plasma[i].vye!=1){
                printf("\nVALORES RAROS! vx: %f vy: %f\n",matriz_plasma[i].vxe,matriz_plasma[i].vye);
                getchar();
            }*/
            //getchar();
            matriz_plasma[i].electrones = 0;
            matriz_plasma[i].vxe = 0;
            matriz_plasma[i].vye = 0;
            /*if(i==nmatriz_plasma[ 62 ][ 62 ]){
                printf("\nceldaobj: %i\n",celdaobj);
            }*/
            /*celdaobj = int(alea()*nceldas)+1;
            if(celdaobj==nceldas+1)celdaobj=nceldas;
            nmov = int( alea_f(1, matriz_plasma[i].electrones ) );
            matriz_plasma[i].electrones-= nmov;
            matriz_plasma[celdaobj].electrones+= nmov;
            matriz_plasma[i].t[0]=distancianormal(i,celdaobj)/matriz_plasma[i].ve;

            celdaobj = int(alea()*nceldas)+1;
            if(celdaobj==nceldas+1)celdaobj=nceldas;
            nmov = int( alea_f(1, matriz_plasma[i].hp ) );
            matriz_plasma[i].hp-= nmov;
            matriz_plasma[celdaobj].hp+= nmov;
            matriz_plasma[i].t[1]=distancianormal(i,celdaobj)/matriz_plasma[i].vhp;

            celdaobj = int(alea()*nceldas)+1;
            if(celdaobj==nceldas+1)celdaobj=nceldas;
            nmov = int( alea_f(1, matriz_plasma[i].h2p ) );
            matriz_plasma[i].h2p-= nmov;
            matriz_plasma[celdaobj].h2p+= nmov;
            matriz_plasma[i].t[2]=distancianormal(i,celdaobj)/matriz_plasma[i].vh2p;*/
        }
        if(matriz_plasma[i].hp!=0){
            celdaobj = celda_dir_rec_hp(i);
            dd = distancianormal(i,celdaobj);
            ndt[i] = dd/sqrt(matriz_plasma[i].vxhp*matriz_plasma[i].vxhp+matriz_plasma[i].vyhp*matriz_plasma[i].vyhp);
            //printf("\nPROTON\tdistancia: %f\tdt: %f",dd,ndt[i]);
            if(dd==0){
                printf("\n\nDD = 0\n");
                celdaobj = celda_dir_rec_error_hp(i);
                getchar();
            }
            vxx = matriz_plasma[i].vxhp;
            vyy = matriz_plasma[i].vyhp;
            matriz_plasma[i].vxhp = vxx*cos(ndt[i]*qe*B/mhp) + vyy*sin(ndt[i]*qe*B/mhp);
            matriz_plasma[i].vyhp = vyy*cos(ndt[i]*qe*B/mhp) - vxx*sin(ndt[i]*qe*B/mhp);
            matriz_plasma2[celdaobj].vxhp += matriz_plasma[i].hp*matriz_plasma[i].vxhp;
            matriz_plasma2[celdaobj].vyhp += matriz_plasma[i].hp*matriz_plasma[i].vyhp;
            matriz_plasma2[celdaobj].hp += matriz_plasma[i].hp;
            /*if(matriz_plasma[i].vye!=1){
                printf("\nVALORES RAROS! vx: %f vy: %f\n",matriz_plasma[i].vxe,matriz_plasma[i].vye);
                getchar();
            }*/
            //getchar();
            matriz_plasma[i].hp = 0;
            matriz_plasma[i].vxhp = 0;
            matriz_plasma[i].vyhp = 0;
        }
    }
    for(i=1;i<=nceldas;i++){
        if(matriz_plasma2[i].electrones!=0){
            /*if(i==nmatriz_plasma[ 62 ][ 62 ]||i==nmatriz_plasma[ 62 ][ 62 ]+1){
                printf("ANTES\n");imprimir_celda_plasma(i);
                getchar();
            }*/
            vxx = (matriz_plasma[i].vxe*matriz_plasma[i].electrones+matriz_plasma2[i].vxe)/(matriz_plasma[i].electrones + matriz_plasma2[i].electrones);
            vyy = (matriz_plasma[i].vye*matriz_plasma[i].electrones+matriz_plasma2[i].vye)/(matriz_plasma[i].electrones + matriz_plasma2[i].electrones);
            /*if(vyy!=1){
                printf("Que show inf");
                printf("mp_i_vx: %f\tmp2_i_vx: %f\tmp_i_e: %f\tmp2_i_e: %f\t",matriz_plasma[i].vxe,matriz_plasma2[i].vxe,matriz_plasma[i].electrones,matriz_plasma2[i].electrones);
                getchar();
            }*/
            matriz_plasma[i].vxe = vxx;
            matriz_plasma[i].vye = vyy;
            matriz_plasma[i].electrones = matriz_plasma2[i].electrones;
            /*if(i==nmatriz_plasma[ 62 ][ 62 ]||i==nmatriz_plasma[ 62 ][ 62 ]+1){
                printf("DESPUES\n");imprimir_celda_plasma(i);
                getchar();
            }*/
        }
        if(matriz_plasma2[i].hp!=0){
            vxx = (matriz_plasma[i].vxhp*matriz_plasma[i].hp+matriz_plasma2[i].vxhp)/(matriz_plasma[i].hp + matriz_plasma2[i].hp);
            vyy = (matriz_plasma[i].vyhp*matriz_plasma[i].hp+matriz_plasma2[i].vyhp)/(matriz_plasma[i].hp + matriz_plasma2[i].hp);
            matriz_plasma[i].vxhp = vxx;
            matriz_plasma[i].vyhp = vyy;
            matriz_plasma[i].hp = matriz_plasma2[i].hp;
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////
void actu_salida(void){
	int i, j, cargatotaaal=0;
	float xx, yy;
    sprintf(salidac,"datos/posiciones%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    for(i=1;i<=nceldas;i++){
        printposiciones(i);
    }
    fclose(dat);

    sprintf(salidac,"datos/electron%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = (int)(matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
            //yy = (int)(matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
            //fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
            //fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].electrones);
            fprintf(dat,"%f\t%f\t%I64d\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
        //}
    }
    fclose(dat);

    sprintf(salidac,"datos/electron_polares%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
        fprintf(dat,"%f\t%f\t%I64d\n", matriz_plasma[i].rho, matriz_plasma[i].phi, matriz_plasma[i].electrones);
    }
    fclose(dat);

    sprintf(salidac,"datos/h2+%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h2+\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi);
        yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi);
        fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].h2p);
    }
    fclose(dat);

    sprintf(salidac,"datos/h20%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h20\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi);
        yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi);
        fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].h20);
    }
    fclose(dat);

    sprintf(salidac,"datos/hp%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"hp\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi);
        yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi);
        fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].hp);
    }
    fclose(dat);

    sprintf(salidac,"datos/todas%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi);
        yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi);
        fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].h2p+matriz_plasma[i].hp+matriz_plasma[i].h20+matriz_plasma[i].electrones);
    }
    fclose(dat);

    sprintf(salidac,"datos/carga%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"carga\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].carga);
        cargatotaaal += matriz_plasma[i].carga;
    }
    fprintf(dat,"carga total: %i\n",cargatotaaal);
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void salida_prom(void){
	int i, j, cargatotaaal=0;
	float xx, yy;

    sprintf(salidac,"datos/promedios/gelectron.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            //xx = (int)(matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
            //yy = (int)(matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
            //fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, gelectron[i]/(p-terma));
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/h2+.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h2+\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, gh2p[i]/(p-terma));
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/h20.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h20\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, gh20[i]/(p-terma));
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/hp.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"hp\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, ghp[i]/(p-terma));
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/todas.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, (gelectron[i]+gh2p[i]+gh20[i]+ghp[i])/(p-terma));
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/carga.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"carga\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%i\n", xx, yy, (gh2p[i]+ghp[i]-gelectron[i])/(p-terma));
            cargatotaaal += matriz_plasma[i].carga;
        }
    }
    fprintf(dat,"carga total: %i\n",cargatotaaal);
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
float energia(void){
    int i,j,dummy=0;
    float energia_total=0,d;
    dat = fopen("datos/energia.dat","a");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=nceldas;j++){
            if((matriz_plasma[i].carga!=0)&&(matriz_plasma[j].carga!=0)){
                if(i!=j){
                    //d = distanciacelda(i,j);
                    d = distancianormal(i,j);
                    energia_total += (qe*qe*matriz_plasma[i].carga*matriz_plasma[j].carga)/(4*pi*epce*epsi*d*esc);
                }
            }
        }
    }
    fprintf(dat,"\n%i\t%e",ienergia+p/actu,energia_total/2.0);
    fclose(dat);
    return(energia_total/2.0);
}
////////////////////////////////////////////////////////////////////////////////////
void beta(void){
    int i,j,jm;
    long long int npi,nei,nii;
    float vol_i,beta_i;
    dat = fopen("datos/beta.dat","w");
    fprintf(dat,"#X\tY\n");
    for(i=1;i<=R;i++){
        npi = nei = nii = 0;
        if(i==1){
            jm = 1;
        }
        else{
            jm = (int)((pi*i)*0.25);
        }
        vol_i = (pi*(pow(i,2)-pow(i-1,2))*esc*esc*1e-6)/8;
        for(j=1;j<=jm;j++){
            npi += matriz_plasma[ nmatriz_plasma[i][j] ].electrones + matriz_plasma[ nmatriz_plasma[i][j] ].h2p +
            matriz_plasma[ nmatriz_plasma[i][j] ].hp + matriz_plasma[ nmatriz_plasma[i][j] ].h20;
            nei += matriz_plasma[ nmatriz_plasma[i][j] ].electrones;
            nii += matriz_plasma[ nmatriz_plasma[i][j] ].h2p + matriz_plasma[ nmatriz_plasma[i][j] ].hp;
        }
        beta_i = (2*muce*kb*8*(nei*tempee+nii*tempei))/(
        vol_i*pow(arreglo[186].absb[i]*1e-4,2));
        if(beta_i<0){
            printf("\nmuce: %e npi: %I64d kb: %e tempee + tempeei: %e vol_i: %e |B_i|: %e",muce,npi,kb,tempee+tempei,vol_i,pow(arreglo[186].absb[i]*1e-4,1));
            //getchar();
        }
        fprintf(dat,"%f\t%e\n",i/10.0,beta_i);
    }
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void leer_campo_magnetico(void){
    int i,j;
    float m,b;
    dat = fopen("campo_magnetico.txt","r");
    fscanf(dat,"      R             Z              Br            Bz            |B|           A           dBz/dr        dBr/dz        Field");
    fscanf(dat,"    (cm)          (cm)             (G)           (G)           (G)         (G-cm)        (G/cm)        (G/cm)        Index");

    for(i=1;i<=371;i++){
        for(j=1;j<=36;j++){
            //fscanf(dat,"   %f       %f      %f  %f  %f  %f  %f  %f  %f\n",&arreglo[i].r[j],&arreglo[i].z,&arreglo[i].br[j],&arreglo[i].bz[j],&arreglo[i].absb[j],&arreglo[i].A[j],&arreglo[i].dbz[j],&arreglo[i].dbr[j],&arreglo[i].fieldindex[j]);
            fscanf(dat,"   %f",&arreglo[i].r[j]);
            fscanf(dat,"       %f",&arreglo[i].z);
            fscanf(dat,"      %f",&arreglo[i].br[j]);
            fscanf(dat,"  %f",&arreglo[i].bz[j]);
            fscanf(dat,"  %f",&arreglo[i].absb[j]);
            fscanf(dat,"  %f",&arreglo[i].A[j]);
            fscanf(dat,"  %f",&arreglo[i].dbz[j]);
            fscanf(dat,"  %f",&arreglo[i].dbr[j]);
            fscanf(dat,"  %f\n",&arreglo[i].fieldindex[j]);
            //printf("   %f       %f      %f  %f  %f  %f  %f  %f  %f\n",arreglo[i].r[j],arreglo[i].z,arreglo[i].br[j],arreglo[i].bz[j],arreglo[i].absb[j],arreglo[i].A[j],arreglo[i].dbz[j],arreglo[i].dbr[j],arreglo[i].fieldindex[j]);
        }
    }
    fclose(dat);
    m = (arreglo[186].absb[36]-arreglo[186].absb[1])/(arreglo[186].r[36]-arreglo[186].r[1]);
    b = arreglo[186].absb[36] - m*(arreglo[186].r[36]);
    for(i=1;i<=24;i++){
        arreglo[186].r[36+i] = 3.6+i/10.0;
        arreglo[186].absb[36+i] = m*arreglo[186].r[36+i]+b;
    }
    m = (arreglo[186].bz[36]-arreglo[186].bz[1])/(arreglo[186].r[36]-arreglo[186].r[1]);
    b = arreglo[186].bz[36] - m*(arreglo[186].r[36]);
    for(i=1;i<=24;i++){
        arreglo[186].r[36+i] = 3.6+i/10.0;
        arreglo[186].bz[36+i] = m*arreglo[186].r[36+i]+b;
    }
    m = (arreglo[186].br[36]-arreglo[186].br[1])/(arreglo[186].r[36]-arreglo[186].r[1]);
    b = arreglo[186].br[36] - m*(arreglo[186].r[36]);
    for(i=1;i<=24;i++){
        arreglo[186].r[36+i] = 3.6+i/10.0;
        arreglo[186].br[36+i] = m*arreglo[186].r[36+i]+b;
    }
    //Asignación de valores del campo magnetico
    /*for(i=1;i<=60;i++){
        campo_magnetico[i].absb = arreglo[186].absb[i];
        campo_magnetico[i].br = arreglo[186].br[i];
        campo_magnetico[i].bz = arreglo[186].bz[i];
    }*/
    //Asignación de valores del campo magnetico constante en z
    for(i=1;i<=60;i++){
        campo_magnetico[i].br = 0.0;
        campo_magnetico[i].bz = arreglo[186].bz[36];
        campo_magnetico[i].absb = fabs(campo_magnetico[i].bz);
    }
    printf("\n\n\n\nCAMPO MAGNETICO\nBZ1: %e", campo_magnetico[1].bz);
    printf("\nBZ2: %e", campo_magnetico[2].bz);
    getchar();
    /*printf("\nm: %f b: %f y(3.5): %f",m, b, m*arreglo[186].r[36]+b);
    for(i=1;i<=24;i++){
        printf("r: %f |B|: %f\n",arreglo[186].r[36+i], arreglo[186].absb[36+i]);
    }
    getchar();*/
}
////////////////////////////////////////////////////////////////////////////////////
void printposiciones(int a){
//fprintf(dat,"\n%4.0i\t%9.0I64d\t%9.0I64d\t%9.0I64d\t%9.0I64d\t%9.0i",i,matriz_plasma[i].electrones,matriz_plasma[i].h20,matriz_plasma[i].hp,matriz_plasma[i].h2p,matriz_plasma[i].carga);
    if(matriz_plasma[a].electrones==0){
        if(matriz_plasma[a].h20==0 ){
            if(matriz_plasma[a].hp==0){
                if(matriz_plasma[a].h2p==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%s\t%s\t%s\t%s",a,cero,cero,cero,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%s\t%s\t%s\t%10.0i",a,cero,cero,cero,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%s\t%s\t%10.0I64d\t%s",a,cero,cero,cero,matriz_plasma[a].h2p,cero);

                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%s\t%s\t%10.0I64d\t%10.0i",a,cero,cero,cero,matriz_plasma[a].h2p,matriz_plasma[a].carga);
                    }
                }
            }
            else{
                if(matriz_plasma[a].h2p==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%s\t%10.0I64d\t%s\t%s",a,cero,cero,matriz_plasma[a].hp,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%s\t%10.0I64d\t%s\t%10.0i",a,cero,cero,matriz_plasma[a].hp,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%s\t%10.0I64d\t%10.0I64d\t%s",a,cero,cero,matriz_plasma[a].hp,matriz_plasma[a].h2p,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%s\t%10.0I64d\t%10.0I64d\t%10.0i",a,cero,cero,matriz_plasma[a].hp,matriz_plasma[a].h2p,matriz_plasma[a].carga);
                    }
                }
            }
        }
        else{
            if(matriz_plasma[a].hp==0){
                if(matriz_plasma[a].h2p==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%10.0I64d\t%s\t%s\t%s",a,cero,matriz_plasma[a].h20,cero,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%10.0I64d\t%s\t%s\t%10.0i",a,cero,matriz_plasma[a].h20,cero,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%10.0I64d\t%s\t%10.0I64d\t%s",a,cero,matriz_plasma[a].h20,cero,matriz_plasma[a].h2p,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%10.0I64d\t%s\t%10.0I64d\t%10.0i",a,cero,matriz_plasma[a].h20,cero,matriz_plasma[a].h2p,matriz_plasma[a].carga);
                    }
                }
            }
            else{
                if(matriz_plasma[a].h2p==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%10.0I64d\t%10.0I64d\t%s\t%s",a,cero,matriz_plasma[a].h20,matriz_plasma[a].hp,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%10.0I64d\t%10.0I64d\t%s\t%10.0i",a,cero,matriz_plasma[a].h20,matriz_plasma[a].hp,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%10.0I64d\t%10.0I64d\t%10.0I64d\t%s",a,cero,matriz_plasma[a].h20,matriz_plasma[a].hp,matriz_plasma[a].h2p,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%10.0I64d\t%10.0I64d\t%10.0I64d\t%10.0i",a,cero,matriz_plasma[a].h20,matriz_plasma[a].hp,matriz_plasma[a].h2p,matriz_plasma[a].carga);
                    }
                }
            }
        }
    }
    else{
        if(matriz_plasma[a].h20==0 ){
            if(matriz_plasma[a].hp==0){
                if(matriz_plasma[a].h2p==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0I64d\t%s\t%s\t%s\t%s",a,matriz_plasma[a].electrones,cero,cero,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0I64d\t%s\t%s\t%s\t%10.0i",a,matriz_plasma[a].electrones,cero,cero,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0I64d\t%s\t%s\t%10.0I64d\t%s",a,matriz_plasma[a].electrones,cero,cero,matriz_plasma[a].h2p,cero);

                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0I64d\t%s\t%s\t%10.0I64d\t%10.0i",a,matriz_plasma[a].electrones,cero,cero,matriz_plasma[a].h2p,matriz_plasma[a].carga);
                    }
                }
            }
            else{
                if(matriz_plasma[a].h2p==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0I64d\t%s\t%10.0I64d\t%s\t%s",a,matriz_plasma[a].electrones,cero,matriz_plasma[a].hp,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0I64d\t%s\t%10.0I64d\t%s\t%10.0i",a,matriz_plasma[a].electrones,cero,matriz_plasma[a].hp,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0I64d\t%s\t%10.0I64d\t%10.0I64d\t%s",a,matriz_plasma[a].electrones,cero,matriz_plasma[a].hp,matriz_plasma[a].h2p,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0I64d\t%s\t%10.0I64d\t%10.0I64d\t%10.0i",a,matriz_plasma[a].electrones,cero,matriz_plasma[a].hp,matriz_plasma[a].h2p,matriz_plasma[a].carga);
                    }
                }
            }
        }
        else{
            if(matriz_plasma[a].hp==0){
                if(matriz_plasma[a].h2p==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0I64d\t%10.0I64d\t%s\t%s\t%s",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,cero,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0I64d\t%10.0I64d\t%s\t%s\t%10.0i",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,cero,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0I64d\t%10.0I64d\t%s\t%10.0I64d\t%s",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,cero,matriz_plasma[a].h2p,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0I64d\t%10.0I64d\t%s\t%10.0I64d\t%10.0i",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,cero,matriz_plasma[a].h2p,matriz_plasma[a].carga);
                    }
                }
            }
            else{
                if(matriz_plasma[a].h2p==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0I64d\t%10.0I64d\t%10.0I64d\t%s\t%s",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,matriz_plasma[a].hp,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0I64d\t%10.0I64d\t%10.0I64d\t%s\t%10.0i",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,matriz_plasma[a].hp,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0I64d\t%10.0I64d\t%10.0I64d\t%10.0I64d\t%s",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,matriz_plasma[a].hp,matriz_plasma[a].h2p,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0I64d\t%10.0I64d\t%10.0I64d\t%10.0I64d\t%10.0i",a,matriz_plasma[a].electrones,matriz_plasma[a].h20,matriz_plasma[a].hp,matriz_plasma[a].h2p,matriz_plasma[a].carga);
                    }
                }
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////
void gdr_plasma(void){
    int i,j;
    for(i=1;i<=nceldas;i++){
        gelectron[i] += matriz_plasma[i].electrones;
        gh20[i] += matriz_plasma[i].h20;
        ghp[i] += matriz_plasma[i].hp;
        gh2p[i] += matriz_plasma[i].h2p;
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir(int a){
    float B = campo_magnetico[int(matriz_plasma[a].rho*10)].bz, normav, normar, rpv, phipv;
    int int_rho, int_phi, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vxe*matriz_plasma[a].vxe + matriz_plasma[a].vye*matriz_plasma[a].vye );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vxe;  v1[1] = matriz_plasma[a].vye; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/

    rpv = (matriz_plasma[a].vxe*matriz_plasma[a].x)/(normar*normav) + (matriz_plasma[a].vye*matriz_plasma[a].y)/(normar*normav);
    phipv = -(matriz_plasma[a].vxe*matriz_plasma[a].y)/(normar*normav) + (matriz_plasma[a].vye*matriz_plasma[a].x)/(normar*normav);

    if(rpv<-(1.0/3.0)){
        int_rho = -1;
    }else if(rpv<1.0/3.0){
        int_rho = 0;
    }else{
        int_rho = 1;
    }
    if(phipv<-(1.0/3.0)){
        int_phi = -1;
    }else if(phipv<1.0/3.0){
        int_phi = 0;
    }else{
        int_phi = 1;
    }
    int_i = (int)( matriz_plasma[a].rho );
    int_j = (int)( (matriz_plasma[a].phi - 0.125*pi/((int)((pi*int_i)*0.25)))/(pi*0.25) )*(int)((pi*int_i)*0.25) + (int)( (matriz_plasma[a].phi - 0.125*pi/((int)((pi*int_i)*0.25)))/((int)((pi*int_i)*0.25)) )*(int)((pi*int_i)*0.25);
    //printf("int_i_0:%i \t int_j_0: %i\tint_rho: %i \t int_phi: %i",int_i, int_j, int_rho,int_phi);
    int_i+=int_rho;
    int_j+=int_phi;
    /*printf("\nint_i_final:%i \t int_j_final: %i matriz[int_i][int_j]: %i",int_i, int_j, nmatriz_plasma[ int_i ][ int_j ]);
    getchar();*/
    return( nmatriz_plasma[ int_i ][ int_j ] );
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec(int a){
    float B = campo_magnetico[int(matriz_plasma[a].rho*10)].bz, normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vxe*matriz_plasma[a].vxe + matriz_plasma[a].vye*matriz_plasma[a].vye );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    /*if(matriz_plasma[a].electrones!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vxe,matriz_plasma[a].vye);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vxe;  v1[1] = matriz_plasma[a].vye; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    //pphi = (atan2(matriz_plasma[a].vye,matriz_plasma[a].vxe)+((atan2(matriz_plasma[a].vye,matriz_plasma[a].vxe))<0?2*pi:0.0));
    //printf("pphi: %f",pphi);
    //getchar();
    if(matriz_plasma[a].vxe/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vxe/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vye/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vye/normav<1.0/3.0){
        int_y = 0;
    }else{
        int_y = 1;
    }
    int_i = floor( matriz_plasma[a].x*reso) + floor(R*reso) + 1 + int_x;
    int_j = floor( matriz_plasma[a].y*reso) + floor(R*reso) + 1 + int_y;
    if(a==10000000000000000){
        printf("\nANTES DE CELDA_ERROR\nCELDA INICIAL: %i i: %i j: %i", a, (int)( matriz_plasma[a].x*reso + R*reso + 1 ), (int)( matriz_plasma[a].y*reso + R*reso + 1 ));
        printf("\nCELDA OBJETIVO: %i i: %i j: %i", nmatriz_plasma[ int_i ][ int_j ], int_i, int_j );
        printf("\nint_x: %i int_y: %i",int_x, int_y);
        printf("\ny*reso: %f int(y*reso): %i\n",matriz_plasma[a].y*reso, (((int) (matriz_plasma[a].y*reso + 32768.5)) - 32768));
        printf("\nR*reso: %f int(R*reso): %i\n",R*reso, int(R*reso));
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vxe,matriz_plasma[a].vye);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vxe,matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vye);
        imprimir_celda_plasma_rec(a);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        getchar();
    }
    /*if(nmatriz_plasma[ int_i ][ int_j ]>=12345||nmatriz_plasma[ int_i ][ int_j ]==0){
        printf("int_i: %i\tint_j: %i",int_i,int_j);
        return( nmatriz_plasma[ int_i ][ int_j ] );*/
    //rrho = sqrt(matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x+matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y);
    if(matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho>60){
        if(a==221099000){
            printf("Que show r: %f",sqrt(int_i*int_i+int_j*int_j));
            getchar();
        }
        return( nmatriz_plasma[ (int)( matriz_plasma[a].x*reso + R*reso + 1 ) ][ (int)( matriz_plasma[a].y*reso + R*reso + 1 ) ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec_error(int a){
    float B = campo_magnetico[int(matriz_plasma[a].rho*10)].bz, normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vxe*matriz_plasma[a].vxe + matriz_plasma[a].vye*matriz_plasma[a].vye );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    /*if(matriz_plasma[a].electrones!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vxe,matriz_plasma[a].vye);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vxe;  v1[1] = matriz_plasma[a].vye; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    //pphi = (atan2(matriz_plasma[a].vye,matriz_plasma[a].vxe)+((atan2(matriz_plasma[a].vye,matriz_plasma[a].vxe))<0?2*pi:0.0));
    //printf("pphi: %f",pphi);
    //getchar();
    if(matriz_plasma[a].vxe/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vxe/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vye/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vye/normav<1.0/3.0){
        int_y = 0;
    }else{
        int_y = 1;
    }
    int_i = floor( matriz_plasma[a].x*reso) + floor(R*reso) + 1 + int_x;
    int_j = floor( matriz_plasma[a].y*reso) + floor(R*reso) + 1 + int_y;
    if(a>=1e30){
        printf("\nANTES DE CELDA_ERROR\nCELDA INICIAL: %i i: %i j: %i", a, (int)( matriz_plasma[a].x*reso + R*reso + 1 ), (int)( matriz_plasma[a].y*reso + R*reso + 1 ));
        printf("\nCELDA OBJETIVO: %i i: %i j: %i", nmatriz_plasma[ int_i ][ int_j ], int_i, int_j );
        printf("\nint_x: %i int_y: %i",int_x, int_y);
        printf("\nx*reso: %f int(x*reso): %i\n",matriz_plasma[a].x*reso, floor(matriz_plasma[a].x*reso));
        printf("\ny*reso: %f int(y*reso): %i\n",matriz_plasma[a].y*reso, floor(matriz_plasma[a].y*reso));
        printf("\nR*reso: %f int(R*reso): %i\n",R*reso, floor(R*reso));
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vxe,matriz_plasma[a].vye);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vxe,matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vye);
        imprimir_celda_plasma_rec(a);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        getchar();
    }
    /*if(nmatriz_plasma[ int_i ][ int_j ]>=12345||nmatriz_plasma[ int_i ][ int_j ]==0){
        printf("int_i: %i\tint_j: %i",int_i,int_j);
        return( nmatriz_plasma[ int_i ][ int_j ] );*/
    //rrho = sqrt(matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x+matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y);
    if(matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho>R){
        if(a==221099000){
            printf("Que show r: %f",sqrt(int_i*int_i+int_j*int_j));
            getchar();
        }
        return( nmatriz_plasma[ (int)( matriz_plasma[a].x*reso + R*reso + 1 ) ][ (int)( matriz_plasma[a].y*reso + R*reso + 1 ) ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec_hp(int a){
    float B = campo_magnetico[int(matriz_plasma[a].rho*10)].bz, normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vxhp*matriz_plasma[a].vxhp + matriz_plasma[a].vyhp*matriz_plasma[a].vyhp );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    /*if(matriz_plasma[a].electrones!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vxe,matriz_plasma[a].vye);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vxe;  v1[1] = matriz_plasma[a].vye; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    //pphi = (atan2(matriz_plasma[a].vye,matriz_plasma[a].vxe)+((atan2(matriz_plasma[a].vye,matriz_plasma[a].vxe))<0?2*pi:0.0));
    //printf("pphi: %f",pphi);
    //getchar();
    if(matriz_plasma[a].vxhp/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vxhp/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vyhp/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vyhp/normav<1.0/3.0){
        int_y = 0;
    }else{
        int_y = 1;
    }
    int_i = floor( matriz_plasma[a].x*reso) + floor(R*reso) + 1 + int_x;
    int_j = floor( matriz_plasma[a].y*reso) + floor(R*reso) + 1 + int_y;
    if(matriz_plasma[a].electrones>=1e30){
        printf("celda inicial: %i celda objetivo: %i int_i: %i int_j:%i\n",a,nmatriz_plasma[ int_i ][ int_j ],int_i,int_j);
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vxe,matriz_plasma[a].vye);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vxe,matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vye);
        imprimir_celda_plasma_rec(a);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        getchar();
    }
    //rrho = sqrt(matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x+matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y);
    if(matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho>R){
            printf("Que show i: %i r: %f",nmatriz_plasma[int_i][int_j],matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho);
            getchar();
        return( nmatriz_plasma[ int_i-int_x ][ int_j-int_y ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec_error_hp(int a){
    float B = campo_magnetico[int(matriz_plasma[a].rho*10)].bz, normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vxhp*matriz_plasma[a].vxhp + matriz_plasma[a].vyhp*matriz_plasma[a].vyhp );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    /*if(matriz_plasma[a].electrones!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vxe,matriz_plasma[a].vye);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vxe;  v1[1] = matriz_plasma[a].vye; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    //pphi = (atan2(matriz_plasma[a].vye,matriz_plasma[a].vxe)+((atan2(matriz_plasma[a].vye,matriz_plasma[a].vxe))<0?2*pi:0.0));
    //printf("pphi: %f",pphi);
    //getchar();
    if(matriz_plasma[a].vxhp/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vxhp/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vyhp/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vyhp/normav<1.0/3.0){
        int_y = 0;
    }else{
        int_y = 1;
    }
    int_i = floor( matriz_plasma[a].x*reso) + floor(R*reso) + 1 + int_x;
    int_j = floor( matriz_plasma[a].y*reso) + floor(R*reso) + 1 + int_y;
    //if(a>=1e30){
        printf("\nANTES DE CELDA_ERROR\nCELDA INICIAL: %i i: %f j: %f", a, floor( matriz_plasma[a].x*reso) + floor(R*reso) + 1, floor( matriz_plasma[a].y*reso) + floor(R*reso) + 1);
        printf("\nCELDA OBJETIVO: %i i: %i j: %i", nmatriz_plasma[ int_i ][ int_j ], int_i, int_j );
        printf("\nint_x: %i int_y: %i",int_x, int_y);
        printf("\nx*reso: %f int(x*reso): %f\n",matriz_plasma[a].x*reso, floor(matriz_plasma[a].x*reso));
        printf("\ny*reso: %f int(y*reso): %f\n",matriz_plasma[a].y*reso, floor(matriz_plasma[a].y*reso));
        printf("\nR*reso: %f int(R*reso): %f\n",R*reso, floor(R*reso));
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vxhp,matriz_plasma[a].vyhp);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vxhp,matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vyhp);
        imprimir_celda_plasma_rec(a);
        printf("\nrho: %f",matriz_plasma[a].rho);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        printf("\nrho: %f",matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].rho);
        getchar();
    //}
    /*if(nmatriz_plasma[ int_i ][ int_j ]>=12345||nmatriz_plasma[ int_i ][ int_j ]==0){
        printf("int_i: %i\tint_j: %i",int_i,int_j);
        return( nmatriz_plasma[ int_i ][ int_j ] );*/
    //rrho = sqrt(matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x+matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y);
    if(matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho>R){
        if(a==221099000){
            printf("Que show r: %f",sqrt(int_i*int_i+int_j*int_j));
            getchar();
        }
        return( nmatriz_plasma[ (int)( matriz_plasma[a].x*reso + R*reso + 1 ) ][ (int)( matriz_plasma[a].y*reso + R*reso + 1 ) ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec_neg(int a){
    float normav, normar, rrho, xx, yy, dt = 2e0, B = me/qe;//campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    float acelneg[2];
    int int_i, int_j, nceldaobj;
    /*if(matriz_plasma[a].electrones!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vxe,matriz_plasma[a].vye);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vxe;  v1[1] = matriz_plasma[a].vye; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    acelneg[0] = -qe*matriz_plasma[a].vye*B/me;
    acelneg[1] = qe*matriz_plasma[a].vxe*B/me;
    xx = matriz_plasma[a].x + matriz_plasma[a].vxe*dt + 0.5*acelneg[0]*dt*dt;
    yy = matriz_plasma[a].y + matriz_plasma[a].vye*dt + 0.5*acelneg[1]*dt*dt;
    printf("x_a_0: %f y_a_0: %f\tx_a_t: %f y_a_t: %f",matriz_plasma[a].x,matriz_plasma[a].y,xx,yy);
    printf("\ndato1: %f dato2: %f dato3: %f dato4: %f", matriz_plasma[a].x,matriz_plasma[a].vxe,dt,acelneg[0] );
    printf("\ndato1: %f dato2: %f dato3: %f dato4: %f", matriz_plasma[a].y,matriz_plasma[a].vye,dt,acelneg[1] );
    getchar();
    int_i = (int)( xx + 61 + 1 );
    int_j = (int)( yy + 61 + 1 );
    nceldaobj = nmatriz_plasma[ int_i ][ int_j ];
    /*if(matriz_plasma[a].electrones>=1e30){
        printf("celda inicial: %i celda objetivo: %i int_i: %i int_j:%i\n",a,nmatriz_plasma[ int_i ][ int_j ],int_i,int_j);
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vxe,matriz_plasma[a].vye);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vxe,matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vye);
        imprimir_celda_plasma_rec(a);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        getchar();
    }*/
    rrho = sqrt(matriz_plasma[ nceldaobj ].x*matriz_plasma[ nceldaobj ].x+matriz_plasma[ nceldaobj ].y*matriz_plasma[ nceldaobj ].y);
    if(rrho>60){
        printf("Que show r: %f",sqrt(int_i*int_i+int_j*int_j));
        getchar();
        return( nmatriz_plasma[ (int)( matriz_plasma[a].x + 61 + 1 ) ][ (int)( matriz_plasma[a].y + 61 + 1 ) ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
vector<float> vect(float a[3]){
    vector<float> resultado;
    resultado.push_back(a[0]);
    resultado.push_back(a[1]);
    resultado.push_back(a[2]);
    return(resultado);
}
////////////////////////////////////////////////////////////////////////////////////
void print_vector(vector<float> a){
    int i;
    for(i=0;i<a.size();i++){
        cout << a[i] << "\t";
    }
    cout << endl;
}
////////////////////////////////////////////////////////////////////////////////////
float punto(vector<float> a, vector<float> b){
    return(a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}
////////////////////////////////////////////////////////////////////////////////////
vector<float> cruz(vector<float> a, vector<float> b){
    vector<float> resultado;
    resultado.push_back(a[1]*b[2]-b[1]*a[2]);
    resultado.push_back(a[2]*b[0]-b[2]*a[0]);
    resultado.push_back(a[0]*b[1]-b[0]*a[1]);
    return(resultado);
}
////////////////////////////////////////////////////////////////////////////////////
void prueba_dinamica(void){
    float x,y,vx,vy,ax,ay;
    float Ex=1.0,Ey=2.0,Bz=1.0,dt=0.000001;
    float x0=0.0, v0x=1.0, y0 = 0.0, v0y = 0.0;
    int i;

    dat=fopen("prueba_dinamica.dat","w");
    fprintf(dat,"#t\tx\tvx\ty\tvy\n");
        fprintf(dat,"%f\t%f\t%f\t%f\t%f\n",0.0,x0,v0x,y0,v0y);
    for(i=1;i<=100000000;i++){
        ax = v0y*Bz;
        ay = -v0x*Bz;
        x = x0 + v0x*dt + 0.5*ax*dt*dt;
        vx = v0x + ax*dt;
        y = y0 + v0y*dt + 0.5*ay*dt*dt;
        vy = v0y + ay*dt;
        if(i%10000==0)fprintf(dat,"%f\t%f\t%f\t%f\t%f\n",i*dt,x,vx,y,vy);
        x0 = x; v0x = vx;
        y0 = y; v0y = vy;
    }
}
////////////////////////////////////////////////////////////////////////////////////
float aaadist_normal(float mu, float sigma){//usando 5 sigmas de intervalo
    float b, fb;
    do{
        b=alea_f(-5*sigma+mu,5*sigma+mu);//sqrt((kb*tempee)/me), 4*sqrt((kb*tempee)/me) );
        fb=exp(-(b-mu)*(b-mu)/(2*sigma*sigma))/(sqrt(2*pi)*sigma);
    }while(alea()/(sqrt(2*pi)*sigma)>fb);
    return(b);
}
////////////////////////////////////////////////////////////////////////////////////
void aaaprueba_dist(void){
    int i,j;
    long int clases[201] = {0};
    float x, mu = 0, sigma = 0.005;

    /*for(i=1;i<=1000000;i++){
        clases[ (int)(200*(aaadist_normal(sigma)+4*sigma)/(8*sigma))+1 ]++;
    }*/
    for(i=1;i<=10000000;i++){
        clases[ int( 200*(aaadist_normal(mu,sigma)+5*sigma-mu)/(10*sigma) )+1]++;
    }

    dat = fopen("aaahistograma.dat","w");
    for(i=1;i<=200;i++){
        fprintf(dat,"%f\t%d\n",-5*sigma+(i-1)*10*sigma/200.0+mu,clases[i]);
    }
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void aaprueba(double &a){
    a = 1.0;
}
	
