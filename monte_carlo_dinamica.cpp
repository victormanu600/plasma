#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
//#include <tgmath.h>
//#include <windows.h>
#include <vector>
#include <fstream>

using namespace std;

#define pi 3.141592653
#define epce 8.854187817e-12
#define epsi 1.0//78.5
#define muce 1.256637061e-6
#define qe 1.602176634e-19
#define kb 1.3806488e-23
#define na 6.022140857e23
#define m_u 1.66053906660e-27

float alea(void);
float alea_f(float a, float b);
int alea_i(int a, int b);
long long int max_alea(void);
void leer_datos_iniciales(void);
void leer_campo_magnetico(void);
void imprimir_datos_iniciales(void);
void imprimir_celda_plasma_rec(int a);
void imprimir_celda_plasma_part(int a, int b);
void condiciones_iniciales(void);
void arreglo_inicial(void);
void importar_part_a_mallado(void);
int part_a_mallado(float a, float b);
void intercambiar(int a, int b, int c );
void crear_matriz_plasma(void);
void vecinos(int a);
void calc_carga(int a);
void hacer_histograma(int a, int b);
void hacer_distribucion(int a);
void hacer_rms(int a);
void crear_archivos_iniciales_en_blanco(void);

float distancia(int a, int b);
int signo(float a);
float norma(float a, float b);
void mover_particulas_part(void);
void pared(void);
void reflejar_pared(int a, int b, float* vx, float* vy);
void metropolis_plasma(void);
void dinamica(void);
void cambiar_a(int* a, float* b);

float de_plasma(void);
float calcular_de_mov(void);
float autoenergia(int a);
void calcular_momento_total(void);
void gdr_plasma(void);
void actu_salida(void);
void actu_salida_carga(void);
void salida_prom(void);
float energia(void);
void beta(void);

void prueba_dinamica(void);
float distribucion_normal_5(float mu, float sigma);
void prueba_dist(void);

void mainsillo(void);


//////////////////////////////Constantes
int pasos, actu, terma, pasoinicial=1, termaanterior, iprint[6], ncomp;
bool solodinamica;
int contadorrr = 0;
float esc, densidad, volumen, nelectronesr, tempee, tempei, tempe[4], B, reso;
long long int nh20, nhp, nh2p, nelectrones, npart[4], naleprom=0, npart3[4] ;
float R;
const int R_reso = 2*(12)*3/1+1, part_plasma_size = 6000001;                      //Primer numero es R, segundo es reso, 2 es por que es de -R a R y +1 para comenzar los arreglos desde 1
float me, mhp, mh2p, mh20, masa[4];
int carga[4];
int n_grupos[4];
///////////////////////////////////////////////////////////////////////Variables globales
int p=0, tipo, ienergia;
int rechazo;
int nceldas, n1, n2, n1i, n2i, alea_part;                                       //numero de celda Nueva/Vieja Inicial/Final
char cero[] = "         0";
int rasdnd;
double dem, de_coulomb = 0, de_autoenergia = 0;
double dt_termico = 0;
long long int part_desplazadas[4] = {0};
int debug=0;
float k_auto = 0, dtau = 0;
bool win_os;
///////////////////////////////////////////////////////////////////////Contadores
int c_mov=0,c_mova=1, c_dina=0;
int c_uno=0,c_unoa=0,c_dos=0,c_dosa=0,c_tres=0,c_tresa=0;
int contador[4],contador_a[4], contador_rechazo=0, contador_nale=0,contador_nalef=0;
int rechazo_neg=0, rechazo_met_mov=0;

long long int gelectron[11260], gh20[11260], gh2p[11260], ghp[11260], gpart[11260][4];

struct smatriz_plasma{
	float rho,phi,x,y;
    float vxe, vye, vxh2p, vyh2p, vxhp, vyhp, vxh20, vyh20;
    float vx[4], vy[4];//, vx_1[4], vy_1[4], vx_2[4], vy_2[4];
    float ve, vhp, vh2p, v[4];
    float anchow, anchol;
	float t[3];
	long long int h20,hp,h2p,electrones,part[4],carga;
	int xi, yi;
	int nvecinos, vecinos[10];
//}matriz_plasma[8*61*50+1];
//}matriz_plasma[7*61*60+1];
//}matriz_plasma[12346+10000];
}matriz_plasma[R_reso*R_reso+1+2];
//}matriz_plasma[3000*3000+1];

struct smatriz_plasma2{
	float rho,phi,x,y;
    float vxe, vye, vxh2p, vyh2p, vxhp, vyhp, vxh20, vyh20;
    float vx[4], vy[4];//, vx_1[4], vy_1[4], vx_2[4], vy_2[4];
    float ve, vhp, vh2p, v[4];
    float anchow, anchol;
	float t[3];
	long long int h20,hp,h2p,electrones,part[4],carga;
	int xi, yi;
	int nvecinos, vecinos[10];
}matriz_plasma2[R_reso*R_reso+1+2];
//}matriz_plasma2[2000*2000+1];

struct spart_plasma{
    float x, y, vx, vy, v;
    long long int part;
    int contador;
//}part_plasma[1000*(R_reso*R_reso+1)][4];
}part_plasma[part_plasma_size][4];

float ndt[R_reso*R_reso+1+2]={0};

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

/*struct _termostato{
    float vx[4][R_reso], vy[4][R_reso];
}termostato[R_reso+1];*/

//int nmatriz_plasma[60+1][7*60];
int nmatriz_plasma[R_reso+1][R_reso+1];
long long int nale;
float tale;
//int nmatriz_plasma[2000+1][2000+1]={0};
////////////////////////////////////////////////////////////////////////Variables para salida
FILE *dat, *dat2;
char salidac[10];

int main(){
    //prueba_dinamica();
    srand((unsigned)time(NULL));

    #ifdef _WIN32 // Includes both 32 bit and 64 bit
	    #ifdef _WIN64
            win_os = true;
	        printf("Windows 64 bit\n");
	    #else
            win_os = true;
	        printf("Windows 32 bit\n");
	    #endif
	#else
        win_os = false;
	    printf("Not a Windows OS\n");
	#endif

    crear_archivos_iniciales_en_blanco();

	if(debug==1)cout << endl << "Condiciones iniciales" << endl;
	condiciones_iniciales();
	if(debug==1)cout << endl << "Crear matriz" << endl;
	crear_matriz_plasma();
	if(debug==1)cout << endl << "Arreglo inicial" << endl;
    arreglo_inicial(); cout << "Energia: " << energia() << endl;
	if(debug==1)cout << endl << "Despues arreglo inicial" << endl;
	//actu_salida_carga();
    //return(0);

    for(p=pasoinicial;p<=pasos;p++){
        if(debug==2)cout << "\nINICIO CICLO\n";
        //if((p%(100*nceldas)!=0||p<=terma)&&(p>=pasoinicial)&&(!solodinamica)){
        //if(alea()<2&&!solodinamica){
        if((dt_termico/(npart[0]+npart[1]))<dtau&&!solodinamica){
            rechazo = 0; c_mov++;

            if(debug==2)cout << "\nANTES MOVER PART\n";
            mover_particulas_part();
            part_plasma[n1][tipo].contador++;

            if(debug==2)cout << "\nANTES PARED\n";
            pared();

            if(debug==2)cout << "\nANTES METROPOLIS\n";
            if(rechazo == 0){
                metropolis_plasma();
            }else{
                rechazo_neg++;
                part_plasma[alea_part][tipo] = part_plasma[0][tipo];
                matriz_plasma[n1] = matriz_plasma[n1i];
                matriz_plasma[n2] = matriz_plasma[n2i];
                //importar_part_a_mallado();
            }
        }else{
            c_dina++;
            if(debug==2)cout << "\nANTES DINAMICA\n";
            dinamica();
            if(debug==2)cout << "\nANTES IMPORTAR PART A MALLADO\n";
            importar_part_a_mallado();
            dt_termico = 0;
        }
        if(p>terma){
            if(debug==1)cout << "\nANTES GDR2\n";
            gdr_plasma();
            if(p%actu==0){
                if(debug==1)cout << "\nANTES ACTUSALIDA2\n";
                actu_salida();
                if(debug==1)cout << "\nANTES SALIDAPROM\n";
                salida_prom();
                for(int i=0; i<ncomp; i++)hacer_histograma(i,p);
            }
        }
        if(p%1000==0||solodinamica){
            if(win_os){
                printf("\rPaso: %9d ",p);
                for(int i=0; i<ncomp; i++)printf("Acep%d: %1.5f ",i,1.0*contador_a[i]/contador[i]);
                printf("dt/npart: %e ",(dt_termico/(npart[0]+npart[1])));
                printf("Rneg: %d uno: %1.5f ",rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0));
                printf("Dinamica: %d ",c_dina);
                if(p%actu==0)printf("Energia: %e ",energia());
            }

            dat = fopen("pruebas/razon_acep.dat","a");
            fprintf(dat,"%d",p);
            for(int i=0; i<ncomp; i++)fprintf(dat,"\t%1.5f",contador_a[i]/(contador[i]*1.0));
            fprintf(dat,"\n");
            fclose(dat);
            //contador_nale=0;naleprom=0,nale_f=0;
        }

        long long int cargatotaaal = 0;
        for(int i=1;i<=nceldas;i++){
            cargatotaaal += matriz_plasma[i].carga;
        }
        if(cargatotaaal!=0){
            printf("\ncargatotal: %lld alea_part: %d tipo: %d celda: %d",cargatotaaal,alea_part,tipo,part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y));
            imprimir_celda_plasma_part(alea_part,tipo);
            printf("\nPaso: %d ",p);
            for(int i=0; i<ncomp; i++)printf("Acep%d: %1.5f ",i,1.0*contador_a[i]/contador[i]);
            printf("dt/npart: %e ",(dt_termico/(npart[0]+npart[1])));
            printf("Rneg: %d uno: %1.5f ",rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0));
            printf("Dinamica: %d",c_dina);
            getchar();
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////FIN DE MAIN
////////////////////////////////////////////////////////////////////////////////////
float alea(void){
    return((float)rand()/RAND_MAX);
}
////////////////////////////////////////////////////////////////////////////////////
float alea_f(float a, float b){
    return( ( a + ( (float)rand()/RAND_MAX )*( b - a ) ) );
}
////////////////////////////////////////////////////////////////////////////////////
int alea_i(int a, int b){//Solo para numeros enteros pequeños
    //return( ( a + int( ( (float)rand()/(RAND_MAX+1) )*( 1 + b - a ) ) ) );
    return(  win_os?( a + max_alea()%( b - a + 1 ) ):( a + rand()%( b - a + 1 ) )   );
    //return(  a + max_alea()%( b - a + 1 )   );
}
///////////////////////////////////////////////////////////////////////////////////////////////
long long int max_alea(void){
    return(rand() + rand()*(RAND_MAX+1));
}
////////////////////////////////////////////////////////////////////////////////////
void leer_datos_iniciales(void){
    float dummy;
    if(dat=fopen("entrada.txt","r")){
        fscanf(dat,"Numero de pasos: %i\n", &pasos);
        fscanf(dat,"Actualizacion: %i\n", &actu);
        fscanf(dat,"Termalizacion: %i\n", &terma);
        fscanf(dat,"Resolucion: %f\n", &reso);
        fscanf(dat,"Componentes: %d\n", &ncomp);
        fscanf(dat,"Masas (kg):");
        for(int i=0; i<ncomp; i++){
            if(i<ncomp-1){
                fscanf(dat," %f,",&masa[i]);
            }else{
                fscanf(dat," %f\n",&masa[i]);
            }
        }
        //masa[0]*=m_u;masa[1]*=m_u;masa[2]*=m_u;masa[3]*=m_u;
        fscanf(dat,"Carga (e):");
        for(int i=0; i<ncomp; i++){
            if(i<ncomp-1){
                fscanf(dat," %d,",&carga[i]);
            }else{
                fscanf(dat," %d\n",&carga[i]);
            }
        }
        fscanf(dat,"R (mm): %f\n",&R);
        fscanf(dat,"Temperaturas (eV):");
        for(int i=0; i<ncomp; i++){
            if(i<ncomp-1){
                fscanf(dat," %f,",&tempe[i]);
            }else{
                fscanf(dat," %f\n",&tempe[i]);
            }
        }
        fscanf(dat,"Campo magnetico (T): %f\n",&B);
        fscanf(dat,"Constante AE: %f\n",&k_auto);
        fscanf(dat,"Tiempo para dinamica (s): %f\n",&dtau);
        if(B<1e-10)B=1e-10;
        fscanf(dat,"Iprint :");
        for(int i=0; i<ncomp+2; i++){
            if(i<ncomp+1){
                fscanf(dat," %d,",&iprint[i]);
            }else{
                fscanf(dat," %d\n",&iprint[i]);
            }
        }
        //fscanf(dat,"(Elec, H+, H2+, H20, todas, carga): %d, %d, %d, %d, %d, %d\n",&iprint[0],&iprint[1],&iprint[2],&iprint[3],&iprint[4],&iprint[5]);
        fscanf(dat,"Solo dinamica (S_1/N_0): %d\n",&solodinamica);
        fclose(dat);
        if(solodinamica){
            printf("Solo dinamica es 1!\n");
            pasos = 1000;
            actu = 1;
            terma = 0;
        }
        printf("Dentro del primer if\n");
    }
    //dat = fopen("datos/energia.dat","r");
    if(dat = fopen("energia_inicial.dat","r")){
        dat2 = fopen("datos/energia.dat","w");
        while(fscanf(dat,"\n%i\t%f",&ienergia,&dummy)!=EOF){
            fprintf(dat2,"\n%i\t%e",ienergia,dummy);
        }
        fclose(dat);fclose(dat2);
        printf("Dentro del segundo if\n");
    }
    //getchar();
}
void imprimir_datos_iniciales(void){
    printf("Numero de pasos: %i\n", pasos);
    printf("Actualizacion: %i\n", actu);
    printf("Termalizacion: %i\n", terma);
    printf("Resolucion: %f\n", reso);
    printf("Componentes: %d\n", ncomp);
    printf("Masas (kg):");
    for(int i=0; i<ncomp; i++){
        if(i<ncomp-1){
            printf(" %e,",masa[i]);
        }else{
            printf(" %e\n",masa[i]);
        }
    }
    //masa[0]*=m_u;masa[1]*=m_u;masa[2]*=m_u;masa[3]*=m_u;
    printf("Carga (e):");
    for(int i=0; i<ncomp; i++){
        if(i<ncomp-1){
            printf(" %d,",carga[i]);
        }else{
            printf(" %d\n",carga[i]);
        }
    }
    printf("R (mm): %f\n",R);
    printf("Temperaturas (K):");
    for(int i=0; i<ncomp; i++){
        if(i<ncomp-1){
            printf(" %e,",tempe[i]);
        }else{
            printf(" %e\n",tempe[i]);
        }
    }
    printf("Campo magnetico (T): %e\n",B);
    printf("Constante AE: %f\n",k_auto);
    printf("Tiempo para dinamica (s): %e\n",dtau);
    printf("Iprint:");
    for(int i=0; i<ncomp+2; i++){
        if(i<ncomp+1){
            printf(" %d,",iprint[i]);
        }else{
            printf(" %d\n",iprint[i]);
        }
    }
    //printf("(Elec, H+, H2+, H20, todas, carga): %d, %d, %d, %d, %d, %d\n",iprint[0],iprint[1],iprint[2],iprint[3],iprint[4],iprint[5]);
    printf("Solo dinamica (S_1/N_0): %d\n\n",solodinamica);
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_iniciales(){
    leer_datos_iniciales();
    esc = 1e-3;
    R -= 2.0;
    densidad = 2e19;
    volumen = pi*R*esc*R*esc*1e-6;
    nelectronesr = densidad*volumen;
    npart[0] = nelectronesr;
    npart[1] = npart[0];
    npart[2] = 9*npart[0];
    /*nelectrones = npart[0] = nelectronesr;
    nhp = npart[1] = nelectrones;
    nh2p = npart[2] = 9*nelectrones;//0;//0.2*nelectrones+1;
    nh20 = npart[3] = 9*nelectrones;*/
    printf("\ndensidad: %e volumen: %e nnr: %e\n",densidad,volumen,nelectronesr);
    for(int i=0; i<ncomp; i++)printf("npart_%d: %lld ",i,npart[i]); cout << endl;
    //getchar();

    //tempe[0] = tempee; tempe[1] = tempei; tempe[2] = tempei; tempe[3] = tempei;
    for(int i=0; i<ncomp; i++)tempe[i] *= qe/kb;
    tempei = tempei*qe/kb;
    tempee = tempee*qe/kb;

    for(int i=0; i<ncomp; i++)printf("T_%d: %f ",i,tempe[i]); cout << endl;

    //printf("\ntempei: %f tempee: %f\n", tempei, tempee);
    //getchar();

    contador[0] = 0; contador[1] = 0; contador[2] = 0;
    contador_a[0] = 0; contador_a[1] = 0; contador_a[2] = 0;

    R += 2.0;

    imprimir_datos_iniciales();
}
////////////////////////////////////////////////////////////////////////////////////
void crear_matriz_plasma(){
    //int i,j,k,contadorm=0,contadorm2=11289,contadorvec=0;
    //int i,j,k,contadorm=0,contadorm2=2821,contadorvec=0;
    int i,j,k,contadorm=0,contadorm2=4053,contadorvec=0;
    //int i,j,k,contadorm=0,contadorm2=317,contadorvec=0;
    float vi, fvi;
    printf("Crear matriz plasma rec\n");
    printf("Asignando posiciones\n");
    /*printf("sizeof(plasma[1]): %i\tsizeof(plasma[2]): %i",sizeof(matriz_plasma[1]),sizeof(matriz_plasma[2]));
    printf("sizeof(plasma): %i\tplasma/plasma[1]: %i\n",sizeof(matriz_plasma),sizeof(matriz_plasma)/sizeof(matriz_plasma[1]));
    printf("sizeof(matriz[1][1]): %i\tsizeof(matriz): %i",sizeof(nmatriz_plasma[1][1]),sizeof(nmatriz_plasma));*/
    //getchar();
    for(i=-int(R*reso);i<=int(R*reso);i++){
        for(j=-int(R*reso);j<=int(R*reso);j++){
            //if((i+0.5)*(i+0.5)+(j+0.5)*(j+0.5)<int(R*reso)*int(R*reso)){
            if(i*i+j*j<=int(R*reso)*int(R*reso)){
                contadorm++;
                nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1]=contadorm;
                matriz_plasma[contadorm].xi = i;
                matriz_plasma[contadorm].yi = j;
                matriz_plasma[contadorm].x = (1.0*(i+0.1))/reso;
                matriz_plasma[contadorm].y = (1.0*(j+0.1))/reso;
                matriz_plasma[contadorm].rho = sqrt(i*i+j*j)/reso;
                //matriz_plasma[contadorm].rho = sqrt((matriz_plasma[contadorm].x*reso-0.1)*(matriz_plasma[contadorm].x*reso-0.1)+(matriz_plasma[contadorm].y*reso-0.1)*(matriz_plasma[contadorm].y*reso-0.1));
                matriz_plasma[contadorm].phi = (atan2(matriz_plasma[contadorm].y,matriz_plasma[contadorm].x)+((atan2(matriz_plasma[contadorm].y,matriz_plasma[contadorm].x))<0?2*pi:0.0));
                matriz_plasma[contadorm].anchow = 1.0*esc/reso;
                matriz_plasma[contadorm].anchol = 1e-6;
            }
            else{
                contadorm2++;
                nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1]=contadorm2;
                matriz_plasma[contadorm2].xi = i;
                matriz_plasma[contadorm2].yi = j;
                matriz_plasma[contadorm2].x = (1.0*(i+0.1))/reso;
                matriz_plasma[contadorm2].y = (1.0*(j+0.1))/reso;
                matriz_plasma[contadorm2].rho = sqrt(i*i+j*j)/reso;
                //matriz_plasma[contadorm2].rho = sqrt(matriz_plasma[contadorm2].x*matriz_plasma[contadorm2].x+matriz_plasma[contadorm2].y*matriz_plasma[contadorm2].y);
                matriz_plasma[contadorm2].phi = (atan2(matriz_plasma[contadorm2].y,matriz_plasma[contadorm2].x)+((atan2(matriz_plasma[contadorm2].y,matriz_plasma[contadorm2].x))<0?2*pi:0.0));
                matriz_plasma[contadorm2].anchow = 1.0*esc/reso;
                matriz_plasma[contadorm2].anchol = 1e-6;
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
    n_grupos[0] = n_grupos[1] = n_grupos[2] = 1400*nceldas;
    //n_grupos[0] = n_grupos[1] = n_grupos[2] = 200*nceldas;
    //n_grupos[0] = n_grupos[1] = n_grupos[2] = 100;
    printf("\nnceldas: %i\tn1i: %i\tn2i: %i\n",nceldas,n1i,n2i);
    npart3[0]=npart[0]/n_grupos[0];
    npart3[1]=npart[1]/n_grupos[1];
    npart3[2]=npart[2]/n_grupos[2];
    for(j=0;j<ncomp;j++)printf("npart3_%d: %lld ",j,npart3[j]); cout << endl;

    for(i=1;i<=nceldas;i++){
        vecinos(i);
    }

    /*for(i=1;i<=nceldas;i++){
        if(i%1000==0)printf("\rCelda actual: %i",i);
        matriz_plasma[i].vx[0] = distribucion_normal_5(0.0,1);
        matriz_plasma[i].vy[0] = distribucion_normal_5(0.0,1);
        matriz_plasma[i].vx[1] = distribucion_normal_5(0.0,1);
        matriz_plasma[i].vy[1] = distribucion_normal_5(0.0,1);
        matriz_plasma[i].vx[2] = distribucion_normal_5(0.0,1);
        matriz_plasma[i].vy[2] = distribucion_normal_5(0.0,1);
        matriz_plasma[i].vx[3] = distribucion_normal_5(0.0,1);
        matriz_plasma[i].vy[3] = distribucion_normal_5(0.0,1);

        matriz_plasma[i].v[0] = sqrt( matriz_plasma[i].vx[0]*matriz_plasma[i].vx[0] + matriz_plasma[i].vy[0]*matriz_plasma[i].vy[0] );
        matriz_plasma[i].v[1] = sqrt( matriz_plasma[i].vx[1]*matriz_plasma[i].vx[1] + matriz_plasma[i].vy[1]*matriz_plasma[i].vy[1] );
        matriz_plasma[i].v[2] = sqrt( matriz_plasma[i].vx[2]*matriz_plasma[i].vx[2] + matriz_plasma[i].vy[2]*matriz_plasma[i].vy[2] );
    }*/
    printf("\rImprimiendo archivos de celda\n");
    if( dat = fopen("pruebas/celdas.dat","w") ){
        fprintf(dat,"#X Y Z\n");
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%e\t%e\t%e\n",matriz_plasma[i].x,matriz_plasma[i].y, matriz_plasma[i].rho );
        }
        fclose(dat);
    }
    if( dat = fopen("pruebas/celdas2.dat","w") ){
        fprintf(dat,"i\tj\tx\ty\tncelda\n");
        for(i=1;i<=2*R*reso+1;i++){
            for(j=1;j<=2*R*reso+1;j++){
                fprintf(dat,"%d\t%d\t%f\t%f\t%d\t%f\t%f\n",i,j,matriz_plasma[nmatriz_plasma[i][j]].x,matriz_plasma[nmatriz_plasma[i][j]].y, nmatriz_plasma[i][j], matriz_plasma[nmatriz_plasma[i][j]].anchow, matriz_plasma[nmatriz_plasma[i][j]].anchol);
            }
        }
        fclose(dat);
    }
    ofstream ovecinos;
    ovecinos.open("pruebas/vecinos.dat");
    for(i=1; i<=nceldas; i++){
        ovecinos << i << "\t" << matriz_plasma[i].nvecinos << "\t";
        for(j=1; j<=matriz_plasma[i].nvecinos; j++){
            ovecinos << matriz_plasma[i].vecinos[j] << "\t";
        }
        ovecinos << endl;
    }
    ovecinos.close();
    cout << "\nACABO VECINOS" << endl;
    //getchar(); getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void vecinos(int a){
    int cont_vecinos = 0, ix = floor(matriz_plasma[a].x*reso)+int(R*reso)+1, iy = floor(matriz_plasma[a].y*reso)+int(R*reso)+1, iobj;
    float xx, yy;
    /*cout << "ix: " << ix << " iy: " << iy << endl;
    getchar();*/
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
            xx = matriz_plasma[nmatriz_plasma[ix+i][iy+j]].x;
            yy = matriz_plasma[nmatriz_plasma[ix+i][iy+j]].y;
            iobj = nmatriz_plasma[ix+i][iy+j];
            //printf("a: %4d obj: %4d ix: %3d iy: %3d i: %3d j: %3d\tx: %3.3f\ty: %3.3f\txx: %3.3f\tyy: %3.3f",a,nmatriz_plasma[ix+i][iy+j], ix,iy,i,j,matriz_plasma[a].x,matriz_plasma[a].y,xx,yy);
            //printf("\ta_i: %3d a_j: %3d obj_i: %3d obj_i: %3d",matriz_plasma[a].xi, matriz_plasma[a].yi,matriz_plasma[iobj].xi, matriz_plasma[iobj].yi);
            if(i!=0||j!=0){
                //printf("\tsi");
                if(nmatriz_plasma[ix+i][iy+j]<=nceldas&&nmatriz_plasma[ix+i][iy+j]>0){
                    //printf("\tsi");
                    cont_vecinos++;
                    matriz_plasma[a].vecinos[cont_vecinos]=nmatriz_plasma[ix+i][iy+j];
                }
            }
            //cout << endl;
        }
    }
    matriz_plasma[a].nvecinos = cont_vecinos;
    if(cont_vecinos>9){
        printf("\nnvecinos_%d: %d",a,matriz_plasma[a].nvecinos);
        getchar();
    }
    //getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void calc_carga(int a){
    //long long int carga_i = matriz_plasma[a].carga;
    for(int j=0; j<ncomp; j++)matriz_plasma[a].carga += matriz_plasma[a].part[j]*carga[j];
    /*if(matriz_plasma[a].carga!=0&&p==pasoinicial){
        printf("\nQue show! a: %d carga_i: %lld 0: %lld 1: %lld 2: %lld 3: %lld",a,carga_i,matriz_plasma[a].part[0],matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].part[3]);
        getchar();
    }*/
}
////////////////////////////////////////////////////////////////////////////////////
void arreglo_inicial(void){
    float xrand, yrand, theta, r, dummy, R_inicial = R-4;
    int part_cont = 0;

    printf("\rAsignando arreglo inicial\n");
    if( dat = fopen("posiciones_inicial.dat","r") ){
        fscanf(dat,"%d %d %d\n",&pasoinicial,&termaanterior,&ncomp);
        fclose(dat);
        printf("\npasoinicial: %d\tterma_anterior: %d\tncomp: %d\n",pasoinicial,termaanterior,ncomp);
    }
    //getchar();
    if((pasoinicial>0)&&(ienergia==pasoinicial)){
        printf("\nContinuando posicion final de corrida anterior");
        //getchar();
        if( dat = fopen("posiciones_inicial.dat","r") ){
            fscanf(dat,"%d %d %d %f %f %f %f %d %d",&pasoinicial,&termaanterior,&ncomp,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
            /*if(termaanterior>pasoinicial)terma=termaanterior-pasoinicial;
            else terma = termaanterior;*/
            terma = termaanterior;
            printf("\nterma: %d\tp: %d",terma,p);
            for(int i=1;i<=nceldas;i++){
                fscanf(dat,"\n%i\t%lld\t%lld\t%lld\t%lld\t%lld",&dummy,&matriz_plasma[i].part[0],&matriz_plasma[i].part[3],&matriz_plasma[i].part[1],&matriz_plasma[i].part[2],&matriz_plasma[i].carga);
                if((matriz_plasma[i].part[0]==0)||(matriz_plasma[i].part[3]==0)||(matriz_plasma[i].part[1]==0)||(matriz_plasma[i].part[2]==0)||(matriz_plasma[i].carga==0)){
                    //imprimir_celda_plasma(i);
                }
            }
            fclose(dat);
            printf("\nenergia inicial: %e",energia());
            //getchar();
        }
    }else{
        printf("Arreglo nuevo\n");
        dat = fopen("datos/energia.dat","w");
        for(int j=0; j<ncomp; j++){
            printf("\r\t\t\tj: %d",j);
            part_cont=0;
            do{
                xrand = alea_f(-R_inicial,R_inicial);
                yrand = alea_f(-R_inicial,R_inicial);
                if(xrand*xrand+yrand*yrand<R_inicial*R_inicial){
                //if(part_a_mallado(xrand,yrand)<=nceldas&&part_a_mallado(xrand,yrand)>=1){
                //if( ( xrand*xrand+yrand*yrand < ( R_inicial - 1 )*( R_inicial - 1 ) )&&( xrand*xrand+yrand*yrand > ( R_inicial - 3 )*( R_inicial - 3 ) ) ){
                    if((10*part_cont)%n_grupos[j]==0)printf("\rpart_cont: %d",part_cont);
                    part_cont++;
                    //for(int j=0; j<ncomp; j++){
                        part_plasma[part_cont][j].x = xrand;
                        part_plasma[part_cont][j].y = yrand;
                        //theta = atan2(yrand,xrand); r = norma(xrand,yrand);
                        part_plasma[part_cont][j].vx = distribucion_normal_5(0.0,sqrt(kb*tempe[j]/masa[j]));
                        part_plasma[part_cont][j].vy = distribucion_normal_5(0.0,sqrt(kb*tempe[j]/masa[j]));
                        //part_plasma[part_cont][j].vx = carga[j]*( r*esc*qe*B/masa[j] )*sin(theta);
                        //part_plasma[part_cont][j].vy = -carga[j]*( r*esc*qe*B/masa[j] )*cos(theta);
                        //part_plasma[part_cont][j].vx = alea_i(-1,1)*2*sqrt(kb*tempe[j]/masa[j]);
                        //part_plasma[part_cont][j].vy = alea_i(-1,1)*2*sqrt(kb*tempe[j]/masa[j]);
                        part_plasma[part_cont][j].v = norma(part_plasma[part_cont][j].vx,part_plasma[part_cont][j].vy);
                        part_plasma[part_cont][j].part = npart3[j];
                        if(j==3){
                            printf("\nnpart_%d: %lld part_%d.part: %lld",j,npart3[j],j,part_plasma[part_cont][j].part);
                            getchar();
                        }
                    //}
                }
            }while(part_cont<n_grupos[j]);
            //n_grupos[j]=part_cont;
        }
        /*for(int j=0; j<ncomp; j++){
            char part_char[20];
            sprintf(part_char,"pruebas/part_%d.dat",j);
            dat = fopen(part_char,"w");
            for(int i=1; i<= n_grupos[j]; i++){
                fprintf(dat,"%d %f %f %f %f\n",i,part_plasma[i][j].x,part_plasma[i][j].y,part_plasma[i][j].vx,part_plasma[i][j].vy);
            }
            fclose(dat);
        }*/
        cout << "\nDespues de part" << endl;
        importar_part_a_mallado();
        cout << "Despues de importar" << endl;
        long long int cargatotaaal = 0;
        for(int i=1;i<=nceldas;i++){
            cargatotaaal += matriz_plasma[i].carga;
        }
        if(debug==2)cout << endl << "QUE SHOW" << endl;
        if(cargatotaaal!=0){
            printf("\ncargatotal: %lld",cargatotaaal);
            getchar();
        }
        p = 0;
        actu_salida();
        for(int j=0; j<ncomp; j++)hacer_histograma(j,0);
        printf("\nterma: %d ASDASD", terma);
        //getchar();
        /*printf("\nCalculando carga");
        for(int i=1;i<=nceldas;i++){
            calc_carga(i);
            if(matriz_plasma[i].carga!=0){
                printf("\nQue show i: %d carga: %lld",i,matriz_plasma[i].carga );
                getchar();
            }
        }*/
    }
}
////////////////////////////////////////////////////////////////////////////////////
void importar_part_a_mallado(void){
    int icelda;
    for(int i=1; i<=nceldas; i++){
        for(int j=0; j<ncomp; j++)matriz_plasma[i].part[j] = 0;
        matriz_plasma[i].carga = 0;
    }
    for(int j=0; j<ncomp; j++){
        for(int i=1;i<=n_grupos[j];i++){
            icelda = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y);
            matriz_plasma[icelda].part[j]+=part_plasma[i][j].part;
            matriz_plasma[icelda].carga+=carga[j]*part_plasma[i][j].part;
            //printf("\ni: %d j: %d p_x: %f p_y: %f",i, j, part_plasma[i][j].x,part_plasma[i][j].y);
            //printf("\nint_x: %d int_y: %d",int_x, int_y);
            //printf("\nnmatriz_i,j: %d mat_x: %f mat_y: %f",nmatriz_plasma[int_x][int_y],matriz_plasma[nmatriz_plasma[int_x][int_y]].x,matriz_plasma[nmatriz_plasma[int_x][int_y]].y);
            //getchar();
            /*if((int_x>R_reso*R_reso+1)||(int_y>R_reso*R_reso+1)||(int_x<=0)||(int_y<=0)){
                printf("\nint_x: %d int_y: %d",int_x,int_y);
                getchar();
            }*/
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////
//int part_a_mallado(int a, int b){
int part_a_mallado(float a, float b){
    int ix, iy;
    ix = floor( a*reso ) + floor(R*reso) + 1;
    iy = floor( b*reso ) + floor(R*reso) + 1;
    return(nmatriz_plasma[ix][iy]);
}
////////////////////////////////////////////////////////////////////////////////////
void hacer_histograma(int a, int b){
    const int int_clases = 200;
    long long int clasesx[int_clases+1]={0},clasesy[int_clases+1]={0};
    int i;
    char nombre[50];
    float vxmin = part_plasma[1][a].vx, vxmax = part_plasma[1][a].vx, vymin = part_plasma[1][a].vy, vymax = part_plasma[1][a].vy;
    for(i=1;i<=int_clases;i++){
        clasesx[i]=clasesy[i]=0;
    }
    /*spritf(nombre,"histogramas/vx%d_%d.dat",a,b);
    dat = fopen(nombre,"w");
    sprintf(nombre,"histogramas/vy%d_%d.dat",a,b);
    dat2 = fopen(nombre,"w");
    fprintf(dat,"#ncelda\tvx\n");
    fprintf(dat2,"#ncelda\tvy\n");*/
    for(i=1;i<=nceldas;i++){
        if(part_plasma[i][a].vx<vxmin)vxmin=part_plasma[i][a].vx;
        if(part_plasma[i][a].vx>vxmax)vxmax=part_plasma[i][a].vx;
        if(part_plasma[i][a].vy<vymin)vymin=part_plasma[i][a].vy;
        if(part_plasma[i][a].vy>vymax)vymax=part_plasma[i][a].vy;
        //fprintf(dat,"%d\t%e\n",i,matriz_plasma[i].vx[a]);
        //fprintf(dat2,"%d\t%e\n",i,matriz_plasma[i].vy[a]);
    }
    vxmin-=1;vxmax+=1;
    vymin-=1;vymax+=1;
    //fclose(dat);fclose(dat2);
    for(i=1;i<=nceldas;i++){
            if(int( int_clases*(matriz_plasma[i].vx[a]-vxmin)/(vxmax-vxmin) )+1<0){
                printf("\nVXX: %e\tvxmin: %e\tvxmax: %e",part_plasma[i][a].vx,vxmin,vxmax);
                getchar();
            }
            if(int( int_clases*(matriz_plasma[i].vy[a]-vymin)/(vymax-vymin) )+1<0){
                printf("\nVYY: %e\tvymin: %e\tvymax: %e",part_plasma[i][a].vy,vymin,vymax);
                getchar();
            }
        clasesx[ int( int_clases*(part_plasma[i][a].vx-vxmin)/(vxmax-vxmin) )+1]+=part_plasma[i][a].part;
        clasesy[ int( int_clases*(part_plasma[i][a].vy-vymin)/(vymax-vymin) )+1]+=part_plasma[i][a].part;
    }
    sprintf(nombre,"histogramas/histv_part_%d_%d.dat",a,b);
    if( dat = fopen(nombre,"w") ){
        fprintf(dat,"#X\tY\n");
        for(i=1;i<=int_clases;i++){
            /*if(clases[i]>300){
                printf("\nQue show! clases[%d]: %d",i,clases[i]);
                getchar();
            }*/
            fprintf(dat,"%e\t%e\t%e\t%e\n", vxmin + (i-1)*(vxmax-vxmin)/int_clases , 1.0*clasesx[i]/n_grupos[a], vymin + (i-1)*(vymax-vymin)/int_clases , 1.0*clasesy[i]/n_grupos[a] );
            //fprintf(dat,"%e\t%e\t%e\t%e\n", vxmin + (i-1)*(vxmax-vxmin)/int_clases , 1.0*clasesx[i], vymin + (i-1)*(vymax-vymin)/int_clases , 1.0*clasesy[i] );
        }
        fclose(dat);
    }
    /*for(i=1;i<=int_clases;i++){
        clasesx[i] = clasesy[i] = 0;
    }
    for(i=1;i<=nceldas;i++){
        clasesx[ int( int_clases*(matriz_plasma[i].vx_1[a]-vxmin)/(vxmax-vxmin) )+1]+=matriz_plasma[i].part[a];
        clasesx[ int( int_clases*(matriz_plasma[i].vx_2[a]-vxmin)/(vxmax-vxmin) )+1]+=matriz_plasma[i].part[a];
        clasesy[ int( int_clases*(matriz_plasma[i].vy_1[a]-vymin)/(vymax-vymin) )+1]+=matriz_plasma[i].part[a];
        clasesy[ int( int_clases*(matriz_plasma[i].vy_2[a]-vymin)/(vymax-vymin) )+1]+=matriz_plasma[i].part[a];
    }
    sprintf(nombre,"histogramas/histv2_%d_%d.dat",a,b);
    if( dat = fopen(nombre,"w") ){
        fprintf(dat,"#X\tY\n");
        for(i=1;i<=int_clases;i++){
            fprintf(dat,"%e\t%e\t%e\t%e\n", vxmin + (i-1)*(vxmax-vxmin)/int_clases , 1.0*clasesx[i]/npart[a], vymin + (i-1)*(vymax-vymin)/int_clases , 1.0*clasesy[i]/npart[a] );
        }
        fclose(dat);
    }*/
}
////////////////////////////////////////////////////////////////////////////////////
void hacer_distribucion( int a ){
    const int int_clases = 100;
    int i,j;
    long long int distrho_carga[int_clases+1] = {0}, distrho_part[int_clases+1] = {0}, distphi_carga[int_clases+1] = {0}, distphi_part[int_clases+1] = {0};
    float rhomin=-1e-3, rhomax=R+1e-3, phimin = -1e-3, phimax = 2*pi+1e-3;
    char nombre[50];

    /*for(i=1;i<=nceldas;i++){
        distrho_carga[ int( int_clases*(matriz_plasma[i].rho-rhomin)/(rhomax-rhomin) )+1 ]+=-matriz_plasma[i].part[0]+matriz_plasma[i].part[1]+matriz_plasma[i].part[2];
        distrho_part[ int( int_clases*(matriz_plasma[i].rho-rhomin)/(rhomax-rhomin) )+1 ]+=matriz_plasma[i].part[0];//+matriz_plasma[i].part[1]+matriz_plasma[i].part[2]+matriz_plasma[i].part[3];
        distphi_carga[ int( int_clases*(matriz_plasma[i].phi-phimin)/(phimax-phimin) )+1 ]+=-matriz_plasma[i].part[0]+matriz_plasma[i].part[1]+matriz_plasma[i].part[2];
        distphi_part[ int( int_clases*(matriz_plasma[i].phi-phimin)/(phimax-phimin) )+1 ]+=matriz_plasma[i].part[0];//+matriz_plasma[i].part[1]+matriz_plasma[i].part[2]+matriz_plasma[i].part[3];
    }*/
    for(i=1;i<=nceldas;i++){
        distrho_carga[ int( int_clases*(matriz_plasma[i].rho-rhomin)/(rhomax-rhomin) )+1 ]+=-gelectron[i]+ghp[i]+gh2p[i];
        distrho_part[ int( int_clases*(matriz_plasma[i].rho-rhomin)/(rhomax-rhomin) )+1 ]+=gelectron[i]+ghp[i]+gh2p[i];//+matriz_plasma[i].part[1]+matriz_plasma[i].part[2]+matriz_plasma[i].part[3];
        distphi_carga[ int( int_clases*(matriz_plasma[i].phi-phimin)/(phimax-phimin) )+1 ]+=-gelectron[i]+ghp[i]+gh2p[i];
        distphi_part[ int( int_clases*(matriz_plasma[i].phi-phimin)/(phimax-phimin) )+1 ]+=gelectron[i]+ghp[i]+gh2p[i];//+matriz_plasma[i].part[1]+matriz_plasma[i].part[2]+matriz_plasma[i].part[3];
    }
    sprintf(nombre,"datos/dist_%d.dat",a);
    if( dat = fopen(nombre,"w") ){
        fprintf(dat,"#rho\tcarga\tparticula\tphi\tcarga\tparticula\n");
        for(i=1;i<=int_clases;i++){
            /*if(clases[i]>300){
                printf("\nQue show! clases[%d]: %d",i,clases[i]);
                getchar();
            }*/
            fprintf(dat,"%f\t%15lld\t%15lld\t", rhomin + (i-1)*(rhomax-rhomin)/int_clases , distrho_carga[i]/(p-terma), distrho_part[i]/(p-terma) );
            fprintf(dat,"%f\t%15lld\t%15lld\n", phimin + (i-1)*(phimax-phimin)/int_clases , distphi_carga[i]/(p-terma), distphi_part[i]/(p-terma) );
        }
        fclose(dat);
    }
    /*sprintf(nombre,"datos/distphi_%d.dat",a);
    if( dat = fopen(nombre,"w") ){
        fprintf(dat,"#phi\tcarga\tparticula\n");
        for(i=1;i<=int_clases;i++){
            //if(clases[i]>300){
            //    printf("\nQue show! clases[%d]: %d",i,clases[i]);
            //    getchar();
            //}
            fprintf(dat,"%f\t%15I64d\t%15I64d\n", phimin + (i-1)*(phimax-phimin)/int_clases , distphi_carga[i]/(p-terma), distphi_part[i]/(p-terma) );
        }
        fclose(dat);
    }*/
    sprintf(nombre,"datos/dens_%d.dat",a);
    if( dat = fopen(nombre,"w") ){
        fprintf(dat,"#rho\tcarga\tparticula\tphi\tcarga\tparticula\n");
        for(i=1;i<=int_clases;i++){
            /*if(clases[i]>300){
                printf("\nQue show! clases[%d]: %d",i,clases[i]);
                getchar();
            }*/
            fprintf(dat,"%f\t%e\t%e\t", rhomin + (i-1)*(rhomax-rhomin)/int_clases , (distrho_carga[i]*int_clases*int_clases)/((p-terma)*pi*(i*i-(i-1)*(i-1))*(rhomax-rhomin)*(rhomax-rhomin)), (distrho_part[i]*int_clases*int_clases)/((p-terma)*pi*(i*i-(i-1)*(i-1))*(rhomax-rhomin)*(rhomax-rhomin)) );
            fprintf(dat,"%f\t%e\t%e\n", phimin + (i-1)*(phimax-phimin)/int_clases , (distphi_carga[i]*int_clases)/((p-terma)*pi*R*R), (distphi_part[i]*int_clases)/((p-terma)*pi*R*R) );
        }
        fclose(dat);
    }
}
////////////////////////////////////////////////////////////////////////////////////
void mover_particulas_part(void){
    int idx, idy;
    int c_pared = 0;
    tipo = alea_i( 0, ncomp-2 );
    //printf("\ntipo: %d contador_tipo: %d",tipo,contador[tipo]);
    contador[tipo]++;
    //printf("\ncontador_tipo: %d",contador[tipo]);
    alea_part = alea_i( 1, n_grupos[tipo] );

    n2 = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y );
    matriz_plasma[n2i]=matriz_plasma[n2];
    part_plasma[0][tipo] = part_plasma[alea_part][tipo];
    do{
        c_pared++;
        if(rechazo==1)part_plasma[alea_part][tipo]=part_plasma[0][tipo];
        do{
            idx = alea_i(-1,1);
            idy = (idx==0)?(-1+2*alea_i(0,1)):(alea_i(-1,1));
            //idy = alea_i(-1,1);
            if(c_pared>100){
                printf("\nidx: %d idy: %d",idx,idy);
            }
        }while(idx==0&&idy==0);
        //printf("\nANTES tipo: %d a_part: %d px_0: %f py_0: %f mx: %f my: %f",tipo,alea_part,part_plasma[0][tipo].x,part_plasma[0][tipo].y,matriz_plasma[n2].x,matriz_plasma[n2].y);
        part_plasma[alea_part][tipo].x += idx/reso;
        part_plasma[alea_part][tipo].y += idy/reso;
        //part_plasma[alea_part][tipo].x = alea_f(-R,R);
        //part_plasma[alea_part][tipo].y = alea_f(-R,R);
        n1 = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y );
        matriz_plasma[n1i]=matriz_plasma[n1];
        rechazo = 0;
        if(norma(part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y)>=(R-2))rechazo=1;
        //if(part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)>nceldas||part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)<1)rechazo=1;
        //pared();
        if(c_pared>100){
            imprimir_celda_plasma_part(0,tipo);
            imprimir_celda_plasma_part(alea_part,tipo);
            printf("\ncontador_pared: %d rechazo: %d norma: %f",c_pared,rechazo, norma(part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y),R);
            getchar();
        }
        if(c_pared>110)rechazo=0;
    }while(rechazo==1);
    if(c_pared>1000){
        printf("\nc_pared: %d tipo: %d a_part: %d x: %f y: %f rho: %f",c_pared,tipo,alea_part,part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y,norma(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y));
        printf("\nx_0: %f y_0: %f rho_0: %f",part_plasma[0][tipo].x,part_plasma[0][tipo].y,norma(part_plasma[0][tipo].x,part_plasma[0][tipo].y));
        getchar();
    }
    c_pared = 0;
    rechazo = 0;
    //printf("\nDESPUES px_f: %f py_f: %f mx_f: %f my_f: %f",part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y,matriz_plasma[n1].x,matriz_plasma[n1].y);
    //getchar();
    matriz_plasma[n1].part[ncomp]-=part_plasma[alea_part][tipo].part;
    matriz_plasma[n1].carga-=carga[tipo]*part_plasma[alea_part][tipo].part;
    matriz_plasma[n2].part[ncomp]+=part_plasma[alea_part][tipo].part;
    matriz_plasma[n2].carga+=carga[tipo]*part_plasma[alea_part][tipo].part;
    //importar_part_a_mallado();
}
////////////////////////////////////////////////////////////////////////////////////
void pared(void){
    float rho = norma(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y);
    //if(part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)>nceldas||part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)<1)rechazo=1;
    if(rho>=(R-2)){
        rechazo=1;
        //printf("\nQue show!");
        //getchar();
    }
}
////////////////////////////////////////////////////////////////////////////////////
void reflejar_pared(int a,int b,float* vx, float* vy){
    float theta, v0x = part_plasma[a][b].vx, v0y = part_plasma[a][b].vy, x0 = part_plasma[a][b].x, y0 = part_plasma[a][b].y;
    /*do{
        v0x = 0.2*alea_f(-1,1); v0y = 0.2*alea_f(-1,1);
        x0 = alea_f(-1,1); y0 = (-1+2*alea_i(0,1))*sqrt(1-x0*x0);
    }while((v0x*cos(theta)+v0y*sin(theta))<0);*/
    theta = atan2(y0,x0);
    //printf("\na: %d b: %d x: %f y: %f v0x: %f v0y: %f",a, b, x0, y0, v0x, v0y);
    //part_plasma[a][b].vx = v0x - 2*(v0x*cos(theta)+v0y*sin(theta))*cos(theta);
    //part_plasma[a][b].vy = v0y - 2*(v0x*cos(theta)+v0y*sin(theta))*sin(theta);
    *vx = v0x - 2*(v0x*cos(theta)+v0y*sin(theta))*cos(theta);
    *vy = v0y - 2*(v0x*cos(theta)+v0y*sin(theta))*sin(theta);
    //printf("\na: %d b: %d vxf: %f vyf: %f",a, b, part_plasma[a][b].vx, part_plasma[a][b].vy);
    //fprintf(dat,"%f\t%f\t%f\t%f\n",x0-v0x/norma(v0x,v0y),y0-v0y/norma(v0x,v0y),v0x/norma(v0x,v0y),v0y/norma(v0x,v0y));
    //fprintf(dat2,"%f\t%f\t%f\t%f\n",x0,y0,part_plasma[a][b].vx/norma(part_plasma[a][b].vx,part_plasma[a][b].vy),part_plasma[a][b].vy/norma(part_plasma[a][b].vx,part_plasma[a][b].vy));
    //getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda_plasma_rec(int a){
    printf("\nCelda: %i x: %f y: %f\nelectrones: %lld h20: %lld hp: %lld h2p: %lld carga: %lld",a,matriz_plasma[a].x,matriz_plasma[a].y,matriz_plasma[a].part[0],matriz_plasma[a].part[3],matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].carga);
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda_plasma_part(int a, int b){
    printf("\nPart: %d Comp: %d x: %f y: %f rho: %f vx: %f vy: %f\npart: %lld",a,b,part_plasma[a][b].x,part_plasma[a][b].y,norma(part_plasma[a][b].x,part_plasma[a][b].y),part_plasma[a][b].vx,part_plasma[a][b].vy, part_plasma[a][b].part);
}
////////////////////////////////////////////////////////////////////////////////////
float distancia(int a, int b){
    float dist;
    dist = sqrt(pow(matriz_plasma[a].x-matriz_plasma[b].x,2)+pow(matriz_plasma[a].y-matriz_plasma[b].y,2));
    return(dist);
}
////////////////////////////////////////////////////////////////////////////////////
int signo(float a){
    /*int sign;
    if(a>=0){
        sign = 1;
    }
    else{
        sign = -1;
    }*/
    return(a>=0?1:-1);
}
////////////////////////////////////////////////////////////////////////////////////
float de_plasma(void){
    int i;
    bool bool1 = false;
    //bool1 = true;
    double ei = 0, ef = 0, d;
    for(i=1;i<=nceldas;i++){
        if(matriz_plasma[i].carga!=0){
            if((i!=n1)&&(i!=n2)){
                //d = distanciacelda(n1i,i);
                d = distancia(n1i,i);
                ei += (qe*qe*matriz_plasma[n1i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                //d = distanciacelda(n1,i);
                d = distancia(n1,i);
                ef += (qe*qe*matriz_plasma[n1].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                /*ei += 8*ecoulomb(n1i,i,1);
                ef += 8*ecoulomb(n1,i,1);*/
                //d = distanciacelda(n2i,i);
                d = distancia(n2i,i);
                ei += (qe*qe*matriz_plasma[n2i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                //d = distanciacelda(n2,i);
                d = distancia(n2,i);
                ef += (qe*qe*matriz_plasma[n2].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                /*ei += 8*ecoulomb(n2i,i,1);
                ef += 8*ecoulomb(n2,i,1);*/
            }
            //getchar();
        }
    }
    //d = distanciacelda(n1i,n2i);
    d = distancia(n1i,n2i);
    ei += (qe*qe*matriz_plasma[n1i].carga*matriz_plasma[n2i].carga)/(4*pi*epce*epsi*d*esc);
    //ei += 8*ecoulomb(n1i,n2i,1);
    //d = distanciacelda(n1,n2);
    d = distancia(n1,n2);
    ef += (qe*qe*matriz_plasma[n1].carga*matriz_plasma[n2].carga)/(4*pi*epce*epsi*d*esc);
    // += 8*ecoulomb(n1,n2,1);
    /*for(i=1;i<=8;i++){
        ei += ecoulomb(n1i,n1i,i+1) + ecoulomb(n2i,n2i,i+1);
        ef += ecoulomb(n1,n1,i+1) + ecoulomb(n2,n2,i+1);
    }*/
    de_coulomb = ef - ei;
    if(p>=terma){
        /*if(fabs(autoenergia(n1i))>fabs(ei)||fabs(autoenergia(n2i))>fabs(ei)||fabs(autoenergia(n1))>fabs(ef)||fabs(autoenergia(n2))>fabs(ef)){
            printf("\ntipo: %i ef: %e ei: %e\naE_n1i: %e aE_n2i: %e aE_n1: %e aE_n2: %e",tipo,ef,ei,autoenergia(n1i),autoenergia(n2i),autoenergia(n1),autoenergia(n2));
            bool1 = true;
        }*/
        /*printf("\npaso: %d tipo: %d n1: %d n2: %d n1.q: %lld n2.q: %lld E_n1n2: %e AE_n1: %e AE_n2: %e",p,tipo,n1,n2,matriz_plasma[n1].carga,matriz_plasma[n2].carga,(qe*qe*matriz_plasma[n1].carga*matriz_plasma[n2].carga)/(4*pi*epce*epsi*distancia(n1,n2)*esc),autoenergia(n1),autoenergia(n2));
        printf("\nn1i.q: %lld n2i.q: %lld",matriz_plasma[n1i].carga,matriz_plasma[n2i].carga);*/
        //bool1 = true;
        ei += autoenergia(n1i);
        ei += autoenergia(n2i);
        ef += autoenergia(n1);
        ef += autoenergia(n2);
        if(bool1){
            printf(" de: %e",ef-ei);
            getchar();
        }
    }
    de_autoenergia = autoenergia(n1) + autoenergia(n2) - autoenergia(n1i) - autoenergia(n2i);
    dem = ef - ei;
    //printf("\nde_c: %e de_ae: %e de_a: %e T: %e kT: %e",de_coulomb,de_autoenergia,dem,tempe[tipo],kb*tempe[tipo]);
    //getchar();
    //printf("\n plasma eic: %e efc: %e de: %e",ei,ef,dem);
    //getchar();
    return(dem);
}
////////////////////////////////////////////////////////////////////////////////////
float autoenergia(int a){
    float aenergia;
    //aenergia = k_auto*(6*matriz_plasma[a].carga*matriz_plasma[a].carga*qe*qe)/(5*4*pi*epce*epsi*matriz_plasma[a].anchow);
    //k_auto = 1.4866047991;//2.0*log(1+sqrt(2))+2.0*(1-sqrt(2))/3.0;
    aenergia = (k_auto*qe*qe*matriz_plasma[a].carga*matriz_plasma[a].carga)/(4.0*pi*epce*epsi*matriz_plasma[a].anchow);
    //printf("\na: %i\t w: %f",a,matriz_plasma[a].anchow);
    //printf("\na: %i\t%I64d\t%e",a,matriz_plasma[a].carga,aenergia);
    //getchar();
    return(aenergia);
}
////////////////////////////////////////////////////////////////////////////////////
void calcular_momento_total(void){
    int i;
    double ptotalx = 0,ptotaly = 0;
    for(i=1;i<=nceldas;i++){
        ptotalx += matriz_plasma[i].part[0]*masa[0]*matriz_plasma[i].vx[0] + matriz_plasma[i].part[1]*masa[1]*matriz_plasma[i].vx[1] + matriz_plasma[i].part[2]*masa[2]*matriz_plasma[i].vx[2] + matriz_plasma[i].part[3]*masa[3]*matriz_plasma[i].vx[3];
        ptotaly += matriz_plasma[i].part[0]*masa[0]*matriz_plasma[i].vy[0] + matriz_plasma[i].part[1]*masa[1]*matriz_plasma[i].vy[1] + matriz_plasma[i].part[2]*masa[2]*matriz_plasma[i].vy[2] + matriz_plasma[i].part[3]*masa[3]*matriz_plasma[i].vy[3];
    }
    if(dat = fopen("momento_total.dat","a")){
        fprintf(dat,"%d\t%e\t%e\n",p-1,ptotalx,ptotaly);
        fclose(dat);
    }
}
////////////////////////////////////////////////////////////////////////////////////
void metropolis_plasma(void){
    float zeta, argexp, emet;
    double de;
    int i, j;
    zeta = alea();
    /*if(a==0){
        argexp = -de_plasma()/(kb*tempe[a]);
    }else if(a==1){
        argexp = -de_plasma()/(kb*tempe[a]);
    }else if(a==2){
        printf("\nQue show:v");
        getchar();
        argexp = -de_plasma()/(kb*tempe[a]);
    }else{
        argexp = -de_plasma()/(kb*tempe[a]);
    }*/
    de = de_plasma();
    argexp = -de/(kb*tempe[tipo]);
    //argexp = 0;
    /*if(a==1){
        argexp = -de_plasma()/(kb*tempee);
    }
    else if(a==2){
        argexp = -de_plasma()/(kb*tempei);
    }
    else{
        argexp = -de_plasma()/(kb*tempee);
        //if(jaja==-1){
            //printf("\ntipo3 de: %e nale: %I64d\n",de_plasma(),nale);
            //printf("n1: %i n2: %i n1i: %i n2i: %i\n",n1,n2,n1i,n2i);
            //printf("\ncarga_1: %I64d carga_2: %I64d carga_1i: %I64d carga_2i: %I64d",matriz_plasma[n1].carga,matriz_plasma[n2].carga,matriz_plasma[n1i].carga,matriz_plasma[n2i].carga);
            //getchar();
        //}
    }*/

    /*if(argexp>=100)
    {
        emet = 2.0;
    }
    else
    {
        emet = exp(argexp);
    }*/

    if((argexp>=100?2.0:exp(argexp))>=zeta){
        contador_a[tipo]++;
        if((100*p)%actu==0){
            dat = fopen("pruebas/denergias_acep.dat","a");
            fprintf(dat,"%d\t%e\t%e\t%e\n",p,de_coulomb/(kb*tempe[tipo]),de_autoenergia/(kb*tempe[tipo]),dem/(kb*tempe[tipo]));
            fclose(dat);
        }
        dt_termico+=(esc*distancia(n1,n2)*part_plasma[alea_part][tipo].part)/part_plasma[alea_part][tipo].v;
        //printf("\np: %d d: %f part_desp: %lld",p,distancia(n1,n2),part_plasma[alea_part][tipo].part);
        //getchar();
        /*if(a==1){
            contador_a[0]++;
        }
        else if(a==2){
            contador_a[1]++;
        }
        else{
            contador_a[2]++;
        }*/
    }
    else{
        //matriz_plasma[n1] = matriz_plasma[n1i];
        //matriz_plasma[n2] = matriz_plasma[n2i];
        if((100*p)%actu==0){
            dat = fopen("pruebas/denergias_rechazo.dat","a");
            fprintf(dat,"%d\t%e\t%e\t%e\n",p,de_coulomb/(kb*tempe[tipo]),de_autoenergia/(kb*tempe[tipo]),dem/(kb*tempe[tipo]));
            fclose(dat);
        }
        part_plasma[alea_part][tipo] = part_plasma[0][tipo];
        matriz_plasma[n1] = matriz_plasma[n1i];
        matriz_plasma[n2] = matriz_plasma[n2i];
        //importar_part_a_mallado();
        rechazo_met_mov++;
    }
}
////////////////////////////////////////////////////////////////////////////////////
void actu_salida(void){
	int max_n_grupos = n_grupos[0];
	long long int cargatotaaal=0,cargamin=matriz_plasma[1].carga,cargamax=matriz_plasma[1].carga;
	if(p%(10*actu)==0&&p<-5){
        for(int i=0; i<ncomp; i++ )if(n_grupos[i]>=max_n_grupos)max_n_grupos = n_grupos[i];
        sprintf(salidac,"datos/posiciones%i.dat",(p-terma)/actu);
        if( dat=fopen(salidac,"w") ){
            fprintf(dat,"%d %d %d %f %f %f %f %d %d",p,terma,ncomp,((contador[0]>0)?1.0*contador_a[0]/contador[0]:0.0),((contador[1]>0)?1.0*contador_a[1]/contador[1]:0.0),((contador[2]>0)?1.0*contador_a[2]/contador[2]:0.0),((contador[3]>0)?1.0*contador_a[3]/contador[3]:0.0),contador_rechazo,rechazo_met_mov);
            for(int i=1;i<=max_n_grupos;i++){
                fprintf(dat, "\n%d ",i);
                for(int j=0; j<ncomp; j++){
                    if(i<=n_grupos[j])fprintf(dat,"%f %f %f %f %lld ",part_plasma[i][j].x, part_plasma[i][j].y, part_plasma[i][j].vx, part_plasma[i][j].vy, part_plasma[i][j].part);
                }
            }
            fclose(dat);
        }
	}
    for(int j=0; j<ncomp; j++){
        sprintf(salidac,"datos/part_%d_%d.dat",j,(p-terma)/actu);
        if( iprint[j]==1 ){
            dat=fopen(salidac,"w");
        //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
            fprintf(dat, "#X\tY\tZ\n");
            for(int i=1;i<=nceldas;i++){
                fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[j]);
            }
            fclose(dat);
        }
    }
    sprintf(salidac,"datos/todas%i.dat",(p-terma)/actu);
    if( iprint[ncomp]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(int i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[0]+matriz_plasma[i].part[1]+matriz_plasma[i].part[2]);
        }
        fclose(dat);
    }

    sprintf(salidac,"datos/carga%i.dat",(p-terma)/actu);
    if( iprint[ncomp+1]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"carga\"\n");
        /*fprintf(dat, "#X\tY\tZ\n");
        for(int i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%e\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].carga*qe/(matriz_plasma[i].anchow*matriz_plasma[i].anchow));
            cargatotaaal += matriz_plasma[i].carga;
        }
        fprintf(dat,"carga total: %lld\n",cargatotaaal);
        fclose(dat);
        dat=fopen("datos/cargatotal.dat","a");
        fprintf(dat,"%d\t%lld\n",p,cargatotaaal);
        fclose(dat);*/


        int ix, iy, numerocelda;
        float xx, yy;
        for(int i=-int(R*reso)-1;i<=int(R*reso);i++){
            for(int j=-int(R*reso)-1;j<=int(R*reso);j++){
                ix = j + int(R*reso) + 1;
                iy = i + int(R*reso) + 1;
                numerocelda = nmatriz_plasma[ix][iy];
                xx = (1.0*(j+0.0))/reso;
                yy = (1.0*(i+0.0))/reso;
                if(ix==0){
                    if(iy==0){
                        fprintf(dat,"0 ");
                    }else{
                        if((iy-1)%(2*int(reso))==0||iy==1||iy==(2*int(R*reso)+1))fprintf(dat,"%1.0f ",yy);
                        else fprintf(dat,"? ");
                    }
                }else{
                    if(iy==0){
                        if((ix-1)%(2*int(reso))==0||ix==1||ix==(2*int(R*reso)+1))fprintf(dat,"%1.0f ",xx);
                        else fprintf(dat,"? ");
                    }else{
                        //if(i*i+j*j<=int(R*reso)*int(R*reso)){
                        if(numerocelda<=nceldas&&numerocelda>=1){
                            fprintf(dat,"%e ",matriz_plasma[ numerocelda ].carga*qe/(matriz_plasma[ numerocelda ].anchow*matriz_plasma[ numerocelda ].anchow));
                            cargatotaaal+=matriz_plasma[ numerocelda ].carga;
                            if(matriz_plasma[ numerocelda ].carga>cargamax)cargamax=matriz_plasma[ numerocelda ].carga;
                            if(matriz_plasma[ numerocelda ].carga<cargamin)cargamin=matriz_plasma[ numerocelda ].carga;
                        }
                        else{
                            fprintf(dat,"? ");
                        }
                    }
                }
            }
            fprintf(dat,"\n");
        }
        fclose(dat);

        sprintf(salidac,"datos/carga%i_minmax.dat",(p-terma)/actu);
        dat = fopen(salidac,"w");
        fprintf(dat,"1 %e\n2 %e",cargamax*qe/(matriz_plasma[ numerocelda ].anchow*matriz_plasma[ numerocelda ].anchow),cargamin*qe/(matriz_plasma[ numerocelda ].anchow*matriz_plasma[ numerocelda ].anchow));
        fclose(dat);
    }
}
////////////////////////////////////////////////////////////////////////////////////
void actu_salida_carga(void){
	long long int cargatotaaal=0;

    sprintf(salidac,"pruebas/prueba_carga.dat");
    dat=fopen(salidac,"w");
//fprintf(dat, "\"x\", \"y\", \"carga\"\n");
    /*fprintf(dat, "#X\tY\tZ\n");
    for(int i=1;i<=nceldas;i++){
        fprintf(dat,"%f\t%f\t%e\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].carga*qe/(matriz_plasma[i].anchow*matriz_plasma[i].anchow));
    }*/
    int ix, iy, numerocelda;
    float xx, yy;
    for(int i=-int(R*reso)-1;i<=int(R*reso);i++){
        for(int j=-int(R*reso)-1;j<=int(R*reso);j++){
            ix = j + int(R*reso) + 1;
            iy = i + int(R*reso) + 1;
            numerocelda = nmatriz_plasma[ix][iy];
            xx = (1.0*(j+0.0))/reso;
            yy = (1.0*(i+0.0))/reso;
            if(ix==0){
                if(iy==0){
                    fprintf(dat,"0 ");
                }else{
                    if((iy)%(3*int(reso))==0||iy==1||iy==(2*int(R*reso)+1))fprintf(dat,"%1.3f ",yy);
                    else fprintf(dat,"? ");
                }
            }else{
                if(iy==0){
                    if((ix)%(3*int(reso))==0||ix==1||ix==(2*int(R*reso)+1))fprintf(dat,"%1.3f ",xx);
                    else fprintf(dat,"? ");
                }else{
                    //if(i*i+j*j<=int(R*reso)*int(R*reso)){
                    if(numerocelda<=nceldas&&numerocelda>=1){
                        fprintf(dat,"%e ",matriz_plasma[ numerocelda ].carga*qe/(matriz_plasma[ numerocelda ].anchow*matriz_plasma[ numerocelda ].anchow));
                        cargatotaaal+=matriz_plasma[ numerocelda ].carga;
                    }
                    else{
                        fprintf(dat,"? ");
                    }
                }
            }
        }
        fprintf(dat,"\n");
    }
    fclose(dat);
    printf("\ncargatotal: %lld",cargatotaaal);
}
////////////////////////////////////////////////////////////////////////////////////
void salida_prom(void){
	//long long int cargatotaaal=0;
	//float xx, yy;
	long long gtotal[11260]={0}, gcarga[11260]={0};
    for(int j=0;j<ncomp;j++){
        sprintf(salidac,"datos/promedios/gpart_%d.dat",j);
        dat=fopen(salidac,"w");
        //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(int i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, gpart[i][j]/(p-terma));
        }
        fclose(dat);
    }

    for(int i=1; i<=nceldas; i++){
        for(int j=0; j<ncomp; j++){
            gtotal[i]+=gpart[i][j];
            gcarga[i]+=carga[j]*gpart[i][j];
        }
    }
    sprintf(salidac,"datos/promedios/todas.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(int i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, gtotal[i]/(p-terma));
        //}
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/carga.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"carga\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(int i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, gcarga[i]/(p-terma));
            //cargatotaaal += matriz_plasma[i].carga;
        //}
    }
    //fprintf(dat,"carga total: %lld\n",cargatotaaal);
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
float energia(void){
    int i,j;
    double energia_total=0, ecoulomb_total=0, aenergia_total=0, energia_i=0, d;
    for(i=1;i<=nceldas;i++){
        energia_i = 0;
        for(j=1;j<=nceldas;j++){
            if((matriz_plasma[i].carga!=0)&&(matriz_plasma[j].carga!=0)){
                if(i!=j){
                    //d = distanciacelda(i,j);
                    d = distancia(i,j);
                    energia_i += (qe*qe*matriz_plasma[i].carga*matriz_plasma[j].carga)/(4*pi*epce*epsi*d*esc);
                }
            }
        }
        aenergia_total += autoenergia(i);
        ecoulomb_total += energia_i;
    }
    energia_total = aenergia_total + 0.5*ecoulomb_total;
    if(p>=pasoinicial){
        dat = fopen("datos/energia.dat","a");
        fprintf(dat,"\n%i\t%e",p,energia_total);
        fclose(dat);
    }
    if(p%actu==0){
        dat = fopen("pruebas/energias.dat","a");
        fprintf(dat,"%d\t%e\t%e\t%e\n",p,ecoulomb_total,aenergia_total,energia_total);
        fclose(dat);
    }
    /*if(energia_total<0){
        //printf("\ncelda: %d ec_t: %e ae_t: %e e_total: %e d_%d_%d: %f",i,ecoulomb_total,aenergia_total,energia_total,dummy_i,dummy_j,distancia(dummy_i,dummy_j));
        imprimir_celda_plasma_rec(dummy_i);
        imprimir_celda_plasma_rec(dummy_j);
        printf("\nae_%d: %e ae_%d: %e ec: %e",dummy_i,autoenergia(dummy_i),dummy_j,autoenergia(dummy_j),(qe*qe*matriz_plasma[dummy_i].carga*matriz_plasma[dummy_j].carga)/(4*pi*epce*epsi*distancia(dummy_i,dummy_j)*esc));
        getchar();
    }*/
    return(energia_total);
}
////////////////////////////////////////////////////////////////////////////////////
void beta(void){
    int i,j,jm;
    long long int npi,nei,nii;
    float vol_i,beta_i;
    dat = fopen("datos/beta.dat","w");
    fprintf(dat,"#X\tY\n");
    for(i=1;i<=int(R*reso);i++){
        npi = nei = nii = 0;
        if(i==1){
            jm = 1;
        }
        else{
            jm = (int)((pi*i)*0.25);
        }
        vol_i = (pi*(pow(i,2)-pow(i-1,2))*esc*esc*1e-6)/8;
        for(j=1;j<=jm;j++){
            npi += matriz_plasma[ nmatriz_plasma[i][j] ].part[0] + matriz_plasma[ nmatriz_plasma[i][j] ].part[2] +
            matriz_plasma[ nmatriz_plasma[i][j] ].part[1] + matriz_plasma[ nmatriz_plasma[i][j] ].part[3];
            nei += matriz_plasma[ nmatriz_plasma[i][j] ].part[0];
            nii += matriz_plasma[ nmatriz_plasma[i][j] ].part[2] + matriz_plasma[ nmatriz_plasma[i][j] ].part[1];
        }
        beta_i = (2*muce*kb*8*(nei*tempee+nii*tempei))/(
        vol_i*pow(arreglo[186].absb[i]*1e-4,2));
        if(beta_i<0){
            printf("\nmuce: %e npi: %lld kb: %e tempee + tempeei: %e vol_i: %e |B_i|: %e",muce,npi,kb,tempee+tempei,vol_i,pow(arreglo[186].absb[i]*1e-4,1));
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
    //getchar();
    /*printf("\nm: %f b: %f y(3.5): %f",m, b, m*arreglo[186].r[36]+b);
    for(i=1;i<=24;i++){
        printf("r: %f |B|: %f\n",arreglo[186].r[36+i], arreglo[186].absb[36+i]);
    }
    getchar();*/
}
////////////////////////////////////////////////////////////////////////////////////
void gdr_plasma(void){
    for(int i=1;i<=nceldas;i++){
        for(int j=0;j<ncomp;j++)gpart[i][j] += matriz_plasma[i].part[j];
    }
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
void dinamica(void){
    float min_mass = masa[0], max_mass = masa[0], max_v;
    int min_charge = abs(carga[0]), max_charge = abs(carga[0]);
    max_v = part_plasma[1][0].v;
    for(int j=0; j<ncomp; j++){
        if(masa[j]<=min_mass)min_mass = masa[j];
        if(masa[j]>=max_mass)max_mass = masa[j];
        if(abs(carga[j])<=min_charge)min_charge = abs(carga[j]);
        if(abs(carga[j])>=max_charge)max_charge = abs(carga[j]);
        for(int i=1; i<=n_grupos[j]; i++){
            if(part_plasma[i][j].v>=max_v)max_v=part_plasma[i][j].v;
        }
    }
    //Zprintf("\nmasa[%d]: %e masa[%d]: %e masa[%d]: %e min_mass: %e",0,masa[0],1,masa[1],2,masa[2],min_mass);
    //printf("\ncarga[%d]: %d carga[%d]: %d carga[%d]: %d max_charge: %d",0,carga[0],1,carga[1],2,carga[2],max_charge);
    //printf("\nperiodo_min: %e dt_sincampo: %e dt(B): %e\n", ((2*pi*min_mass)/(100.0*max_charge*qe*B)), (esc/(100.0*max_v*reso)), B>1e-5?((2*pi*min_mass)/(100.0*max_charge*qe*B)):(esc/(100.0*max_v*reso)) );
    //getchar();
    //float dtau = B>1e-5?((2*pi*min_mass)/(100.0*max_charge*qe*B)):(esc/(100.0*max_v*reso));
    if(debug==3)printf("\nANTES CICLO\n");
    for(int j=0; j<ncomp-1; j++){
        //for(int i=0; i<ncomp; i++){
            //printf("\nnpart2_%d: %d",i,n_grupos[i]);
            //getchar();
        //}
        if(n_grupos[j]>0){
            for(int i=1; i<=n_grupos[j]; i++){

                if(debug==3&&i%100==0)printf("\rDENTRO");
                //printf("\nANTES DE DINAMICA i: %d j: %d\n",i,j);
                //if((10*i)%n_grupos[j]==0)printf("\ri: %d\n",i);
                //dinamica();
                float x,y,vx,vy,ax,ay;
                float mass = masa[j], dt;
                //dt = dtau/100.0;
                if(solodinamica){
                    dt = 0.01*( B>1e-5?((2*pi*min_mass)/(100.0*max_charge*qe*B)):(esc/(1.0*max_v*reso)) );
                }else{
                    dt = dt_termico/(100.0*(npart[0]+npart[1]));
                }
                int charge = carga[j];
                float Ex=0.0,Ey=0.0,Bz=B;
                float x0=part_plasma[i][j].x*esc, v0x=part_plasma[i][j].vx, y0 = part_plasma[i][j].y*esc, v0y = part_plasma[i][j].vy;
                float periodo = (2*pi*mass)/(abs(charge)*qe*Bz);
                int ktau = 100, contadork = 0;
                if(j>=ncomp-1)dt=0;
                //printf("\nq: %e m: %e q/m: %e",charge*qe,mass*m_u,charge*qe/mass/m_u);
                //printf("\nipart: %d p_x: %f p_y: %f x0: %f y0: %f",ipart,part_plasma[1][0].x,part_plasma[1][0].y,x0,y0);
                //printf("\ncarga: %d masa: %e charge: %d mass: %e",carga[0],masa[0],charge,mass);

                //dat=fopen("dinamica.dat","w");
                //fprintf(dat,"#t\tx\tvx\ty\tvy\n");
                //fprintf(dat,"%f\t%f\t%f\t%f\t%f\n",0.0,x0,v0x,y0,v0y);
                //for(int i=1;i<=int(periodo/dt)+1;i++){
                //long long int k;
                for(int k=1;k<=ktau;k++){
                    if(debug==3&&k%1==0){
                        contadork++;
                        printf("\rDENTRO CICLO K j: %d i: %d rho: %f\t\tcontador: %d",j,i,norma(x,y)/esc,contadork);
                    }
                    ax = ((charge*qe)/(mass))*(Ex+v0y*Bz);
                    ay = ((charge*qe)/(mass))*(Ey-v0x*Bz);
                    x = x0 + v0x*dt + 0.5*ax*dt*dt;
                    vx = v0x + ax*dt;
                    y = y0 + v0y*dt + 0.5*ay*dt*dt;
                    vy = v0y + ay*dt;
                    //printf("\ndx: %e dy: %e dt: %e",x-x0,y-y0,dt);
                    //getchar();
                    /*printf("ax: %f ay: %f\n%f\t%f\t%f\t%f\t%f\n",ax,ay,i*dt,x,vx,y,vy);
                    getchar();*/
                    /*if(i%itau==0){
                        fprintf(dat,"%e\t%f\t%f\t%f\t%f\n",i*dt,x,vx,y,vy);
                        //if(jcomp==0)printf("\rp: %d j: %d itau: %d i: %lld",p,jcomp, itau, i);
                    }*/
                    //printf("\nx0: %e y0: %e vx0: %e vy0: %e",x0,y0,v0x,v0y);
                    //printf("\nx: %e y: %e vx: %e vy: %e ax: %e ay: %e",x,y,vx,vy,ax,ay);
                    //getchar();
                    //if(norma(x,y)>R*esc){
                    //if(part_a_mallado(x/esc,y/esc)>nceldas||part_a_mallado(x/esc,y/esc)<1){
                    if( norma(x,y)>=(R-2)*esc ){
                        /*dat = fopen("vectores.dat","w");
                        dat2 = fopen("vectoresfinales.dat","w");
                        printf("\nANTES v0x: %f v0y: %f",v0x, v0y);
                        reflejar_pared(i,j,&v0x,&v0y);
                        printf("\nDESPUES v0x: %f v0y: %f",v0x, v0y);

                        fclose(dat);
                        fclose(dat2);
                        printf("\nQue show");getchar();*/
                        k-=1;
                        if((v0x*cos( atan2(y0,x0) )+v0y*sin( atan2(y0,x0) ))>=0){
                            reflejar_pared(i,j,&v0x,&v0y);
                        }else{
                            x0 = x; v0x = vx;
                            y0 = y; v0y = vy;
                        }
                        //intercambiar(i,n_grupos[j],j);
                        //n_grupos[j]--;
                    }else{
                        x0 = x; v0x = vx;
                        y0 = y; v0y = vy;
                    }
                }
                //dinamica(i,j,B>1e-5?((2*pi*min_mass)/(100.0*max_charge*qe*B)):(esc/(100.0*max_v*reso)) );
                //dinamica(i,j,1e-7 );
                part_plasma[i][j].x = x/esc;
                part_plasma[i][j].y = y/esc;
                part_plasma[i][j].vx = vx;
                part_plasma[i][j].vy = vy;
            }
            for(int i=1; i<=n_grupos[j]; i++){
                //printf("\nANTES DE DINAMICA i: %d j: %d\n",i,j);
                //if(norma(part_plasma[i][j].x,part_plasma[i][j].y)>R){
                if(part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y)>nceldas||part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y)<1){
                    printf("\ni: %d n_grupos[%d]: %d rho2: %f",i,j,n_grupos[j],norma(part_plasma[i][j].x,part_plasma[i][j].y));
                    getchar();
                }
            }
        }else{
            printf("\nn_grupos[%d]: %d",j,n_grupos[j]);
            getchar();
        }
    }
    if(debug==2)printf("\nDESPUES CICLO");
    //fprintf(dat,"radio:\t%e",mass*norma(v0x,v0y)/(fabs(charge)*qe*Bz));
    //fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void intercambiar(int a, int b, int c){
    part_plasma[0][c] = part_plasma[a][c];
    part_plasma[a][c] = part_plasma[b][c];
    part_plasma[b][c] = part_plasma[0][c];
}
////////////////////////////////////////////////////////////////////////////////////
float distribucion_normal_5(float mu, float sigma){//usando 5 sigmas de intervalo
    float b, fb;
    do{
        b=alea_f(-5*sigma+mu,5*sigma+mu);//sqrt((kb*tempee)/me), 4*sqrt((kb*tempee)/me) );
        fb=exp(-(b-mu)*(b-mu)/(2*sigma*sigma))/(sqrt(2*pi)*sigma);
    }while(alea()/(sqrt(2*pi)*sigma)>fb);
    return(b);
}
////////////////////////////////////////////////////////////////////////////////////
void prueba_dist(void){
    int i,j;
    long int clases[201] = {0};
    float x, mu = 0, sigma = 0.005;

    /*for(i=1;i<=1000000;i++){
        clases[ (int)(200*(distribucion_normal_5(sigma)+4*sigma)/(8*sigma))+1 ]++;
    }*/
    for(i=1;i<=10000000;i++){
        clases[ int( 200*(distribucion_normal_5(mu,sigma)+5*sigma-mu)/(10*sigma) )+1]++;
    }

    dat = fopen("aaahistograma.dat","w");
    for(i=1;i<=200;i++){
        fprintf(dat,"%f\t%d\n",-5*sigma+(i-1)*10*sigma/200.0+mu,clases[i]);
    }
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
float norma(float a, float b){
    return(sqrt(a*a+b*b));
}
////////////////////////////////////////////////////////////////////////////////////
void hacer_rms(int a){
    double ms_x = 0, ms_y = 0, rms_x = 0, rms_y = 0;
    //CALCULAR EL MS EN X y Y
    for(int i=1;i<=nceldas;i++){
        ms_x += matriz_plasma[i].x*matriz_plasma[i].x*matriz_plasma[i].part[0];
        ms_y += matriz_plasma[i].y*matriz_plasma[i].y*matriz_plasma[i].part[0];
        /*printf("\nx: %e\ty: %e\telectrones: %lld\tprodx: %e\tprody: %e",matriz_plasma[i].x,matriz_plasma[i].y,matriz_plasma[i].part[0],matriz_plasma[i].x*matriz_plasma[i].x*matriz_plasma[i].part[0],matriz_plasma[i].y*matriz_plasma[i].y*matriz_plasma[i].part[0]);
        getchar();*/
    }
    rms_x = sqrt(1.0*ms_x/npart[0]);
    rms_y = sqrt(1.0*ms_x/npart[0]);
    if(a<=pasoinicial){
        dat = fopen("rms.dat","w");
    }else{
        dat = fopen("rms.dat","a");
        fprintf(dat,"%d\t%e\t%e\n",a,rms_x,rms_y);
    }
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void crear_archivos_iniciales_en_blanco(void){
    system("mkdir datos");
    system("mkdir datos/promedios");
    system("mkdir datos/PNG");
    system("mkdir histogramas");
    system("mkdir pruebas");

    dat = fopen("pruebas/naleprom.dat","w");
    fclose(dat);
    dat = fopen("pruebas/nale_f.dat","w");
    fclose(dat);
    dat = fopen("pruebas/energias.dat","w");
    fclose(dat);
    dat = fopen("pruebas/denergias_acep.dat","w");
    fclose(dat);
    dat = fopen("pruebas/denergias_rechazo.dat","w");
    fclose(dat);
    dat = fopen("pruebas/razon_acep.dat","w");
    fclose(dat);
    dat = fopen("pruebas/rms_part.dat","w");
    fclose(dat);
    dat=fopen("datos/cargatotal.dat","w");
    fclose(dat);


    for(int i=1;i<=R_reso*R_reso+1;i++){
        //matriz_plasma[i].part[0]=matriz_plasma[i].part[3]=matriz_plasma[i].part[2]=matriz_plasma[i].part[1]=matriz_plasma[i].carga = 0;
        matriz_plasma[i].carga=matriz_plasma[i].part[0]=matriz_plasma[i].part[3]=matriz_plasma[i].part[2]=matriz_plasma[i].part[1]=0;
        matriz_plasma[i].phi=matriz_plasma[i].rho=matriz_plasma[i].ve=matriz_plasma[i].vh2p=matriz_plasma[i].vhp=matriz_plasma[i].vx[0]=matriz_plasma[i].vxh20=matriz_plasma[i].vxh2p=matriz_plasma[i].vxhp=matriz_plasma[i].vye=matriz_plasma[i].vyh20=matriz_plasma[i].vyh2p=matriz_plasma[i].vyhp=matriz_plasma[i].x=matriz_plasma[i].y=0;
        for(int j=0;j<=3;j++){
            matriz_plasma[i].t[j]=0;
            matriz_plasma[i].v[j]=0;
            matriz_plasma[i].vy[j]=0;
            matriz_plasma[i].vx[j]=0;
            matriz_plasma[i].part[j]=0;
        }
    }
    for(int i=0; i<part_plasma_size; i++){
        for(int j=0; j<4; j++){
            part_plasma[i][j].contador = 0;
            part_plasma[i][j].part = 0;
            part_plasma[i][j].v = part_plasma[i][j].vx = part_plasma[i][j].vy = part_plasma[i][j].x = part_plasma[i][j].y = 0;
        }
    }
    for(int i=1; i<=nceldas; i++)for(int j=0;j<ncomp;j++)gpart[i][j]=0;
}
////////////////////////////////////////////////////////////////////////////////////
void cambiar_a(int* a, float* b){
    *a+=1;
    *b = *b + 12.5;
}
////////////////////////////////////////////////////////////////////////////////////
