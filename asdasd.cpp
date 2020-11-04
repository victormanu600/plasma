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
int part_a_mallado(float a, float b, int c, int d, int e, int f, int g);
void intercambiar(int a, int b, int c );
void crear_matriz_plasma(void);
void vecinos(int a);
void calc_carga(int a);
void calc_contador(void);
void calc_dens(int a);
void hacer_histograma(int a, int b);
void hacer_distribucion(int a);
void hacer_rms(int a);
void crear_archivos_iniciales_en_blanco(void);

float distancia(int a, int b);
int signo(float a);
float norma(float a, float b);
void mover_particulas_part(void);
void pared(void);
//void reflejar_pared(int a, int b, float* vx, float* vy);
void reflejar_pared(float x0, float y0, float *vx, float *vy);
void metropolis_plasma(void);
void dinamica(void);
void cambiar_a(int* a, float* b);

double de_plasma(void);
double de_plasma_lineal(void);
float calcular_de_mov(void);
float autoenergia(int a);
void calcular_momento_total(void);
void gdr_plasma(void);
void actu_salida(void);
void actu_salida_carga(void);
void actu_salida_dina(int a);
void salida_prom(void);
double energia(void);
void beta(void);

void prueba_dinamica(void);
float distribucion_normal_5(float mu, float sigma);
void prueba_dist(void);


//////////////////////////////Constantes
int pasos, actu, terma, pasoinicial=1, termaanterior, iprint[6], ncomp, separacion;
bool solodinamica, b_dinamica, blineal=true;
int contadorrr = 0;
float esc, densidad[4], volumen, nelectronesr, tempee, tempei, tempe[4], B, dx, dy;
float reso;
long long int nh20, nhp, nh2p, nelectrones, npart[4], naleprom=0, part_por_grupo[4] ;
float R;
const int R_reso = (2*(110))*11/11+1, part_plasma_size = 6000001;                      //Primer numero es R, segundo es reso, 2 es por que es de -R a R y +1 para comenzar los arreglos desde 1
float me, mhp, mh2p, mh20, masa[4];
int carga[4];
int n_grupos[4];
///////////////////////////////////////////////////////////////////////Variables globales
int p=0, tipo, ienergia=0;
int rechazo;
int nceldas, n1, n2, n1i, n2i, alea_part;                                       //numero de celda Nueva/Vieja Inicial/Final
char cero[] = "         0";
int rasdnd;
double dem, de_coulomb = 0, de_autoenergia = 0/*, energia_p = 0*/, energia_con_de = 0;
double dt_termico = 0;//, energia_promedio = 0, desviacion2_energia = 0, desviacion = 0, energia_prom[10+1] = {0};
long long int part_desplazadas[4] = {0};
int debug=0;
float k_auto = 0, dtau = 0, de_frac=0;
bool win_os;
float min_mass, max_mass;
int min_charge, max_charge;
float x00[4], rmin[2]={0},rmax[2]={0};
time_t inicio, final;
//int celdas_mc1[8]={0}, celdas_mc2[8]={0};
float dxx, dyy;
///////////////////////////////////////////////////////////////////////Contadores
int c_mov=0,c_mova=1, c_dina=0;
int c_uno=0,c_unoa=0,c_dos=0,c_dosa=0,c_tres=0,c_tresa=0;
int contador[4],contador_a[4], contador_rechazo=0, contador_nale=0,contador_nalef=0;
int rechazo_neg=0, rechazo_met_mov=0;

long long int gelectron[R_reso*R_reso+1+2], gh20[R_reso*R_reso+1+2], gh2p[R_reso*R_reso+1+2], ghp[R_reso*R_reso+1+2], gpart[R_reso*R_reso+1+2][4];

struct smatriz_plasma{
	float rho,phi,x,y;
    float vx[4], vy[4], v[4];
    float anchow, anchol;
	//float t[3];
	long long int part[4],carga;
	int xi, yi;
	int nvecinos, vecinos[10], mov_a;
}matriz_plasma[R_reso*R_reso+1+10];

struct sceldas_mc{
    int size;
    int iniciales[8], finales[8];
}celdas_mc;

/*struct smatriz_plasma2{
	float rho,phi,x,y;
    float vxe, vye, vxh2p, vyh2p, vxhp, vyhp, vxh20, vyh20;
    float vx[4], vy[4];//, vx_1[4], vy_1[4], vx_2[4], vy_2[4];
    float ve, vhp, vh2p, v[4];
    float anchow, anchol;
	float t[3];
	long long int h20,hp,h2p,electrones,part[4],carga;
	int xi, yi;
	int nvecinos, vecinos[10];
}matriz_plasma2[R_reso*R_reso+1+2];*/
//}matriz_plasma2[2000*2000+1];

struct spart_plasma{
    float x, y, vx, vy, v;
    long long int part;
    int contador, celda_antes;
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
float superficie(float a, float b){
    return(exp(-(a*a)/5));
}

int main(){
    //prueba_dinamica();
    srand((unsigned)time(NULL));
/*
        dat = fopen("asd.dat","w");
        for(int asd=0; asd<1000; asd++){
            fprintf(dat,"%d %f %f %f\n",asd,asd/1000.0,asd/2000.0,asd/3000.0);
        }
        fclose(dat);
        pasoinicial=200000000;
        int paso_asd = 0;
        float float_asd[3]={0};
            //fprintf(dat,"%d",p);
            //for(int i=0; i<ncomp; i++)fprintf(dat,"\t%1.5f",contador_a[i]/(contador[i]*1.0));
            //fprintf(dat,"\n");
        if( dat = fopen("razon_acep_inicial.dat","r") ){
            while( (fscanf(dat,"%d\t%f\t%f\t-nan\n",&paso_asd,&float_asd[0],&float_asd[1])!=EOF) ){
                printf("%d %1.5f %1.5f %1.5f\n",paso_asd,float_asd[0],float_asd[1],float_asd[2]);
            }
        }
        printf("\nDESPUES CICLO V:");
        getchar();

    return(0);*/
    /*int xreso = 50, yreso = 50;
    dat = fopen("datos.dat","w");
    for(int i=-5*xreso;i<=5*xreso;i++){
        for(int j=-5*yreso;j<=5*yreso;j++){
            float x_i=1.0*i/xreso, y_i=1.0*j/yreso;
            fprintf(dat,"%f\t%f\t%f\n",x_i,y_i,superficie(x_i,y_i) );
        }
        fprintf(dat,"\n");
    }
    fclose(dat);*/
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

	if(debug>=1)cout << endl << "Condiciones iniciales" << endl;
	condiciones_iniciales();
	if(debug>=1)cout << endl << "Crear matriz" << endl;
	crear_matriz_plasma();
	if(debug>=1)cout << endl << "Arreglo inicial" << endl;
    arreglo_inicial();
    energia_con_de = energia();
    cout << "Energia: " << energia_con_de << endl;
	if(debug>=1)cout << endl << "Despues arreglo inicial" << endl;
	//actu_salida_carga();
    //return(0);

    time(&inicio);
    for(p=pasoinicial;p<=pasos;p++){
        if(debug>=2)cout << "\nINICIO CICLO\n";
        b_dinamica=false;
        if(solodinamica)b_dinamica=true;
        //if(p>=(pasos-1000))b_dinamica=true;
        if((dt_termico/(npart[0]+npart[1]))>=dtau)b_dinamica=true;
        //if(p==pasoinicial)b_dinamica=true;
        //if((desviacion/energia_promedio)<de_frac&&bool_energia_prom)b_dinamica=true;
        if(!b_dinamica){
            rechazo = 0; c_mov++;

            if(debug>=2)cout << "\nANTES MOVER PART\n";
            mover_particulas_part();
            part_plasma[n1][tipo].contador++;

            if(debug>=2)cout << "\nANTES PARED\n";
            pared();

            if(debug>=2)cout << "\nANTES METROPOLIS\n";
            if(rechazo == 0){
                metropolis_plasma();
            }else{

                part_plasma[alea_part][tipo] = part_plasma[0][tipo];
                if(blineal){
                    for(int i=0; i<8; i++)matriz_plasma[celdas_mc.finales[i]]=matriz_plasma[celdas_mc.finales[i]];
                }else{
                    matriz_plasma[n1] = matriz_plasma[n1i];
                    matriz_plasma[n2] = matriz_plasma[n2i];
                }

                //importar_part_a_mallado();
                rechazo_neg++;

                //printf("\nQue show v:2");
                //getchar();
                //importar_part_a_mallado();
            }
        }else{
            c_dina++;
            if(debug>=2)cout << "\nANTES DINAMICA\n";
            //if(!solodinamica)actu_salida_dina(0);
            dinamica();
            if(debug>=2)cout << "\nANTES IMPORTAR PART A MALLADO\n";
            importar_part_a_mallado();
            //if(!solodinamica)actu_salida_dina(1);
        }
        if((b_dinamica||p%actu==0)&&!solodinamica){
            dat = fopen("pruebas/energia_dinamica.dat","a");
            fprintf(dat, "%d %e",p,energia_con_de);
            fclose(dat);
            energia_con_de = energia();
            dat = fopen("pruebas/energia_dinamica.dat","a");
            fprintf(dat, " %e\n",energia_con_de);
            fclose(dat);
        }
        if(p>terma){
            if(debug>=1)cout << "\nANTES GDR\n";
            gdr_plasma();
            if(p%actu==0){
                if(debug>=1)cout << "\nANTES ACTUSALIDA2\n";
                actu_salida();
                if(debug>=1)cout << "\nANTES SALIDAPROM\n";
                salida_prom();
                for(int i=0; i<ncomp; i++)hacer_histograma(i,p);
                calc_contador();
                calc_dens(1);
            }
        }
        if(p%1000==0||b_dinamica){
            if(win_os){
                printf("\rPaso: %9d ",p);
                for(int i=0; i<ncomp; i++)printf("Acep%d: %1.5f ",i,1.0*contador_a[i]/contador[i]);
                printf("dt/npart: %e ",(dt_termico/(npart[0]+npart[1])));
                printf("Rneg: %d uno: %1.5f ",rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0));
                printf("Dinamica: %d ",c_dina);
                if(p%actu==0)printf("Energia: %e ",energia_con_de);
            }
            dat = fopen("pruebas/razon_acep.dat","a");
            fprintf(dat,"%d",p);
            for(int i=0; i<ncomp; i++)fprintf(dat,"\t%1.5f",contador_a[i]/(contador[i]*1.0));
            fprintf(dat,"\n");
            fclose(dat);
            //contador_nale=0;naleprom=0,nale_f=0;
        }

        /*long long int cargatotaaal = 0, cargatotalpart = 0, min_carga = matriz_plasma[1].carga, max_carga = matriz_plasma[1].carga, maxx, carga_cero=0;
        if(p%100000==0){
            for(int i=1;i<=nceldas;i++){
                if(fabs( matriz_plasma[i].x)>(R+2)||fabs( matriz_plasma[i].y)>(R+2)){
                    printf("\ni: %d x: %e y: %e rho: %e",i,matriz_plasma[i].x,matriz_plasma[i].y,matriz_plasma[i].rho);
                    getchar();
                }
                if(matriz_plasma[i].rho<=(R)){
                    cargatotaaal += matriz_plasma[i].carga;
                    if(matriz_plasma[i].carga>max_carga)max_carga=matriz_plasma[i].carga;
                    if(matriz_plasma[i].carga<min_carga)min_carga=matriz_plasma[i].carga;
                }
            }
            for(int j=0; j<ncomp; j++){
                for(int i=1; i<=n_grupos[j]; i++){
                    cargatotalpart+=part_plasma[i][j].part*carga[j];
                }
            }
            if( abs(min_carga)>abs(max_carga) ){
                maxx = abs(min_carga);
            }else{
                maxx = abs(max_carga);
            }
            for(int i=1;i<=nceldas;i++){
                if( ( abs( matriz_plasma[i].carga ) < 0.01*maxx )&&( matriz_plasma[i].rho<=(R-12) ) )carga_cero++;
            }
            dat = fopen("pruebas/carga_cero.dat","a");
            fprintf(dat,"%d %lld\n",p,carga_cero);
            fclose(dat);
            dat = fopen("pruebas/carga.dat","a");
            fprintf(dat,"%d %lld %lld\n",p,cargatotaaal,cargatotalpart);
            fclose(dat);
        }
        if(cargatotaaal!=0){
            printf("\ncargatotal: %lld alea_part: %d tipo: %d celda: %d",cargatotaaal,alea_part,tipo,part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y,2,alea_part,tipo,0,0));
            imprimir_celda_plasma_part(alea_part,tipo);
            printf("\nPaso: %d ",p);
            for(int i=0; i<ncomp; i++)printf("Acep%d: %1.5f ",i,1.0*contador_a[i]/contador[i]);
            printf("dt/npart: %e ",(dt_termico/(npart[0]+npart[1])));
            printf("Rneg: %d uno: %1.5f ",rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0));
            printf("Dinamica: %d",c_dina);
            getchar();
        }*/
        if((p%1000==0)||(solodinamica&&(p%10==0))){
            time(&final);
            dat = fopen("pruebas/tiempo.dat","a");
            fprintf(dat,"%d %f\n",p,difftime(final,inicio));
            fclose(dat);
            //getchar();
            time(&inicio);
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
///////////////////////////////////////////////////////////////////////////////////////////////
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
        fscanf(dat,"Densidades (mt-3):");
        for(int i=0; i<ncomp; i++){
            if(i<ncomp-1){
                fscanf(dat," %f,",&densidad[i]);
            }else{
                fscanf(dat," %f\n",&densidad[i]);
            }
        }
        fscanf(dat,"Campo magnetico (T): %f\n",&B);
        fscanf(dat,"Constante AE: %f\n",&k_auto);
        fscanf(dat,"dx (mm): %f\n",&dx);
        fscanf(dat,"dy (mm): %f\n",&dy);
        fscanf(dat,"fraccion desv: %f\n",&de_frac);
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
        fscanf(dat,"Separacion (mm): %d\n",&separacion);
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
    //getchar();
}
///////////////////////////////////////////////////////////////////////////////////////////////
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
    printf("Densidades (mt-3):");
    for(int i=0; i<ncomp; i++){
        if(i<ncomp-1){
            printf(" %e,",densidad[i]);
        }else{
            printf(" %e\n",densidad[i]);
        }
    }
    printf("Campo magnetico (T): %e\n",B);
    printf("Constante AE: %f\n",k_auto);
    printf("dx (mm): %e\n",dx);
    printf("dy (mm): %e\n",dy);
    printf("fraccion desv: %f\n",de_frac);
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
    printf("Solo dinamica (S_1/N_0): %d\n",solodinamica);
    printf("Separacion (mm): %d\n\n",separacion);
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_iniciales(){
    leer_datos_iniciales();
    esc = 1e-3;
    R -= 2.0;
    //densidad[0] = 2e19;
    volumen = pi*R*esc*R*esc*1e-6;
    //nelectronesr = densidad*volumen;
    for(int i=0; i<ncomp; i++)npart[i] = densidad[i]*volumen;
    /*npart[1] = densidad[1]*volumen;
    npart[2] = densidad[2]*volumen;*/

    //npart[1] = npart[0];
    //npart[2] = 9*npart[0];
    /*nelectrones = npart[0] = nelectronesr;
    nhp = npart[1] = nelectrones;
    nh2p = npart[2] = 9*nelectrones;//0;//0.2*nelectrones+1;
    nh20 = npart[3] = 9*nelectrones;*/
    printf("\nvolumen: %e\n",volumen);
    for(int i=0; i<ncomp; i++)printf("densidad_%d: %e ",i,densidad[i]); cout << endl;
    for(int i=0; i<ncomp; i++)printf("npart_%d: %lld ",i,npart[i]); cout << endl;
    //getchar();getchar();

    //tempe[0] = tempee; tempe[1] = tempei; tempe[2] = tempei; tempe[3] = tempei;
    for(int i=0; i<ncomp; i++)tempe[i] *= qe/kb; //Conversion T[eV]->T[K]
    tempei = tempei*qe/kb;
    tempee = tempee*qe/kb;

    for(int i=0; i<ncomp; i++)printf("T_%d: %f ",i,tempe[i]); cout << endl;

    printf("\ndistancia: %f\n",1.0*esc/reso);
    for(int i=0; i<ncomp; i++)printf("vprom_%d: %e dt_%d: %e\n",i,sqrt((8*kb*tempe[i])/(pi*masa[i])),i,1.0*esc/(reso*sqrt((8*kb*tempe[i])/(pi*masa[i]))));cout << endl;
    //getchar();

    //printf("\ntempei: %f tempee: %f\n", tempei, tempee);
    //getchar();

    for(int i=0; i<ncomp; i++)contador[i] = 0;
    for(int i=0; i<ncomp; i++)contador_a[i] = 0;

    min_mass = masa[0]; max_mass = masa[0];
    min_charge = abs(carga[0]); max_charge = abs(carga[0]);
    for(int j=0; j<ncomp; j++){
        if(masa[j]<=min_mass)min_mass = masa[j];
        if(masa[j]>=max_mass)max_mass = masa[j];
        if(abs(carga[j])<=min_charge)min_charge = abs(carga[j]);
        if(abs(carga[j])>=max_charge)max_charge = abs(carga[j]);
    }

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
            if( (i*i+j*j<=int(R*reso)*int(R*reso))||contadorm2>0 ){         //SE ACEPTAN TODOS DE MOMENTO
                contadorm++;
                nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1]=contadorm;
                matriz_plasma[contadorm].xi = i;
                matriz_plasma[contadorm].yi = j;
                matriz_plasma[contadorm].x = (1.0*(i))/reso;
                matriz_plasma[contadorm].y = (1.0*(j))/reso;
                matriz_plasma[contadorm].rho = sqrt(i*i+j*j)/reso;
                matriz_plasma[contadorm].phi = (atan2(matriz_plasma[contadorm].y,matriz_plasma[contadorm].x)+((atan2(matriz_plasma[contadorm].y,matriz_plasma[contadorm].x))<0?2*pi:0.0));
                matriz_plasma[contadorm].anchow = 1.0*esc/reso;
                matriz_plasma[contadorm].anchol = 1e-6;

                if(fabs( matriz_plasma[contadorm].x)>(R+1)||fabs( matriz_plasma[contadorm].y)>(R+1)){
                    dat = fopen("pruebas/getchares.dat","a");
                    fprintf(dat,"Dentro de crear_matriz_plasma, dentro de if en ciclo. paso: %d\n",p);
                    fprintf(dat,"\ni: %d x: %e y: %e rho: %e",contadorm,matriz_plasma[contadorm].x,matriz_plasma[contadorm].y,matriz_plasma[contadorm].rho);
                    fclose(dat);
                    getchar();
                }
            }else{
                contadorm2++;
                nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1]=contadorm2;
                matriz_plasma[contadorm2].xi = i;
                matriz_plasma[contadorm2].yi = j;
                matriz_plasma[contadorm2].x = (1.0*(i+0.1))/reso;
                matriz_plasma[contadorm2].y = (1.0*(j+0.1))/reso;
                matriz_plasma[contadorm2].rho = sqrt(i*i+j*j)/reso;
                matriz_plasma[contadorm2].phi = (atan2(matriz_plasma[contadorm2].y,matriz_plasma[contadorm2].x)+((atan2(matriz_plasma[contadorm2].y,matriz_plasma[contadorm2].x))<0?2*pi:0.0));
                matriz_plasma[contadorm2].anchow = 1.0*esc/reso;
                matriz_plasma[contadorm2].anchol = 1e-6;
            }
        }
    }
    nceldas=contadorm;

    n1i = nceldas+1;
    n2i = nceldas+2;
    celdas_mc.iniciales[0]=nceldas+3;
    celdas_mc.iniciales[1]=nceldas+4;
    celdas_mc.iniciales[2]=nceldas+5;
    celdas_mc.iniciales[3]=nceldas+6;
    celdas_mc.iniciales[4]=nceldas+7;
    celdas_mc.iniciales[5]=nceldas+8;
    celdas_mc.iniciales[6]=nceldas+9;
    celdas_mc.iniciales[7]=nceldas+10;
    //n_grupos[0] = n_grupos[1] = n_grupos[2] = 1100*nceldas;
    //n_grupos[0] = n_grupos[1] = n_grupos[2] = 140*nceldas;
    //n_grupos[0] = n_grupos[1] = n_grupos[2] = 200*nceldas;
    n_grupos[0] = n_grupos[1] = n_grupos[2] = 60000;
    //n_grupos[0] = n_grupos[1] = n_grupos[2] = 2;
    printf("\nnceldas: %i\tn1i: %i\tn2i: %i R_reso: %d\n",nceldas,n1i,n2i,R_reso);
    for(j=0;j<ncomp;j++)printf("n_grupos_%d: %d ",j,n_grupos[j]); printf("\n");
    part_por_grupo[0]=npart[0]/n_grupos[0];
    part_por_grupo[1]=npart[1]/n_grupos[1];
    part_por_grupo[2]=npart[2]/n_grupos[2];
    for(j=0;j<ncomp;j++)printf("Part_por_grupo_%d: %lld ",j,part_por_grupo[j]); printf("\n");
    for(i=1;i<=nceldas;i++){
        vecinos(i);
    }

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
    dat = fopen("pruebas/vecinos.dat","w");
    for(i=1; i<=nceldas; i++){
        fprintf(dat,"%d\t%d\t",i,matriz_plasma[i].nvecinos);
        //ovecinos << i << "\t" << matriz_plasma[i].nvecinos << "\t";
        for(j=1; j<=matriz_plasma[i].nvecinos; j++){
            fprintf(dat,"%d\t",matriz_plasma[i].vecinos[j]);
            //ovecinos << matriz_plasma[i].vecinos[j] << "\t";
        }
        fprintf(dat,"\n");
        //ovecinos << endl;
    }
    fclose(dat);
    //ovecinos.close();
    printf("\nACABO VECINOS\n");
    //getchar(); getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void vecinos(int a){
    int cont_vecinos = 0, ix = floor(matriz_plasma[a].x*reso)+int(R*reso)+1, iy = floor(matriz_plasma[a].y*reso)+int(R*reso)+1, iobj;
    //float xx, yy;
    /*cout << "ix: " << ix << " iy: " << iy << endl;
    getchar();*/
    for(int i=-1;i<=1;i++){
        for(int j=-1;j<=1;j++){
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
    if(cont_vecinos>=9){
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
void calc_contador(void){
    sprintf(salidac,"datos/contador_%d.dat",p);
    int contador_mov=0;
    dat = fopen(salidac,"w");
    for(int i=1; i<=nceldas; i++){
        fprintf(dat,"%f\t%f\t%d\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].mov_a);
        contador_mov+=matriz_plasma[i].mov_a;
    }
    fprintf(dat,"contador_mov: %d\t suma_mov: %d",contador_mov,contador_a[0]+contador_a[1]+contador_a[2]);
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void calc_dens(int a){
    int clasesr = 200;
    float rmin=0, rmax=R;
    double anchor;
    long long int I_total[5]={0}, gdr[5001][5];
    for(size_t k1 = 0; k1 < 5001; k1++){
        for(size_t k2 = 0; k2 < 5; k2++){
            gdr[k1][k2]=0;
        }
    }
    double r_i;
    anchor = (rmax - rmin)/clasesr;

    for(int j=0; j<ncomp; j++){
        for( int i = 1; i <= n_grupos[j]; i++ ) {
            r_i = norma( part_plasma[i][j].x, part_plasma[i][j].y );
            gdr[(int)((r_i-rmin)/(anchor))+1][j]+=part_plasma[i][j].part;
            I_total[ j ] += part_plasma[i][j].part;
            //printf("\n%f\t%lld\t%lld",)
        }
    }
    char densi[50];
    if(a==0){
        sprintf(densi,"datos/densidad_inicial.dat");
    }else{
        sprintf(densi,"datos/densidad_%d.dat",p);
    }
    dat = fopen(densi,"w");
    for(int i=1;i<=clasesr;i++){
        fprintf(dat,"%f",rmin+(i-0.5)*anchor);
        for(int j=0; j<ncomp; j++)fprintf(dat,"\t%e",1.0*gdr[i][ j ]/I_total[j] );
        fprintf(dat,"\n");
        /*printf("%f",rmin+i*anchor);
        for(int k=0; k<ncomp; k++)printf("\t%e",1.0*gdr[i][ k ]/I_total[k] );
        //printf("\n");
        getchar();*/
    }
    fclose(dat);
    //printf("\nQUE SHOW TERMINO DENS");
}
////////////////////////////////////////////////////////////////////////////////////
void arreglo_inicial(void){
    float xrand, yrand, theta, r, dummy, R_inicial = R-12;
	int max_n_grupos = n_grupos[0];
    int part_cont = 0;
    int dinaanterior;

    printf("\rAsignando arreglo inicial\n");
    if( dat = fopen("posiciones_inicial.dat","r") ){
        fscanf(dat,"%d %d %d\n",&pasoinicial,&termaanterior,&ncomp);
        fclose(dat);
        printf("\npasoinicial: %d\tterma_anterior: %d\tncomp: %d\n",pasoinicial,termaanterior,ncomp);
    }
    if(dat = fopen("energia_inicial.dat","r")){
        dat2 = fopen("datos/energia.dat","w");
        while( (ienergia<pasoinicial)&&(fscanf(dat,"\n%i\t%f",&ienergia,&dummy)!=EOF) ){
            fprintf(dat2,"\n%i\t%e",ienergia,dummy);
        }
        fclose(dat);fclose(dat2);
        printf("Dentro del segundo if\n");
    }
    //getchar();
    if((pasoinicial>0)&&(ienergia==pasoinicial)){
        printf("\nContinuando posicion final de corrida anterior");
        //getchar();
        if( dat = fopen("posiciones_inicial.dat","r") ){

            fscanf(dat,"%d %d %d ",&pasoinicial,&termaanterior,&ncomp);
            for(int i=0;i<ncomp;i++)fscanf(dat,"%f ",&dummy );
            fscanf(dat,"%d",&dinaanterior);

            //fscanf(dat,"%d %d %d %f %f %f %f %d %d",&pasoinicial,&termaanterior,&ncomp,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
            /*if(termaanterior>pasoinicial)terma=termaanterior-pasoinicial;
            else terma = termaanterior;*/
            terma = termaanterior;
            printf("\nterma: %d\tp: %d\n",terma,p);
            for(int i=0; i<ncomp; i++ )if(n_grupos[i]>=max_n_grupos)max_n_grupos = n_grupos[i];
            for(int i=1;i<=max_n_grupos;i++){
                if((10*i)%(max_n_grupos)==0)printf("\ri: %d",i);
                fscanf(dat, "\n%d ",&dummy);
                for(int j=0; j<ncomp; j++){
                    if(i<=n_grupos[j])fscanf(dat,"%f %f %f %f %lld ",&part_plasma[i][j].x, &part_plasma[i][j].y, &part_plasma[i][j].vx, &part_plasma[i][j].vy, &part_plasma[i][j].part);
                }
            }
            /*for(int i=1;i<=nceldas;i++){
                fscanf(dat,"\n%i\t%lld\t%lld\t%lld\t%lld\t%lld",&dummy,&matriz_plasma[i].part[0],&matriz_plasma[i].part[3],&matriz_plasma[i].part[1],&matriz_plasma[i].part[2],&matriz_plasma[i].carga);
                if((matriz_plasma[i].part[0]==0)||(matriz_plasma[i].part[3]==0)||(matriz_plasma[i].part[1]==0)||(matriz_plasma[i].part[2]==0)||(matriz_plasma[i].carga==0)){
                    //imprimir_celda_plasma(i);
                }
            }*/
            fclose(dat);
            importar_part_a_mallado();
            //printf("\nenergia inicial: %e",energia());
            //getchar();
            pasos+=pasoinicial;
        }
        int paso_racep=0;
        float float_part[3]={0};
        if( dat = fopen("razon_acep_inicial.dat","r") ){
            dat2 = fopen("pruebas/razon_acep.dat","w");
            while( (paso_racep<pasoinicial)&&(fscanf(dat,"%d\t%f\t%f\t-nan\n",&paso_racep,&float_part[0],&float_part[1])!=EOF) ){
                fprintf(dat2,"%d\t%1.5f\t%1.5f\t-nan\n",paso_racep,float_part[0],float_part[1]);
            }
            fclose(dat);
            fclose(dat2);
        }
        contador_a[0] = 0.5*(pasoinicial-dinaanterior)*float_part[0];
        contador_a[1] = 0.5*(pasoinicial-dinaanterior)*float_part[1];
        contador[0] = 0.5*(pasoinicial-dinaanterior);
        contador[1] = 0.5*(pasoinicial-dinaanterior);
        rechazo_met_mov = pasoinicial - contador_a[0] - contador_a[1];
        c_mov = pasoinicial;

        printf("%d %f %f\n",paso_racep,float_part[0],float_part[1]);
        printf("racep_0: %1.5f racep_1: %1.5f\n",1.0*contador_a[0]/contador[0],1.0*contador_a[1]/contador[1]);
        printf("DESPUES CICLO V:");
    }else{
        printf("Arreglo nuevo\n");
        dat = fopen("datos/energia.dat","w");
        fclose(dat);
        for(int j=0; j<ncomp; j++){
            printf("\r\t\t\tj: %d",j);
            part_cont=0;
            do{
                //xrand = alea_f(-R_inicial,R_inicial);
                //yrand = alea_f(-R_inicial,R_inicial);
                //xrand = distribucion_normal_5(0.0,R_inicial/5.0);
                //yrand = distribucion_normal_5(0.0,R_inicial/5.0);
                if(j==0){
                    //xrand = separacion;
                    //yrand = 0;
                    xrand = distribucion_normal_5(-separacion,10);
                    yrand = distribucion_normal_5(0,10);
                }else if(j==1){
                    //xrand = -separacion;
                    //yrand = 0;
                    xrand = distribucion_normal_5(separacion,10);
                    yrand = distribucion_normal_5(0,10);
                }
                if(xrand*xrand+yrand*yrand<R_inicial*R_inicial){
                //if(part_a_mallado(xrand,yrand)<=nceldas&&part_a_mallado(xrand,yrand)>=1){
                //if( ( xrand*xrand+yrand*yrand < ( R_inicial - 1 )*( R_inicial - 1 ) )&&( xrand*xrand+yrand*yrand > ( R_inicial - 3 )*( R_inicial - 3 ) ) ){
                    if((10*part_cont)%n_grupos[j]==0)printf("\rpart_cont: %d",part_cont);
                    part_cont++;
                    //for(int j=0; j<ncomp; j++){
                        part_plasma[part_cont][j].x = xrand;
                        part_plasma[part_cont][j].y = yrand;
                        //theta = atan2(yrand,xrand); r = norma(xrand,yrand);
                        part_plasma[part_cont][j].vx = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[j]/masa[j]));
                        part_plasma[part_cont][j].vy = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[j]/masa[j]));
                        //part_plasma[part_cont][j].vx = carga[j]*( r*esc*qe*B/masa[j] )*sin(theta);
                        //part_plasma[part_cont][j].vy = -carga[j]*( r*esc*qe*B/masa[j] )*cos(theta);
                        //part_plasma[part_cont][j].vx = alea_i(-1,1)*2*sqrt(kb*tempe[j]/masa[j]);
                        //part_plasma[part_cont][j].vy = alea_i(-1,1)*2*sqrt(kb*tempe[j]/masa[j]);
                        part_plasma[part_cont][j].v = norma(part_plasma[part_cont][j].vx,part_plasma[part_cont][j].vy);
                        part_plasma[part_cont][j].part = part_por_grupo[j];
                        if(j==3){
                            printf("\nnpart_%d: %lld part_%d.part: %lld",j,part_por_grupo[j],j,part_plasma[part_cont][j].part);
                            getchar();
                        }
                    //}
                }
            }while(part_cont<n_grupos[j]);
            //n_grupos[j]=part_cont;
        }
        rmin[0]=part_plasma[1][0].x;
        rmin[1]=part_plasma[1][0].y;
        rmax[0]=part_plasma[1][0].x;
        rmax[1]=part_plasma[1][0].y;
        /*int celdaa = part_a_mallado(-5,0,6,1,0);
        celdaa=1;
        part_plasma[celdaa][0].x=-5;
        part_plasma[celdaa][0].y=0;
        part_plasma[celdaa][0].vx = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[0]/masa[0]));
        part_plasma[celdaa][0].vy = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[0]/masa[0]));
        part_plasma[celdaa][0].v = norma(part_plasma[celdaa][0].vx,part_plasma[celdaa][0].vy);
        part_plasma[celdaa][0].part=npart[0]/2;

        //celdaa = part_a_mallado(0,-5,8,2,0);
        celdaa=2;
        part_plasma[celdaa][0].x=0;
        part_plasma[celdaa][0].y=-5;
        part_plasma[celdaa][0].vx = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[0]/masa[0]));
        part_plasma[celdaa][0].vy = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[0]/masa[0]));
        part_plasma[celdaa][0].v = norma(part_plasma[celdaa][0].vx,part_plasma[celdaa][0].vy);
        part_plasma[celdaa][0].part=npart[0]/2;

        //celdaa = part_a_mallado(5,0,7,1,1);
        celdaa=1;
        part_plasma[celdaa][1].x=5;
        part_plasma[celdaa][1].y=0;
        part_plasma[celdaa][1].vx = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[1]/masa[1]));
        part_plasma[celdaa][1].vy = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[1]/masa[1]));
        part_plasma[celdaa][1].v = norma(part_plasma[celdaa][1].vx,part_plasma[celdaa][1].vy);
        part_plasma[celdaa][1].part=npart[1]/2;

        //celdaa = part_a_mallado(0,5,7,2,1);
        celdaa=2;
        part_plasma[celdaa][1].x=0;
        part_plasma[celdaa][1].y=5;
        part_plasma[celdaa][1].vx = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[1]/masa[1]));
        part_plasma[celdaa][1].vy = 1e3*distribucion_normal_5(0.0,sqrt(kb*tempe[1]/masa[1]));
        part_plasma[celdaa][1].v = norma(part_plasma[celdaa][1].vx,part_plasma[celdaa][1].vy);
        part_plasma[celdaa][1].part=npart[1]/2;*/
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
    cout << "\nDespues de importar" << endl;
    calc_dens(0);
    printf("\nDESPUES DE CALC_DENS");
    long long int cargatotaaal = 0;
    for(int i=1;i<=nceldas;i++){
        cargatotaaal += matriz_plasma[i].carga;
    }
    if(debug==2)cout << endl << "QUE SHOW" << endl;
    /*if(cargatotaaal!=0){
        dat = fopen("pruebas/getchares.dat","a");
        fprintf(dat,"\ncargatotal: %lld ",cargatotaaal);
        fprintf(dat,"Dentro de arreglo_inicial, dentro de if cargatotaal.\n");
        fclose(dat);
        getchar();
    }*/
    p = 0;
    actu_salida();
    for(int j=0; j<ncomp; j++)hacer_histograma(j,0);
    calc_contador();
    printf("terma: %d ASDASD\n", terma);
    //getchar();getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void importar_part_a_mallado(void){
    int icelda_;
    int icelda[4]={0};
    long long int fraccion[4]={0};
    double area[4]={0};
    for(int i=1; i<=nceldas; i++){
        for(int j=0; j<ncomp; j++)matriz_plasma[i].part[j] = 0;
        matriz_plasma[i].carga = 0;
    }
    for(int j=0; j<ncomp; j++){
        for(int i=1;i<=n_grupos[j];i++){
            /*icelda_ = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,0,0);
            matriz_plasma[icelda_].part[j]+=part_plasma[i][j].part;
            matriz_plasma[icelda_].carga+=carga[j]*part_plasma[i][j].part;*/

            icelda[0] = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,0,0);
            icelda[1] = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,1,0);
            icelda[2] = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,0,1);
            icelda[3] = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,1,1);
            area[0] = (matriz_plasma[icelda[3]].x-part_plasma[i][j].x)*(matriz_plasma[icelda[3]].y-part_plasma[i][j].y);
            fraccion[0] = reso*reso*part_plasma[i][j].part*area[0];
            area[1] = (part_plasma[i][j].x-matriz_plasma[icelda[2]].x)*(matriz_plasma[icelda[2]].y-part_plasma[i][j].y);
            fraccion[1] = reso*reso*part_plasma[i][j].part*area[1];
            area[2] = (matriz_plasma[icelda[1]].x-part_plasma[i][j].x)*(part_plasma[i][j].y-matriz_plasma[icelda[1]].y);
            fraccion[2] = reso*reso*part_plasma[i][j].part*area[2];
            area[3] = (part_plasma[i][j].x-matriz_plasma[icelda[0]].x)*(part_plasma[i][j].y-matriz_plasma[icelda[0]].y);
            fraccion[3] = reso*reso*part_plasma[i][j].part*area[3];
            /*if(1=2){
                if(fraccion[0]+fraccion[1]+fraccion[2]+fraccion[3]!=part_plasma[i][j].part){
                    int fraccion_max = 0;
                    for(int k=1; k<4; k++){
                        if(fraccion[k]>fraccion_max)fraccion_max=k;
                    }
                    for(int k=0; k<4; k++)fraccion[k]+=(part_plasma[i][j].part-fraccion[0]-fraccion[1]-fraccion[2]-fraccion[3])/4;
                    fraccion[fraccion_max]+=(part_plasma[i][j].part-fraccion[0]-fraccion[1]-fraccion[2]-fraccion[3])%4;
                }
                //fraccion[3] = (part_plasma[i][j].part>( fraccion[0]+fraccion[1]+fraccion[2] ))?(part_plasma[i][j].part-fraccion[0]-fraccion[1]-fraccion[2]):0;
                if(part_plasma[i][j].part!=(fraccion[0]+fraccion[1]+fraccion[2]+fraccion[3])){
                    printf("\npart_%d_%d: %lld frac_0: %lld frac_1: %lld frac_2: %lld frac_3: %lld suma: %lld",i,j,part_plasma[i][j].part,fraccion[0],fraccion[1],fraccion[2],fraccion[3],fraccion[0]+fraccion[1]+fraccion[2]+fraccion[3]);
                    //getchar();
                }
            }
            for(int k=0; k<4; k++){
                if(fraccion[k]<0){
                    printf("\nfraccion_%d: %lld",k,fraccion[k]);
                    for(int l=0; l<4; l++)imprimir_celda_plasma_rec(icelda[l]);
                    imprimir_celda_plasma_part(i,j);
                    //getchar();
                }
            }*/
            for(int k=0; k<4; k++){
                matriz_plasma[icelda[k]].part[j]+=fraccion[k];
                matriz_plasma[icelda[k]].carga+=carga[j]*fraccion[k];
            }
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
void importar_part_a_mallado_lineal(void){
    int icelda;
    int icelda_0, icelda_1, icelda_2, icelda_3;
    long long int fraccion_0, fraccion_1, fraccion_2;
    for(int i=1; i<=nceldas; i++){
        for(int j=0; j<ncomp; j++)matriz_plasma[i].part[j] = 0;
        matriz_plasma[i].carga = 0;
    }
    for(int j=0; j<ncomp; j++){
        for(int i=1;i<=n_grupos[j];i++){
            icelda_0 = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,0,0);
            icelda_1 = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,1,0);
            icelda_2 = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,0,1);
            icelda_3 = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,1,1);
            fraccion_0 = part_plasma[i][j].part*(matriz_plasma[icelda_3].x-part_plasma[i][j].x)*(matriz_plasma[icelda_3].y-part_plasma[i][j].y)/(matriz_plasma[icelda_3].anchow*matriz_plasma[icelda_3].anchow);
            fraccion_1 = part_plasma[i][j].part*(part_plasma[i][j].x-matriz_plasma[icelda_2].x)*(matriz_plasma[icelda_2].y-part_plasma[i][j].y)/(matriz_plasma[icelda_2].anchow*matriz_plasma[icelda_2].anchow);
            fraccion_2 = part_plasma[i][j].part*(matriz_plasma[icelda_1].x-part_plasma[i][j].x)*(part_plasma[i][j].y-matriz_plasma[icelda_1].y)/(matriz_plasma[icelda_3].anchow*matriz_plasma[icelda_3].anchow);
            printf("\npart_%d_%d: %lld frac_0: %lld frac_1: %lld frac_2: %lld",i,j,part_plasma[i][j].part,fraccion_0,fraccion_1,fraccion_2);
            getchar();
            //icelda = part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,3,i,j,0,0);
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
int part_a_mallado(float a, float b, int c, int d, int e, int f, int g){
    int ix, iy;
    ix = floor( a*reso ) + floor(R*reso) + f + 1;
    iy = floor( b*reso ) + floor(R*reso) + g + 1;
    if( ( ix <= 0 )||( iy <= 0 )||( ix > (R_reso+1) )||( iy > (R_reso+1) ) ){
        if(c>=20){
            ix = 10;
            iy = 10;
            rechazo=1;
        }else{
            dat = fopen("pruebas/error_mallado.dat","w");
            fprintf(dat,"ix: %d iy: %d a: %e b: %e c: %d",ix,iy,a,b,c);
            fclose(dat);
            printf("\nix: %d iy: %d a: %e b: %e c: %d d: %d e: %d x: %e y: %e",ix,iy,a,b,c,d,e,part_plasma[d][e].x,part_plasma[d][e].y);
            printf("\nx0: %e y0: %e vx0: %e vy0: %e",x00[0],x00[1],x00[2],x00[3]);
            dat = fopen("pruebas/getchares.dat","a");
            fprintf(dat,"Dentro de part_a_mallado, unico if viene de %d. paso: %d\n",c,p);
            fclose(dat);
            getchar();
        }
    }
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
        if(int( int_clases*(part_plasma[i][a].vx-vxmin)/(vxmax-vxmin) )+1<0){
            dat = fopen("pruebas/getchares.dat","a");
            fprintf(dat,"\nVXX: %e\tvxmin: %e\tvxmax: %e entero: %d",part_plasma[i][a].vx,vxmin,vxmax,int( int_clases*(part_plasma[i][a].vx-vxmin)/(vxmax-vxmin) )+1);
            fprintf(dat,"Dentro de hacer_histograma, dentro de if entrada < 0 x de v. paso: %d\n",p);
            fclose(dat);
            getchar();
        }
        if(int( int_clases*(part_plasma[i][a].vy-vymin)/(vymax-vymin) )+1<0){
            dat = fopen("pruebas/getchares.dat","a");
            fprintf(dat,"\nVYY: %e\tvymin: %e\tvymax: %e entero: %d",part_plasma[i][a].vy,vymin,vymax,int( int_clases*(part_plasma[i][a].vy-vymin)/(vymax-vymin) )+1);
            fprintf(dat,"Dentro de hacer_histograma, dentro de if entrada < 0 y de v. paso: %d\n",p);
            fclose(dat);
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
    tipo = alea_i( 0, ncomp-1 );
    //printf("\ntipo: %d contador_tipo: %d",tipo,contador[tipo]);
    contador[tipo]++;
    //printf("\ncontador_tipo: %d",contador[tipo]);
    alea_part = alea_i( 1, n_grupos[tipo] );

    if(blineal){

        celdas_mc.finales[0] = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 20, alea_part, tipo, 0, 0 );
        celdas_mc.finales[1] = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 21, alea_part, tipo, 1, 0 );
        celdas_mc.finales[2] = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 22, alea_part, tipo, 0, 1 );
        celdas_mc.finales[3] = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 23, alea_part, tipo, 1, 1 );
        for(int i=0; i<4; i++){
            for(int j=i+1; j<4; j++){
                if(celdas_mc.finales[i]==celdas_mc.finales[j]){
                    imprimir_celda_plasma_part(alea_part,tipo);
                    imprimir_celda_plasma_rec(celdas_mc.finales[i]);
                    imprimir_celda_plasma_rec(celdas_mc.finales[j]);
                }
            }
        }
        for(int i=0; i<4; i++)matriz_plasma[celdas_mc.iniciales[i]]=matriz_plasma[celdas_mc.finales[i]];
    }else{
        n2 = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 0, alea_part, tipo, 0, 0 );
        matriz_plasma[n2i]=matriz_plasma[n2];
    }
    part_plasma[0][tipo] = part_plasma[alea_part][tipo];

    do{
        if(rechazo==1)part_plasma[alea_part][tipo]=part_plasma[0][tipo];
        if(blineal){
            dxx = alea_f(-dx,dx);
            dyy = alea_f(-dy,dy);
            part_plasma[alea_part][tipo].x += dxx;
            part_plasma[alea_part][tipo].y += dyy;

            celdas_mc.finales[4] = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 24, alea_part, tipo, 0, 0 );
            celdas_mc.finales[5] = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 25, alea_part, tipo, 1, 0 );
            celdas_mc.finales[6] = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 26, alea_part, tipo, 0, 1 );
            celdas_mc.finales[7] = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 27, alea_part, tipo, 1, 1 );
            for(int i=0; i<4; i++)matriz_plasma[celdas_mc.iniciales[4+i]]=matriz_plasma[celdas_mc.finales[4+i]];
            if(rechazo==1)break;
            //if(norma(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)>=(R-12))rechazo=1;
        }else{
            do{
                if(c_pared>100){
                    printf("\nidx: %d idy: %d",idx,idy);
                }
                idx = alea_i(-int(dx*reso),int(dx*reso));
                idy = alea_i(-int(dy*reso),int(dy*reso));
            }while(idx==0&&idy==0);
            //}while(idx<(1.0/reso)&&idy<(1.0/reso));
            //printf("\nANTES tipo: %d a_part: %d px_0: %f py_0: %f mx: %f my: %f",tipo,alea_part,part_plasma[0][tipo].x,part_plasma[0][tipo].y,matriz_plasma[n2].x,matriz_plasma[n2].y);
            part_plasma[alea_part][tipo].x += idx/reso;
            part_plasma[alea_part][tipo].y += idy/reso;

            //part_plasma[alea_part][tipo].x = alea_f(-R,R);
            //part_plasma[alea_part][tipo].y = alea_f(-R,R);

            int ix = floor( part_plasma[alea_part][tipo].x*reso ) + floor(R*reso) + 1, iy = floor( part_plasma[alea_part][tipo].y*reso ) + floor(R*reso) + 1;
            if( !(( ix <= 0 )||( iy <= 0 )||( ix > (R_reso+1) )||( iy > (R_reso+1) )) ){
                n1 = part_a_mallado( part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y, 8, alea_part, tipo, 0, 0);
                matriz_plasma[n1i]=matriz_plasma[n1];
                rechazo = 0;
                if(norma(part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y)>=(R-12))rechazo=1;
            }else{
                rechazo=1;
            }
        }
        //if(part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)>nceldas||part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)<1)rechazo=1;
        //pared();
        if(c_pared>100){
            imprimir_celda_plasma_part(0,tipo);
            imprimir_celda_plasma_part(alea_part,tipo);
            dat = fopen("pruebas/getchares.dat","a");
            fprintf(dat,"\ncontador_pared: %d rechazo: %d norma: %f",c_pared,rechazo, norma(part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y),R);
            fprintf(dat,"Dentro de mover_particulas_plasma, dentro de if de c_pared. paso: %d\n",p);
            fclose(dat);
            getchar();
        }
        if(c_pared>110)rechazo=0;
        c_pared++;
    }while(rechazo==1);
    if(c_pared>110){
        dat = fopen("pruebas/getchares.dat","a");
        fprintf(dat,"\nc_pared: %d tipo: %d a_part: %d x: %f y: %f rho: %f",c_pared,tipo,alea_part,part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y,norma(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y));
        fprintf(dat,"\nx_0: %f y_0: %f rho_0: %f",part_plasma[0][tipo].x,part_plasma[0][tipo].y,norma(part_plasma[0][tipo].x,part_plasma[0][tipo].y));
        fprintf(dat,"Dentro de mover_particulas_plasma, dentro de if de c_pared fuera de do-while. paso: %d\n",p);
        fclose(dat);
        getchar();
    }
    c_pared = 0;
    rechazo = 0;
    //printf("\nDESPUES px_f: %f py_f: %f mx_f: %f my_f: %f",part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y,matriz_plasma[n1].x,matriz_plasma[n1].y);
    //getchar();

    //printf("\ntipo: %d alea_part: %d dx: %f dy: %f",tipo,alea_part,idx/reso,idy/reso);
    //imprimir_celda_plasma_part(0,tipo);
    //imprimir_celda_plasma_part(alea_part,tipo);
    //imprimir_celda_plasma_rec(n2i);
    //imprimir_celda_plasma_rec(n1i);
    if(1==2){
        for(int i=0; i<4; i++)imprimir_celda_plasma_rec(celdas_mc.finales[i]);
        for(int i=0; i<4; i++)imprimir_celda_plasma_rec(celdas_mc.finales[4+i]);
    }
    if(rechazo==0){
        if(blineal){
            float area=0;
            long long int fraccion1[4]={0}, fraccion2[4]={0};

            area = (matriz_plasma[celdas_mc.finales[3]].x-part_plasma[0][tipo].x)*(matriz_plasma[celdas_mc.finales[3]].y-part_plasma[0][tipo].y);
            fraccion1[0] = reso*reso*part_plasma[0][tipo].part*area;
            area = (part_plasma[0][tipo].x-matriz_plasma[celdas_mc.finales[2]].x)*(matriz_plasma[celdas_mc.finales[2]].y-part_plasma[0][tipo].y);
            fraccion1[1] = reso*reso*part_plasma[0][tipo].part*area;
            area = (matriz_plasma[celdas_mc.finales[1]].x-part_plasma[0][tipo].x)*(part_plasma[0][tipo].y-matriz_plasma[celdas_mc.finales[1]].y);
            fraccion1[2] = reso*reso*part_plasma[0][tipo].part*area;
            area = (part_plasma[0][tipo].x-matriz_plasma[celdas_mc.finales[0]].x)*(part_plasma[0][tipo].y-matriz_plasma[celdas_mc.finales[0]].y);
            fraccion1[3] = reso*reso*part_plasma[0][tipo].part*area;

            area = (matriz_plasma[celdas_mc.finales[7]].x-part_plasma[alea_part][tipo].x)*(matriz_plasma[celdas_mc.finales[7]].y-part_plasma[alea_part][tipo].y);
            fraccion2[0] = reso*reso*part_plasma[alea_part][tipo].part*area;
            area = (part_plasma[alea_part][tipo].x-matriz_plasma[celdas_mc.finales[6]].x)*(matriz_plasma[celdas_mc.finales[6]].y-part_plasma[alea_part][tipo].y);
            fraccion2[1] = reso*reso*part_plasma[alea_part][tipo].part*area;
            area = (matriz_plasma[celdas_mc.finales[5]].x-part_plasma[alea_part][tipo].x)*(part_plasma[alea_part][tipo].y-matriz_plasma[celdas_mc.finales[5]].y);
            fraccion2[2] = reso*reso*part_plasma[alea_part][tipo].part*area;
            area = (part_plasma[alea_part][tipo].x-matriz_plasma[celdas_mc.finales[4]].x)*(part_plasma[alea_part][tipo].y-matriz_plasma[celdas_mc.finales[4]].y);
            fraccion2[3] = reso*reso*part_plasma[alea_part][tipo].part*area;

            long long int fraccion_total1 = 0, fraccion_total2 = 0;
            fraccion_total1=fraccion1[0]+fraccion1[1]+fraccion1[2]+fraccion1[3];
            fraccion_total2=fraccion2[0]+fraccion2[1]+fraccion2[2]+fraccion2[3];

            /*printf("\n");
            for(int i=0; i<4; i++)printf("frac1_%d: %lld\t",i,fraccion1[i]);
            printf("\tfractot1: %lld\n",fraccion_total1);
            for(int i=0; i<4; i++)printf("frac2_%d: %lld\t",i,fraccion2[i]);
            printf("\tfractot2: %lld",fraccion_total2);*/

            for(int i=0; i<4; i++){
                matriz_plasma[celdas_mc.finales[i]].part[tipo]-=fraccion1[i];
                matriz_plasma[celdas_mc.finales[i]].carga-=carga[tipo]*fraccion1[i];
                matriz_plasma[celdas_mc.finales[4+i]].part[tipo]+=fraccion2[i];
                matriz_plasma[celdas_mc.finales[4+i]].carga+=carga[tipo]*fraccion2[i];
            }
        }else{
            matriz_plasma[n1].part[tipo]+=part_plasma[alea_part][tipo].part;
            matriz_plasma[n1].carga+=carga[tipo]*part_plasma[alea_part][tipo].part;
            matriz_plasma[n2].part[tipo]-=part_plasma[alea_part][tipo].part;
            matriz_plasma[n2].carga-=carga[tipo]*part_plasma[alea_part][tipo].part;
        }
        celdas_mc.size=8;
        /*printf("\nsize: %d\n",celdas_mc.size);
        for(int i=0; i<celdas_mc.size; i++)printf("c_f_%d: %d\t",i,celdas_mc.finales[i]);
        printf("\n");
        for(int i=0; i<celdas_mc.size; i++)printf("c_i_%d: %d\t",i,celdas_mc.iniciales[i]);*/
        for(int i=0; i<celdas_mc.size-4; i++){
            for(int j=0; j<4; j++){
                if(celdas_mc.finales[4+i]==celdas_mc.finales[j]){
                    //printf("\nLas celdas c_%d y c_%d son iguales.",j,4+i);
                    int i_intermedio = celdas_mc.finales[4+i], i_intermedio_2 = celdas_mc.iniciales[4+i];
                    celdas_mc.size--;
                    celdas_mc.finales[4+i]=celdas_mc.finales[celdas_mc.size];
                    celdas_mc.finales[celdas_mc.size]=i_intermedio;
                    celdas_mc.iniciales[4+i]=celdas_mc.iniciales[celdas_mc.size];
                    celdas_mc.iniciales[celdas_mc.size]=i_intermedio_2;
                    i=0;
                }
            }
        }
    }
    /*printf("\nsize: %d\n",celdas_mc.size);
    for(int k=0; k<celdas_mc.size; k++)printf("c_f_%d: %d\t",k,celdas_mc.finales[k]);
    printf("\n");
    for(int i=0; i<celdas_mc.size; i++)printf("c_i_%d: %d\t",i,celdas_mc.iniciales[i]);

    long long int dpart_total=0;
    for(int i=0; i<celdas_mc.size; i++)dpart_total+=matriz_plasma[celdas_mc.finales[i]].carga-matriz_plasma[celdas_mc.iniciales[i]].carga;
    printf("\ndpart: %lld",dpart_total);
    getchar();*/
    //imprimir_celda_plasma_rec(n2);
    //imprimir_celda_plasma_rec(n1);
    //getchar();

    //importar_part_a_mallado();
}
////////////////////////////////////////////////////////////////////////////////////
void pared(void){
    float rho = norma(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y);
    //if(part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)>nceldas||part_a_mallado(part_plasma[alea_part][tipo].x,part_plasma[alea_part][tipo].y)<1)rechazo=1;
    if(rho>=(R-12)){
        rechazo=1;
        //printf("\nQue show!");
        //getchar();
    }
}
////////////////////////////////////////////////////////////////////////////////////
//void reflejar_pared(int a,int b,float* vx, float* vy){
void reflejar_pared(float x0, float y0, float *vx, float *vy){
    float theta, v0x = *vx, v0y = *vy;
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
double de_plasma(void){
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
double de_plasma_lineal(void){
    int i, j, k;
    bool bool1 = false, b_no_nmc=true;
    //bool1 = true;
    double ei = 0, ef = 0, d;
    de_autoenergia=0;
    for(i=1;i<=nceldas;i++){
        if(matriz_plasma[i].carga!=0){
            b_no_nmc = true;
            for(j=0; j<celdas_mc.size; j++){
                if(i==celdas_mc.finales[j]){
                    b_no_nmc=false;
                    break;
                }
            }
            if(b_no_nmc){
                for(j=0; j<celdas_mc.size; j++){
                    d = distancia(celdas_mc.finales[j],i);
                    if(d<0.5){
                        printf("\nd<0.5! interaccion entre las otras celdas con las finales i: %d j: %d",i,j);
                        getchar();
                    }
                    ef += (qe*qe*matriz_plasma[celdas_mc.finales[j]].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                    d = distancia(celdas_mc.iniciales[j],i);
                    if(d<0.5){
                        printf("\nd<0.5! interaccion entre las otras celdas con las iniciales i: %d j: %d",i,j);
                        getchar();
                    }
                    ei += (qe*qe*matriz_plasma[celdas_mc.iniciales[j]].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                }
            }
        }
    }
    for(i=0; i<celdas_mc.size; i++){
        for(j=i+1; j<celdas_mc.size; j++){
            d = distancia(celdas_mc.finales[i],celdas_mc.finales[j]);
            if(d<0.5){
                printf("\nd<0.5! interaccion entre las celdas finales i: %d j: %d",i,j);
                getchar();
            }
            ef += (qe*qe*matriz_plasma[celdas_mc.finales[i]].carga*matriz_plasma[celdas_mc.finales[j]].carga)/(4*pi*epce*epsi*d*esc);
            d = distancia(celdas_mc.iniciales[i],celdas_mc.iniciales[j]);
            if(d<0.5){
                printf("\nd<0.5! interaccion entre las celdas iniciales i: %d j: %d",i,j);
                getchar();
            }
            ei += (qe*qe*matriz_plasma[celdas_mc.iniciales[i]].carga*matriz_plasma[celdas_mc.iniciales[j]].carga)/(4*pi*epce*epsi*d*esc);
        }
    }
    /*for(i=0; i<4; i++){
        for(k=i+1; k<4; k++){
            d = distancia(celdas_mc1[4+i],celdas_mc2[4+k]);
            ei += (qe*qe*matriz_plasma[celdas_mc1[4+i]].carga*matriz_plasma[celdas_mc2[4+k]].carga)/(4*pi*epce*epsi*d*esc);
            d = distancia(celdas_mc1[i],celdas_mc2[k]);
            if(d<0.5){
                printf("\nd_%d_%d: %e",celdas_mc1[i],celdas_mc2[k],d);
                imprimir_celda_plasma_rec(celdas_mc1[i]);
                imprimir_celda_plasma_rec(celdas_mc2[k]);
                imprimir_celda_plasma_part(alea_part,tipo);
                getchar();
            }
            ef += (qe*qe*matriz_plasma[celdas_mc1[i]].carga*matriz_plasma[celdas_mc2[k]].carga)/(4*pi*epce*epsi*d*esc);
        }
    }
    for(i de las primeras){
        for(de las primeras sin la i){
            energias
        }
        for(j de las segundas){
            if(j!=(todas las primeras)){
                energias
            }
        }
    }
    for(i de las segundas){
        if(i!=(todas las primeras)){
            for(j de las segundas){
                if(j!=(todas las primeras)&&j!=i){
                    energias
                }
            }
        }
    }
    int arregloceldas_inicial[8];
    for(i en arregloceldas_inicial){
        for(k=i-1;k>=0;k--)if(k==i)continue;
        for(j en arregloceldas_inicial, j=i+1){
            for(k=j-1;k>=0;k--){
                if(k==j){
                    boola = false;
                    break;
                }
            }
            if(boola){
                sumarenergia(i,j);
                contador++;
            }
        }
    }*/

    de_coulomb = ef - ei;
    /*if(p>=terma){
        for(i=0; i<celdas_mc.size; i++){
            ef += autoenergia(celdas_mc.finales[i]);
            ei += autoenergia(celdas_mc.iniciales[i]);
            de_autoenergia += autoenergia(celdas_mc.finales[i]);
            de_autoenergia -= autoenergia(celdas_mc.iniciales[i]);
        }
        if(bool1){
            printf(" de: %e\ndxx: %e dyy: %e",ef-ei,dxx,dyy);
            getchar();
        }
    }*/
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
    float zeta, argexp;
    double de;
    int i, j;

    zeta = alea();
    //de = -kb*tempe[tipo];
    if(blineal){
        de = de_plasma_lineal();
    }else{
        de = de_plasma();
    }

    argexp = -de/(kb*tempe[tipo]);
    //argexp = 200;

    if((argexp>=100?2.0:exp(argexp))>=zeta){
        contador_a[tipo]++;
        if(p%(actu/100)==0){
            dat = fopen("pruebas/denergias_acep.dat","a");
            fprintf(dat,"%d\t%e\t%e\t%e\n",p,de_coulomb/(kb*tempe[tipo]),de_autoenergia/(kb*tempe[tipo]),dem/(kb*tempe[tipo]));
            fclose(dat);
        }
        //dt_termico+=(esc*distancia(n1,n2)*part_plasma[alea_part][tipo].part)/part_plasma[alea_part][tipo].v;
        float distancia_termica = norma(part_plasma[alea_part][tipo].x-part_plasma[0][tipo].x,part_plasma[alea_part][tipo].y-part_plasma[0][tipo].y);
        //dt_termico+=(esc*distancia(n1,n2)*part_plasma[alea_part][tipo].part)/sqrt( (8*kb*tempe[tipo])/(pi*masa[tipo]) );
        dt_termico+=(esc*distancia_termica*part_plasma[alea_part][tipo].part)/sqrt( (8*kb*tempe[tipo])/(pi*masa[tipo]) );
        //printf("\ndttermico: %e\tdist: %e\tpart: %lld\tn1: %d\tn2: %d\talea_part: %d",dt_termico,distancia(n1,n2),part_plasma[alea_part][tipo].part,n1,n2,alea_part);
        //getchar();
        energia_con_de+=de;
        //printf("\np: %d d: %f part_desp: %lld",p,distancia(n1,n2),part_plasma[alea_part][tipo].part);
        //getchar();
        matriz_plasma[n2].mov_a++;
    }else{
        if(p%(actu/100)==0){
            dat = fopen("pruebas/denergias_rechazo.dat","a");
            fprintf(dat,"%d\t%e\t%e\t%e\n",p,de_coulomb/(kb*tempe[tipo]),de_autoenergia/(kb*tempe[tipo]),dem/(kb*tempe[tipo]));
            fclose(dat);
        }
        part_plasma[alea_part][tipo] = part_plasma[0][tipo];
        if(blineal){
            for(i=0; i<8; i++)matriz_plasma[celdas_mc.finales[i]]=matriz_plasma[celdas_mc.finales[i]];
        }else{
            matriz_plasma[n1] = matriz_plasma[n1i];
            matriz_plasma[n2] = matriz_plasma[n2i];
        }

        //importar_part_a_mallado();
        rechazo_met_mov++;
    }
}
////////////////////////////////////////////////////////////////////////////////////
//FORMATO posiciones.dat:
//paso terma ncomp racep(0,1,2,...) cdina
//info part_plasma
void actu_salida(void){
	int max_n_grupos = n_grupos[0];
	long long int cargatotaaal=0,cargamin=matriz_plasma[1].carga,cargamax=matriz_plasma[1].carga;
	if(p%1000000==0&&p>1){
        for(int i=0; i<ncomp; i++ )if(n_grupos[i]>=max_n_grupos)max_n_grupos = n_grupos[i];
        sprintf(salidac,"datos/posiciones.dat");
        if( dat=fopen(salidac,"w") ){
            fprintf(dat,"%d %d %d ",p,terma,ncomp);
            for(int i=0;i<ncomp;i++)fprintf(dat,"%f ",(contador[i]>0)?1.0*contador_a[i]/contador[i]:0.0 );
            fprintf(dat,"%d",c_dina);
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
            /*for(int i=1;i<=nceldas;i++){
                fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[j]);
            }*/

            for(int i=-int(R*reso);i<=int(R*reso);i++){
                for(int k=-int(R*reso);k<=int(R*reso);k++){
                    //if(abs(i)<=(R-5)*reso&&abs(k)<=(R-5)*reso){
                    if(i*i+k*k<=(R-1)*reso*(R-1)*reso){
                        int n_celda_actu = nmatriz_plasma[i+int(R*reso)+1][k+int(R*reso)+1];
                        fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[n_celda_actu].x, matriz_plasma[n_celda_actu].y, matriz_plasma[n_celda_actu].part[j]);
                    }else{
                        int n_celda_actu = nmatriz_plasma[i+int(R*reso)+1][k+int(R*reso)+1];
                        //fprintf(dat,"%f\t%f\t?\n", matriz_plasma[n_celda_actu].x, matriz_plasma[n_celda_actu].y);
                    }
                }
                fprintf(dat,"\n");
            }

            fclose(dat);
        }
    }
    if( iprint[ncomp]==1 ){
        sprintf(salidac,"datos/todas%i.dat",(p-terma)/actu);
        dat=fopen(salidac,"w");
        fprintf(dat, "#X\tY\tZ\n");
        for(int i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[0]+matriz_plasma[i].part[1]+matriz_plasma[i].part[2]);
        }
        fclose(dat);
    }

    if( iprint[ncomp+1]==1 ){
        sprintf(salidac,"datos/carga%i.dat",(p-terma)/actu);
        dat=fopen(salidac,"w");
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
void actu_salida_dina(int a){
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
        if(a==0)sprintf(salidac,"datos/dina0_part_%d_%d.dat",j,p);
        else sprintf(salidac,"datos/dina1_part_%d_%d.dat",j,p);
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
    if(a==0)sprintf(salidac,"datos/dina0_todas%i.dat",p);
    else sprintf(salidac,"datos/dina1_todas%i.dat",p);
    if( iprint[ncomp]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(int i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[0]+matriz_plasma[i].part[1]+matriz_plasma[i].part[2]);
        }
        fclose(dat);
    }

    if(a==0)sprintf(salidac,"datos/dina0_carga%i.dat",p);
    else sprintf(salidac,"datos/dina1_carga%i.dat",p);
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

        if(a==0)sprintf(salidac,"datos/dina0_carga%i_m.dat",p);
        else sprintf(salidac,"datos/dina1_carga%i_m.dat",p);
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
	long long int gtotal[R_reso*R_reso+1+2]={0}, gcarga[R_reso*R_reso+1+2]={0};
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
        for(int j=0; j<ncomp-1; j++){
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
double energia(void){
    int i,j;
    double energia_total=0, ecoulomb_total=0, aenergia_total=0, energia_i=0, d;
    for(i=1;i<=nceldas;i++){
        energia_i = 0;
        for(j=i+1;j<=nceldas;j++){
            if((matriz_plasma[i].carga!=0)&&(matriz_plasma[j].carga!=0)){
                if(i!=j){
                    //d = distanciacelda(i,j);
                    d = distancia(i,j);
                    energia_i += (qe*qe*matriz_plasma[i].carga*matriz_plasma[j].carga)/(4*pi*epce*epsi*d*esc);
                }
            }
        }
        //aenergia_total += autoenergia(i);
        ecoulomb_total += energia_i;
    }
    energia_total = aenergia_total + ecoulomb_total;
    if(p>=pasoinicial||p==0){
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
    float x,y,vx,vy,ax,ay;
    float x0, v0x, y0, v0y;
    float Ex=0.0,Ey=0.0,Bz=B;
    float dt, normar;
    int ktau = 200, contadork, contador_dina, cambios_de_celda[2]={0};
    float suma_dist[2]={0};

    if(solodinamica){
        if(B>1e-5){
            dt = ((2*pi*min_mass)/(100.0*max_charge*qe*B));
        }else{
            float max_v = part_plasma[1][0].v;
            for(int j=0; j<ncomp; j++){
                for(int i=1; i<=n_grupos[j]; i++){
                    if(part_plasma[i][j].v>=max_v)max_v=part_plasma[i][j].v;
                }
            }
            dt = (1.0/(1.0*max_v*reso));
            printf("\tdt: %e",dt);
        }
    }else{
        if(dt_termico>0){
            dt = dt_termico/(npart[0]+npart[1]);
            //printf("\nTERMICO\tp: %d dt: %e dt_termico: %e",p,dt,dt_termico);
            //getchar();
            dt_termico=0;
        }else{
            dt = dtau;
            //printf("\nDTAU\tp: %d dt: %e dt_termico: %e",p,dt,dt_termico);
            //getchar();
        }
    }
    dt *= 0.005;
    //Zprintf("\nmasa[%d]: %e masa[%d]: %e masa[%d]: %e min_mass: %e",0,masa[0],1,masa[1],2,masa[2],min_mass);
    //printf("\ncarga[%d]: %d carga[%d]: %d carga[%d]: %d max_charge: %d",0,carga[0],1,carga[1],2,carga[2],max_charge);
    //printf("\nperiodo_min: %e dt_sincampo: %e dt(B): %e\n", ((2*pi*min_mass)/(100.0*max_charge*qe*B)), (esc/(100.0*max_v*reso)), B>1e-5?((2*pi*min_mass)/(100.0*max_charge*qe*B)):(esc/(100.0*max_v*reso)) );
    //getchar();
    //float dtau = B>1e-5?((2*pi*min_mass)/(100.0*max_charge*qe*B)):(esc/(100.0*max_v*reso));
    if(debug==3)printf("\nANTES CICLO\n");
    for(int j=0; j<ncomp; j++){
        //for(int i=0; i<ncomp; i++){
            //printf("\nnpart2_%d: %d",i,n_grupos[i]);
            //getchar();
        //}
        if(n_grupos[j]>0){
            for(int i=1; i<=n_grupos[j]; i++){
                contador_dina=0;
                contadork = 0;
                if(debug==3&&i%100==0)printf("\rDENTRO");
                //printf("\nANTES DE DINAMICA i: %d j: %d\n",i,j);
                //if((10*i)%n_grupos[j]==0)printf("\ri: %d\n",i);
                //dinamica();
                //dt = dtau/100.0;
                x00[0] = part_plasma[i][j].x;
                x00[1] = part_plasma[i][j].y;
                x00[2] = part_plasma[i][j].vx;
                x00[3] = part_plasma[i][j].vy;

                part_plasma[i][j].celda_antes=part_a_mallado(x00[0],x00[1],15,i,j,0,0);

                x0=part_plasma[i][j].x;
                v0x=part_plasma[i][j].vx;
                y0 = part_plasma[i][j].y;
                v0y = part_plasma[i][j].vy;


                //float periodo = (2*pi*mass)/(abs(charge)*qe*Bz);
                //printf("\nq: %e m: %e q/m: %e",charge*qe,mass*m_u,charge*qe/mass/m_u);
                //printf("\nipart: %d p_x: %f p_y: %f x0: %f y0: %f",ipart,part_plasma[1][0].x,part_plasma[1][0].y,x0,y0);
                //printf("\ncarga: %d masa: %e charge: %d mass: %e",carga[0],masa[0],charge,mass);

                //dat=fopen("dinamica.dat","w");
                //fprintf(dat,"#t\tx\tvx\ty\tvy\n");
                //fprintf(dat,"%f\t%f\t%f\t%f\t%f\n",0.0,x0,v0x,y0,v0y);
                //for(int i=1;i<=int(periodo/dt)+1;i++){
                //long long int k;
                for(int k=1;k<=ktau;k++){
                    if(debug>=4&&k%1==0){
                        contadork++;
                        printf("\rDENTRO CICLO K j: %d i: %d rho: %f\t\tcontador: %d",j,i,norma(x,y)/esc,contadork);
                    }
                    /*if(contador_dina==100){
                        x0=part_plasma[i][j].x*esc;
                        v0x=part_plasma[i][j].vx;
                        y0 = part_plasma[i][j].y*esc;
                        v0y = part_plasma[i][j].vy;
                        dat = fopen("puntos.dat","w");
                        fclose(dat);
                    }*/
                    if(contador_dina>=100){
                        dat = fopen("pruebas/contador_dina.dat","a");
                        fprintf(dat,"\ni: %d j: %d cont: %d",i,j,contador_dina);
                        fclose(dat);
                        k-=1;
                        contador_dina=0;
                        x0 *= 0.9;
                        y0 *= 0.9;
                    }
                    Ex=0.0; Ey=0.0; Bz=B;
                    //normar = norma(x0,y0)/esc;
                    //Ex=0.0; Ey=0.0; Bz=4.0*B*( (normar>5&&normar<6)?(-2*pow(normar-5,3)+3*pow(normar-5,2)):0+(normar>6)?1.0:0 );
                    //Ex=0.0; Ey=0.0; Bz=B*( 1+tanh((normar-4)) );
                    //Ex=0.0; Ey=0.0; Bz=2*B*sqrt( pow((x0+(R-2)*esc),2)+pow((y0+(R-2)*esc),2) )/(R*esc);
                    //Ex=carga[j]*x0*1000.0/(2*(R-2)*esc*esc); Ey=carga[j]*y0*1000.0/(2*(R-2)*esc*esc); Bz=B;
                    //Ex=0.0; Ey=p<50?50.0/(2*(R-2)*esc*esc):0.0; Bz=B;
                    ax = ((qe*carga[j])/masa[j])*(Ex+v0y*Bz);
                    ay = ((qe*carga[j])/masa[j])*(Ey-v0x*Bz);
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
                    /*if(contador_dina>100){
                        printf("\nx_0: %f y_0: %f vx_0: %f vy_0: %f",part_plasma[i][j].x*esc,part_plasma[i][j].y*esc,part_plasma[i][j].vx,part_plasma[i][j].vy);
                        printf("\nx: %f y: %f rho: %f R-2: %f vx: %f vy: %f vr: %f",x,y, norma(x,y), (R-2)*esc,vx,vy,(vx*cos( atan2(y,x) )+vy*sin( atan2(y,x) )) );
                        dat = fopen("puntos.dat","a");
                        fprintf(dat,"%f %f 1 %f %f\n",part_plasma[i][j].x*esc,part_plasma[i][j].y*esc,part_plasma[i][j].vx/part_plasma[i][j].v,part_plasma[i][j].vy/part_plasma[i][j].v);
                        fprintf(dat,"%f %f 1 %f %f\n",x,y,vx/norma(vx,vy),vy/norma(vx,vy));
                        fclose(dat);
                        //getchar();
                    }*/
                    if( norma(x,y)>=(R-12) ){
                        /*dat = fopen("vectores.dat","w");
                        dat2 = fopen("vectoresfinales.dat","w");
                        printf("\nANTES v0x: %f v0y: %f",v0x, v0y);
                        reflejar_pared(i,j,&v0x,&v0y);
                        printf("\nDESPUES v0x: %f v0y: %f",v0x, v0y);

                        fclose(dat);
                        fclose(dat2);
                        printf("\nQue show");getchar();*/
                        k-=1;
                        contador_dina++;
                        if((v0x*cos( atan2(y0,x0) )+v0y*sin( atan2(y0,x0) ))>=0){
                            reflejar_pared(x0,y0,&v0x,&v0y);
                            if((v0x*cos( atan2(y0,x0) )+v0y*sin( atan2(y0,x0) ))>=0){
                                printf("\nVr positiva: %e",(v0x*cos( atan2(y0,x0) )+v0y*sin( atan2(y0,x0) )));
                                dat = fopen("pruebas/getchares.dat","a");
                                fprintf(dat,"Dentro de dinamica, segundo if vr positiva. paso: %d\n",p);
                                fclose(dat);
                                getchar();
                            }
                            /*if(contador_dina>101){
                                printf("\nSI SE REFLEJO");
                                printf("\nx: %f y: %f vx: %f vy: %f vr: %f",x0,y0,v0x,v0y, v0x*cos(atan2(y0,x0))+v0y*sin(atan2(y0,x0)) );
                                getchar();
                            }*/
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
                /*printf("\ni: %d dx: %f dy: %f distancia: %f",i,fabs( part_plasma[i][j].x - x/esc ), fabs( part_plasma[i][j].y - y/esc ), norma( part_plasma[i][j].x - x/esc, part_plasma[i][j].y - y/esc ) );
                printf("\nvx_i: %f vy_i: %f vx_f: %f vy_f: %f", part_plasma[i][j].vx, part_plasma[i][j].vy, vx, vy );
                getchar();*/
                part_plasma[i][j].x = x;
                part_plasma[i][j].y = y;
                part_plasma[i][j].vx = vx;
                part_plasma[i][j].vy = vy;
                part_plasma[i][j].v = norma(vx,vy);

                suma_dist[j] += norma(x00[0]-x,x00[1]-y);
                int celda_despues = part_a_mallado(x,y,16,i,j,0,0);
                if(celda_despues!=part_plasma[i][j].celda_antes){
                    cambios_de_celda[j]++;
                    /*printf("\ni: %d j: %d cel_antes: %d cel_despues: %d",i,j,part_plasma[i][j].celda_antes,celda_despues);
                    imprimir_celda_plasma_rec(part_plasma[i][j].celda_antes);
                    imprimir_celda_plasma_rec(celda_despues);
                    printf("\nx0: %f y0: %f x: %f y: %f",x00[0],x00[1],x,y);
                    getchar();*/
                }
            }
            for(int i=1; i<=n_grupos[j]; i++){
                //printf("\nANTES DE DINAMICA i: %d j: %d\n",i,j);
                //if(norma(part_plasma[i][j].x,part_plasma[i][j].y)>R){
                if(part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,1,i,j,0,0)>nceldas||part_a_mallado(part_plasma[i][j].x,part_plasma[i][j].y,1,i,j,0,0)<1){
                    dat = fopen("pruebas/getchares.dat","a");
                    fprintf(dat,"Dentro de dinamica, else de if ngrupos. paso: %d\n",p);
                    fprintf(dat,"i: %d n_grupos[%d]: %d rho2: %f\n",i,j,n_grupos[j],norma(part_plasma[i][j].x,part_plasma[i][j].y));
                    fclose(dat);
                    getchar();
                }
            }
        }else{
            printf("\nn_grupos[%d]: %d",j,n_grupos[j]);
            dat = fopen("pruebas/getchares.dat","a");
            fprintf(dat,"Dentro de dinamica, else de if ngrupos. paso: %d\n",p);
            fclose(dat);
            getchar();
        }
    }
    if(part_plasma[1][0].x<rmin[0])rmin[0]=part_plasma[1][0].x;
    if(part_plasma[1][0].y<rmin[1])rmin[1]=part_plasma[1][0].y;
    if(part_plasma[1][0].x>rmax[0])rmax[0]=part_plasma[1][0].x;
    if(part_plasma[1][0].y>rmax[1])rmax[1]=part_plasma[1][0].y;
    if(debug==2)printf("\nDESPUES CICLO");
    dat = fopen("pruebas/cambios_celda.dat","a");
    fprintf(dat,"%e\t%d\t%d\t%d\t%d\t%f\t%f\n",ktau*dt,p,n_grupos[0],cambios_de_celda[0],cambios_de_celda[1],1.0*suma_dist[0]/n_grupos[0],1.0*suma_dist[1]/n_grupos[1]);
    fclose(dat);
    dat = fopen("pruebas/r_larmor.dat","a");
    fprintf(dat,"%d %f %f %f %f\n",p,rmin[0],rmax[0],rmin[1],rmax[1]);
    fclose(dat);
    if(debug==2)printf("\nDESPUES CAMBIO DE CELDA");
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
    if(win_os)system("mkdir datos\\promedios");
    else system("mkdir datos/promedios");
    if(win_os)system("mkdir datos\\PNG");
    else system("mkdir datos//PNG");
    system("mkdir histogramas");
    system("mkdir pruebas");

    /*dat = fopen("pruebas/naleprom.dat","w");
    fclose(dat);
    dat = fopen("pruebas/nale_f.dat","w");
    fclose(dat);*/
    dat = fopen("pruebas/energias.dat","w");
    fclose(dat);
    dat = fopen("pruebas/denergias_acep.dat","w");
    fclose(dat);
    dat = fopen("pruebas/denergias_rechazo.dat","w");
    fclose(dat);
    dat = fopen("pruebas/razon_acep.dat","w");
    fclose(dat);
    /*dat = fopen("pruebas/rms_part.dat","w");
    fclose(dat);*/
    dat = fopen("pruebas/carga_cero.dat","w");
    fclose(dat);
    dat = fopen("datos/cargatotal.dat","w");
    fclose(dat);
    dat = fopen("pruebas/carga.dat","w");
    fclose(dat);
    dat = fopen("pruebas/getchares.dat","w");
    fclose(dat);
    dat = fopen("pruebas/cambios_celda.dat","w");
    fprintf(dat,"#dt\tpaso\tngrupos\tcambios_0\tcambios_1\td0_prom\td1_prom\n");
    fclose(dat);
    dat = fopen("pruebas/energia_dinamica.dat","w");
    fclose(dat);
    dat = fopen("pruebas/r_larmor.dat","w");
    fclose(dat);
    dat = fopen("pruebas/tiempo.dat","w");
    fclose(dat);

    for(int i=1;i<=R_reso*R_reso+1;i++){
        matriz_plasma[i].carga=matriz_plasma[i].part[0]=matriz_plasma[i].part[1]=matriz_plasma[i].part[2]=matriz_plasma[i].part[3]=matriz_plasma[i].mov_a=0;
        matriz_plasma[i].phi=matriz_plasma[i].rho=matriz_plasma[i].x=matriz_plasma[i].y=0;
        for(int j=0;j<=3;j++){
            //matriz_plasma[i].t[j]=0;
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
    int int_a = *a;
    float float_b = *b;
    printf("\nint: %d float: %f",int_a,float_b);
    *a=int_a+1;
    *b=float_b+1;
    printf("\nint: %d float: %f",*a,*b);
}
////////////////////////////////////////////////////////////////////////////////////
