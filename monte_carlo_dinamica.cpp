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
void leer_datos_iniciales2(void);
void leer_campo_magnetico(void);
void imprimir_datos_iniciales(void);
void imprimir_datos_iniciales2(void);
void imprimir_celda_plasma(int a);
void imprimir_celda_plasma_rec(int a);
void imprimir_celda_plasma_part(int a, int b);
void imprimir_celda_plasma_vec(int a, int b);
void imprimir_celda_plasma_vec2(int a, int b);
void condiciones_iniciales(void);
void condiciones_iniciales2(void);
void arreglo_inicial(void);
void arreglo_inicial_part(void);
void importar_part_a_mallado(void);
int part_a_mallado(int a, int b);
void intercambiar(int a, int b, int c );
void crear_matriz_plasma(void);
void crear_matriz_plasma_rec(void);
void crear_matriz_plasma_rec2(void);
void vecinos(int a);
void calc_carga(int a);
void hacer_histograma(int a, int b);
void hacer_histograma2(int a, int b);
void hacer_distribucion(int a);
void hacer_rms(int a);
void vel_txt(char a[10],int b);
void vel_txt2(char a[10],int b);
void crear_archivos_iniciales_en_blanco(void);

float distanciacelda(int a, int b);
float distancianormal(int a, int b);
float dist_ima(int a, int b, int c);
int signo(float a);
float norma(float a, float b);
void mover_particulas(void);
void mover_particulas_new(void);
void mover_particulas_part(void);
void pared(void);
void metropolis_plasma(void);
void dinamica(void);
void dinamica2(void);
void retermalizacion(void);

float de_plasma(void);
float calcular_de_mov(void);
float autoenergia(int a);
void calcular_momento_total(void);
void gdr_plasma(void);
void gdr_plasma2(void);
void actu_salida(void);
void actu_salida2(void);
void salida_prom(void);
void salida_prom2(void);
float energia(void);
float energia_cinetica(int a);
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
int celda_dir_rec_tipob(int a, int b);
void celda_dir_rec_tipob_error(int a, int b);
void mainsillo(void);

void printposiciones(int a);

//////////////////////////////Constantes
int pasos, actu, terma, pasoinicial, termaanterior, iprint[6], ncomp;
bool solodinamica;
int contadorrr = 0;
float esc, densidad, volumen, nelectronesr, tempee, tempei, tempe[4], B, reso;
long long int nh20, nhp, nh2p, nelectrones, npart[4], naleprom=0, npart3[4] ;
double nale_f;
float R;
const int R_reso = 2*10*3/1+1, part_plasma_size = 6000001;                   //Primer numero es R, segundo es reso, 2 es por que es de -R a R y +1 para comenzar los arreglos desde 1
float me, mhp, mh2p, mh20, masa[4];
int carga[4];
int n_grupos[4];
///////////////////////////////////////////////////////////////////////Variables globales
int p=0, tipo, ienergia;
int rechazo;
int nceldas, n1, n2, n1i, n2i, alea_part;                                //numero de celda Nueva/Vieja Inicial/Final
char cero[] = "         0";
int rasdnd;
double dem, de_coulomb = 0, de_autoenergia = 0;
double dt_termico = 0;
long long int part_desplazadas[4] = {0};
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
    crear_archivos_iniciales_en_blanco();
    printf("\nRAND_MAX: %d\n",RAND_MAX);
    int debug=0;
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
	if(debug==1)cout << endl << "Condiciones iniciales" << endl;
	condiciones_iniciales2();
	if(debug==1)cout << endl << "Crear matriz" << endl;
	crear_matriz_plasma_rec2();
	if(debug==1)cout << endl << "Arreglo inicial" << endl;
    arreglo_inicial_part();
	if(debug==1)cout << endl << "Despues arreglo inicial" << endl;
    //mover_particulas_part();
    /*cout << "\nANTES MOVER" << endl;
    for(p=1; p<=1000000; p++){
        mover_particulas_part();
        pared();
        if(rechazo==1){
            part_plasma[alea_part][tipo] = part_plasma[0][tipo];
            importar_part_a_mallado();
        }
        if(p%1000==0){
            cout<<"\r" << p;
            actu_salida2();
        }
    }
    cout << "\nTERMINO MOVER" << endl;
    getchar();
    return(0);*/
    contador[0] = 0; contador[1] = 0; contador[2] = 0;
    contador_a[0] = 0; contador_a[1] = 0; contador_a[2] = 0;
    for(p=pasoinicial;p<=pasos;p++){
        if(debug==2)cout << "\nINICIO CICLO\n";
        //calcular_momento_total();
        //retermalizacion();
        //if((p%(100*nceldas)!=0||p<=terma)&&(p>=pasoinicial)&&(!solodinamica)){
        //if(alea()<2&&!solodinamica){
        if((dt_termico/(npart[0]+npart[1]))<5e-9&&!solodinamica){
            rechazo = 0;
            c_mov++;
            /*printf("\n<Antes mov>");
            for(int ia=1;i<=nceldas;ia++){
                if(matriz_plasma[ia].carga!=0){
                    printf("\nQue show");
                    getchar();
                }
            }*/
            //mover_particulas();
            //mover_particulas_new();
            if(debug==2)cout << "\nANTES MOVER PART\n";
            mover_particulas_part();
            part_plasma[n1][tipo].contador++;
            //printf("\nQUE SHOW 3 tipo: %d\n", tipo);getchar();
            if(debug==2)cout << "\nANTES PARED\n";
            pared();
            if(tipo!=2){
                naleprom+=nale;
                //nale_f+=1.0*nale/matriz_plasma[n2].part[tipo];
                nale_f+=(rasdnd/(1.0*RAND_MAX));
                contador_nale++;
                contadorrr++;
                if(contadorrr%100==0){
                    dat = fopen("pruebas/nale_f.dat","a");
                    fprintf(dat,"%lld %f %lld\n",nale,(rasdnd/(1.0*RAND_MAX)),matriz_plasma[n2].part[tipo]);
                    fclose(dat);
                    contadorrr = 0;
                }
            }
            if(p%1000==0){
                float ms[2] = {0}, rms[2] = {0};
                for(int i=1;i<=nceldas;i++){
                    ms[0] += pow(matriz_plasma[i].part[0],2);
                    ms[1] += pow(matriz_plasma[i].part[1],2);
                }
                rms[0] = sqrt(1.0*ms[0]/nceldas);
                rms[1] = sqrt(1.0*ms[1]/nceldas);
                dat = fopen("pruebas/rms_part.dat","a");
                fprintf(dat,"%d\t%e\t%e\n",p,rms[0],rms[1]);
                fclose(dat);
            }
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
            dinamica2();
            if(debug==2)cout << "\nANTES IMPORTAR PART A MALLADO\n";
            importar_part_a_mallado();
            dt_termico = 0;
        }
        /*if(p%actu==0){
            hacer_histograma(0,p);
            energia_cinetica(p);
            if(p%(1*actu/10)==0&&p>=1*actu/10){
                retermalizacion();
                hacer_histograma(0,p+1);
            }
        }*/
        if(p>terma){
            if(debug==1)cout << "\nANTES GDR2\n";
            gdr_plasma2();
            if(p%actu==0){
                if(debug==1)cout << "\nANTES ACTUSALIDA2\n";
                actu_salida2();
                if(debug==1)cout << "\nANTES SALIDAPROM\n";
                salida_prom2();
                for(int i=0; i<ncomp; i++)hacer_histograma2(i,p);
                //hacer_distribucion(p);
                //hacer_rms(p);
                //beta();
            }
        }
        /*long long int nn[4]={0};
        for(i=1;i<=nceldas;i++){
            nn[0] += matriz_plasma[i].part[0];
            nn[1] += matriz_plasma[i].part[1];
            nn[2] += matriz_plasma[i].part[2];
            nn[3] += matriz_plasma[i].part[3];
        }
        if(nn[0]!=nelectrones){
            printf("\nnelectrones mal\tsuma: %lld\tnelectrones: %lld\n",nn[0],nelectrones);
            getchar();
        }
        if(nn[1]!=nhp){
            printf("\nnhp mal\tsuma: %lld\tnhp: %lld\n",nn[1],nhp);
            getchar();
        }
        if(nn[2]!=nh2p){
            printf("\nnh2p mal\tsuma: %lld\tnh2p: %lld\n",nn[2],nh2p);
            getchar();
        }
        if(nn[3]!=nh20){
            printf("\nnh20 mal\tsuma: %lld\tnh20: %lld\n",nn[3],nh20);
            getchar();
        }*/
        if(p%1000==0&&p>pasoinicial||solodinamica){
            if(p%actu==0){
                printf("\rPaso: %d ",p);
                for(int i=0; i<ncomp; i++)printf("Acep%d: %1.5f ",i,1.0*contador_a[i]/contador[i]);
                printf("dt/npart: %e ",(dt_termico/(npart[0]+npart[1])));
                printf("Rneg: %d uno: %1.5f ",rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0));
                printf("Dinamica: %d Energia: %e",c_dina,energia());
                //printf("\rPaso: %i Acept. Mov: %1.5f Rechazo negativas: %1.5f Rechazo met: %1.5f Energia: %e",p,c_mova/(c_mov*1.0),rechazo_neg/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0),energia());
                /*if(solodinamica){
                    printf("\rPaso: %i Acep0: %1.5f Acep1: %1.5f Acep3: %1.5f Rneg: %d uno: %1.5f Dinamica: %i Energia: %e",p,contador_a[0]/(contador[0]*1.0),contador_a[1]/(contador[1]*1.0),contador_a[2]/(contador[2]*1.0),rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0),c_dina,energia());
                }else{
                    printf("\rPaso: %i Acep0: %1.5f Acep1: %1.5f Acep3: %1.5f Rneg: %d uno: %1.5f Dinamica: %i\tnaleprom: %lld nale_f: %f Energia: %e",p,contador_a[0]/(contador[0]*1.0),contador_a[1]/(contador[1]*1.0),contador_a[2]/(contador[2]*1.0),rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0),c_dina,naleprom/(contador_nale),1.0*nale_f/contador_nale,energia());
                }*/
            }else{
                printf("\rPaso: %d ",p);
                for(int i=0; i<ncomp; i++)printf("Acep%d: %1.5f ",i,1.0*contador_a[i]/contador[i]);
                printf("dt/npart: %e ",(dt_termico/(npart[0]+npart[1])));
                printf("Rneg: %d uno: %1.5f ",rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0));
                printf("Dinamica: %d",c_dina);
                //printf("\rPaso: %i Acept. Mov: %1.5f Rechazo negativas: %1.5f Rechazo met: %1.5f",p,c_mova/(c_mov*1.0),rechazo_neg/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0));
                /*if(solodinamica){
                    printf("\rPaso: %i Acep0: %1.5f Acep1: %1.5f Acep3: %1.5f Rneg: %d uno: %1.5f Dinamica: %i Energia: %e",p,contador_a[0]/(contador[0]*1.0),contador_a[1]/(contador[1]*1.0),contador_a[2]/(contador[2]*1.0),rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0),c_dina,energia());
                }else{
                    printf("\rPaso: %i Acep0: %1.5f Acep1: %1.5f Acep3: %1.5f Rneg: %d uno: %1.5f Dinamica: %i\tnaleprom: %lld nale_f: %f",p,contador_a[0]/(contador[0]*1.0),contador_a[1]/(contador[1]*1.0),contador_a[2]/(contador[2]*1.0),rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0),c_dina,naleprom/(contador_nale),1.0*nale_f/contador_nale);
                }*/
            }
            dat = fopen("pruebas/razon_acep.dat","a");
            fprintf(dat,"%d",p);
            for(int i=0; i<ncomp; i++)fprintf(dat,"\t%1.5f",contador_a[i]/(contador[i]*1.0));
            fprintf(dat,"\n");
            fclose(dat);
            if(!solodinamica){
                dat = fopen("pruebas/naleprom.dat","a");
                fprintf(dat,"%d\t%lld\t%f\n",p,naleprom/contador_nale,1.0*nale_f/contador_nale);
                fclose(dat);
            }
            //contador_nale=0;naleprom=0,nale_f=0;
        }
    }

    /*for(p=pasoinicial;p<=pasos;p++){
        dinamica2();
        if(p%10==0)printf("\rp: %d npart[%d]: %d npart[%d]: %d npart[%d]: %d energia: %e",p,0,n_grupos[0],1,n_grupos[1],2,n_grupos[2],energia());
        //printf("\nANTES DE IMPORTAR\n");
        importar_part_a_mallado();
        //printf("\nANTES DE ACTUSALIDA\n");
        actu_salida2();
        if(p%10==0){
            for(int i=0; i<ncomp ; i++)hacer_histograma2(i,p);
        }
        hacer_rms(p);
    }*/
    printf("LISTO");
}
void mainsillo(void){
    int i,j,k;
    srand((unsigned)time(NULL));
    dat = fopen("momento_total.dat","w");
    fprintf(dat,"#paso\tpx\tpy\n");
    fclose(dat);
    dat = fopen("energia_cinetica.dat","w");
    fprintf(dat,"#paso\tT\n");
    fclose(dat);
    dat = fopen("naleprom.dat","w");
    fclose(dat);
    dat = fopen("rms_part.dat","w");
    fclose(dat);
    dat = fopen("nale_f.dat","w");
    fclose(dat);
    /*bool a[4]={true, true, true, true};
    int x1;
    x1 = (a[0]?1:0)+(a[1]?2:0)+(a[2]?3:0)+(a[3]?4:0);
    printf("\nx1: %d\tRAND_MAX: %d",x1,RAND_MAX);
    getchar();*/

    /*dat = fopen("histograma_aleaii.dat","w");
    long int histclase[3]={0};
    long int li;
    for(li=1;li<=1e9;li++){
        if(li%10000000==0)printf("\rpaso: %ld",li);
        histclase[alea_i(-1,1)+1]++;
        //fprintf(dat,"%d\t%d\n",i,alea_i(1,15));
    }
    for(i=0;i<=2;i++){
        fprintf(dat,"%d\t%ld\n",i-1,histclase[i]);
    }
    fclose(dat);
    getchar();*/

    /*system("mkdir datos");
    system("mkdir datos/datitos");*/
    //CreateDirectory ("datos", NULL);
    //CreateDirectory ("datos/promedios", NULL);
    printf("R_reso*R_reso+1: %i\n",R_reso*R_reso+1);
    /*int histograaaama[11]={0};
    for(i=1;i<=1000000000;i++){
        if(i%100000==0)printf("\ri: %i\talea_i: %i",i,alea_i(1,10));
        histograaaama[ alea_i(1,10) ]++;
    }
    dat = fopen("histograma_aleai.txt","w");
    for(i=1;i<=10;i++){
        fprintf(dat,"%i\t%i\n",i,histograaaama[i]);
    }
    fclose(dat);*/
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
	matriz_plasma[1].vx[0] = 0; matriz_plasma[1].vy[0] = 1;
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
        //matriz_plasma[i].part[0]=matriz_plasma[i].part[3]=matriz_plasma[i].part[2]=matriz_plasma[i].part[1]=matriz_plasma[i].carga = 0;
        matriz_plasma[i].carga=matriz_plasma[i].part[0]=matriz_plasma[i].part[3]=matriz_plasma[i].part[2]=matriz_plasma[i].part[1]=0;
        matriz_plasma[i].phi=matriz_plasma[i].rho=matriz_plasma[i].ve=matriz_plasma[i].vh2p=matriz_plasma[i].vhp=matriz_plasma[i].vx[0]=matriz_plasma[i].vxh20=matriz_plasma[i].vxh2p=matriz_plasma[i].vxhp=matriz_plasma[i].vye=matriz_plasma[i].vyh20=matriz_plasma[i].vyh2p=matriz_plasma[i].vyhp=matriz_plasma[i].x=matriz_plasma[i].y=0;
        for(j=0;j<=3;j++){
            matriz_plasma[i].t[j]=0;
            matriz_plasma[i].v[j]=0;
            matriz_plasma[i].vy[j]=0;
            //matriz_plasma[i].vy_1[j]=0;
            //matriz_plasma[i].vy_1[j]=0;
            matriz_plasma[i].vx[j]=0;
            //matriz_plasma[i].vx_2[j]=0;
            //matriz_plasma[i].vx_2[j]=0;
            matriz_plasma[i].part[j]=0;
        }
    }
	condiciones_iniciales();
    printf("int(R*reso): %i\n",int(R*reso));
    printf("R_reso: %i\n",R_reso);
	//crear_matriz_plasma();
	crear_matriz_plasma_rec();
	printf("Chequeo de celdas con rho>R\n");
	for(i=1;i<=nceldas;i++){
        if(i%100==0)printf("\rCelda actual: %i",i);
        if(matriz_plasma[i].rho>R){
            imprimir_celda_plasma_rec(i);
            getchar();
        }
	}
	arreglo_inicial();
	//leer_campo_magnetico();
	p = 1 + terma;
    hacer_distribucion(0);
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
    /*for(i=1;i<=nceldas;i++){
        //matriz_plasma[nmatriz_plasma[ 62 ][ 62 ]].vx[0] = 0;
        //matriz_plasma[nmatriz_plasma[ 62 ][ 62 ]].vy[0] = 1;
        matriz_plasma[i].vx[0] = 0;
        matriz_plasma[i].vy[0] = 2*sqrt(2*kb*tempee/(pi*masa[0]));
        matriz_plasma[i].vx[1] = 0;
        matriz_plasma[i].vy[1] = 2*sqrt(2*kb*tempei/(pi*masa[1]));
        matriz_plasma[i].vx[2] = 0;
        matriz_plasma[i].vy[2] = 2*sqrt(2*kb*tempei/(pi*masa[2]));
        matriz_plasma[i].vx[3] = 0;
        matriz_plasma[i].vy[3] = 2*sqrt(2*kb*tempei/(pi*masa[3]));
        matriz_plasma[i].v[0] = norma(matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
        matriz_plasma[i].v[1] = norma(matriz_plasma[i].vx[1],matriz_plasma[i].vy[1]);
        matriz_plasma[i].v[2] = norma(matriz_plasma[i].vx[2],matriz_plasma[i].vy[2]);
        matriz_plasma[i].v[3] = norma(matriz_plasma[i].vx[3],matriz_plasma[i].vy[3]);
        if(matriz_plasma[i].vy[0]<2*sqrt(2*kb*tempee/(pi*masa[0]))){
            printf("\nvx: %f vy: %f\n",matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
            getchar();
        }
    }*/
    /*matriz_plasma[nmatriz_plasma[ 801 ][ 901 ]].vx[0] = 1;
    matriz_plasma[nmatriz_plasma[ 801 ][ 901 ]].vy[0] = 0;
    matriz_plasma[nmatriz_plasma[ 811 ][ 901 ]].vx[0] = -2;
    matriz_plasma[nmatriz_plasma[ 811 ][ 901 ]].vy[0] = 0;*/
    /*printf("INICIANDO DINAMICA\n");
    FILE *dat3, *dat4;
    sprintf(salidac,"datos/electron.dat",i);
    dat=fopen(salidac,"w");
    fprintf(dat, "#X\tY\tZ\n");
    sprintf(salidac,"datos/hp.dat",i);
    dat2=fopen(salidac,"w");
    fprintf(dat2, "#X\tY\tZ\n");
    sprintf(salidac,"datos/h2p.dat",i);
    dat3=fopen(salidac,"w");
    fprintf(dat3, "#X\tY\tZ\n");
    sprintf(salidac,"datos/h20.dat",i);
    dat4=fopen(salidac,"w");
    fprintf(dat4, "#X\tY\tZ\n");
    fclose(dat);
    fclose(dat2);
    fclose(dat3);
    fclose(dat4);
    dat=fopen("datos/electron.dat","a");
    dat2=fopen("datos/hp.dat","a");
    dat3=fopen("datos/h2p.dat","a");
    dat4=fopen("datos/h20.dat","a");
    for(i=1;i<=1000;i++){
        printf("\rdinamica_i: %i",i);
        dinamica();
        if(i%10==0){
            //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
            for(j=1;j<=nceldas;j++){
                if(matriz_plasma[j].part[0]!=0){
                    fprintf(dat,"%f\t%f\t%I64d\t%i\n", matriz_plasma[j].x, matriz_plasma[j].y, matriz_plasma[j].part[0], i);
                }
                if(matriz_plasma[j].part[1]!=0){
                    fprintf(dat2,"%f\t%f\t%I64d\t%i\n", matriz_plasma[j].x, matriz_plasma[j].y, matriz_plasma[j].part[1], i);
                }
                if(matriz_plasma[j].part[2]!=0){
                    fprintf(dat3,"%f\t%f\t%I64d\t%i\n", matriz_plasma[j].x, matriz_plasma[j].y, matriz_plasma[j].part[2], i);
                }
                if(matriz_plasma[j].part[3]!=0){
                    fprintf(dat4,"%f\t%f\t%I64d\t%i\n", matriz_plasma[j].x, matriz_plasma[j].y, matriz_plasma[j].part[3], i);
                }
            }
        }
        //getchar();
    }
    fclose(dat);
    fclose(dat2);
    fclose(dat3);
    fclose(dat4);
    printf("\nDESPUES DINAMICA\n");
    printf("Radio de giro_ electron: %f\thp: %f\th2p: %f\n",masa[0]*2*sqrt(2*kb*tempee/(pi*masa[0]))/(qe*0.01),masa[1]*2*sqrt(2*kb*tempee/(pi*masa[1]))/(qe*0.01),masa[2]*2*sqrt(2*kb*tempee/(pi*masa[2]))/(qe*0.01));
    getchar();getchar();getchar();getchar();getchar();*/
    /*for(i=1;i<=npart;i++){
        printf("\ni:%i,x:%f,y:%f,z:%f,q:%f",i,part[i].x,part[i].y,part[i].z,part[i].carga);
    }
    getchar();*/
    //terma = 0;
    //Inicio
    hacer_rms(0);
    for(p=pasoinicial;p<=pasos;p++){
        //calcular_momento_total();
        //retermalizacion();
        if((p%(100*nceldas)!=0||p<=terma)&&(p>=pasoinicial)&&(!solodinamica)){
            rechazo = 0;
            c_mov++;
            /*printf("\n<Antes mov>");
            for(int ia=1;i<=nceldas;ia++){
                if(matriz_plasma[ia].carga!=0){
                    printf("\nQue show");
                    getchar();
                }
            }*/
            mover_particulas();
            //mover_particulas_new();
            if(tipo!=3){
                naleprom+=nale;
                //nale_f+=1.0*nale/matriz_plasma[n2].part[tipo];
                nale_f+=(rasdnd/(1.0*RAND_MAX));
                contador_nale++;
                contadorrr++;
                if(contadorrr%100==0){
                    dat = fopen("nale_f.dat","a");
                    fprintf(dat,"%lld %f %lld\n",nale,(rasdnd/(1.0*RAND_MAX)),matriz_plasma[n2].part[tipo]);
                    fclose(dat);
                    contadorrr = 0;
                }
            }
            if(p%1000==0){
                float ms[2] = {0}, rms[2] = {0};
                for(i=1;i<=nceldas;i++){
                    ms[0] += pow(matriz_plasma[i].part[0],2);
                    ms[1] += pow(matriz_plasma[i].part[1],2);
                }
                rms[0] = sqrt(1.0*ms[0]/nceldas);
                rms[1] = sqrt(1.0*ms[1]/nceldas);
                dat = fopen("rms_part.dat","a");
                fprintf(dat,"%d\t%e\t%e\n",p,rms[0],rms[1]);
                fclose(dat);
            }
            if(rechazo == 0){
                metropolis_plasma();
            }else{
                rechazo_neg++;
                matriz_plasma[n1] = matriz_plasma[n1i];
                matriz_plasma[n2] = matriz_plasma[n2i];
            }
        }else{
            c_dina++;
            dinamica();
        }
        if(p%actu==0){
            hacer_histograma(0,p);
            energia_cinetica(p);
            if(p%(1*actu/10)==0&&p>=1*actu/10){
                retermalizacion();
                hacer_histograma(0,p+1);
            }
        }
        if(p>terma){
            gdr_plasma();
            if(p%actu==0){
                actu_salida();
                salida_prom();
                hacer_distribucion(p);
                hacer_rms(p);
                beta();
            }
        }
        /*long long int nn[4]={0};
        for(i=1;i<=nceldas;i++){
            nn[0] += matriz_plasma[i].part[0];
            nn[1] += matriz_plasma[i].part[1];
            nn[2] += matriz_plasma[i].part[2];
            nn[3] += matriz_plasma[i].part[3];
        }
        if(nn[0]!=nelectrones){
            printf("\nnelectrones mal\tsuma: %lld\tnelectrones: %lld\n",nn[0],nelectrones);
            getchar();
        }
        if(nn[1]!=nhp){
            printf("\nnhp mal\tsuma: %lld\tnhp: %lld\n",nn[1],nhp);
            getchar();
        }
        if(nn[2]!=nh2p){
            printf("\nnh2p mal\tsuma: %lld\tnh2p: %lld\n",nn[2],nh2p);
            getchar();
        }
        if(nn[3]!=nh20){
            printf("\nnh20 mal\tsuma: %lld\tnh20: %lld\n",nn[3],nh20);
            getchar();
        }*/
        if(p%1000==0&&p>pasoinicial||solodinamica){
            if(p%actu==0){
                //printf("\rPaso: %i Acept. Mov: %1.5f Rechazo negativas: %1.5f Rechazo met: %1.5f Energia: %e",p,c_mova/(c_mov*1.0),rechazo_neg/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0),energia());
                if(solodinamica){
                    printf("\rPaso: %i Acep0: %1.5f Acep1: %1.5f Acep3: %1.5f Rneg: %d uno: %1.5f Dinamica: %i Energia: %e",p,contador_a[0]/(contador[0]*1.0),contador_a[1]/(contador[1]*1.0),contador_a[3]/(contador[3]*1.0),rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+contador_a[3]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0),c_dina,energia());
                }else{
                    printf("\rPaso: %i Acep0: %1.5f Acep1: %1.5f Acep3: %1.5f Rneg: %d uno: %1.5f Dinamica: %i\tnaleprom: %lld nale_f: %f Energia: %e",p,contador_a[0]/(contador[0]*1.0),contador_a[1]/(contador[1]*1.0),contador_a[3]/(contador[3]*1.0),rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+contador_a[3]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0),c_dina,naleprom/(contador_nale),1.0*nale_f/contador_nale,energia());
                }
            }else{
                //printf("\rPaso: %i Acept. Mov: %1.5f Rechazo negativas: %1.5f Rechazo met: %1.5f",p,c_mova/(c_mov*1.0),rechazo_neg/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0));
                if(solodinamica){
                    printf("\rPaso: %i Acep0: %1.5f Acep1: %1.5f Acep3: %1.5f Rneg: %d uno: %1.5f Dinamica: %i Energia: %e",p,contador_a[0]/(contador[0]*1.0),contador_a[1]/(contador[1]*1.0),contador_a[3]/(contador[3]*1.0),rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+contador_a[3]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0),c_dina,energia());
                }else{
                    printf("\rPaso: %i Acep0: %1.5f Acep1: %1.5f Acep3: %1.5f Rneg: %d uno: %1.5f Dinamica: %i\tnaleprom: %lld nale_f: %f",p,contador_a[0]/(contador[0]*1.0),contador_a[1]/(contador[1]*1.0),contador_a[3]/(contador[3]*1.0),rechazo_neg,(contador_a[0]+contador_a[1]+contador_a[2]+contador_a[3]+rechazo_neg+rechazo_met_mov)/(c_mov*1.0),c_dina,naleprom/(contador_nale),1.0*nale_f/contador_nale);
                }
            }
            if(!solodinamica){
                dat = fopen("naleprom.dat","a");
                fprintf(dat,"%d\t%lld\t%f\n",p,naleprom/contador_nale,1.0*nale_f/contador_nale);
                fclose(dat);
            }
            //contador_nale=0;naleprom=0,nale_f=0;
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
    //return(  a + rand()%( b - a + 1 )   );
    return(  a + max_alea()%( b - a + 1 )   );
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
        if(B<1e-10)B=1e-10;
        fscanf(dat,"(Elec, H+, H2+, H20, todas, carga): %d, %d, %d, %d, %d, %d\n",&iprint[0],&iprint[1],&iprint[2],&iprint[3],&iprint[4],&iprint[5]);
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
////////////////////////////////////////////////////////////////////////////////////
void leer_datos_iniciales2(void){
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
            pasos = 10000;
            actu = 10;
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
////////////////////////////////////////////////////////////////////////////////////
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
    printf("(Elec, H+, H2+, H20, todas, carga): %d, %d, %d, %d, %d, %d\n",iprint[0],iprint[1],iprint[2],iprint[3],iprint[4],iprint[5]);
    printf("Solo dinamica (S_1/N_0): %d\n\n",solodinamica);
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_datos_iniciales2(void){
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
    printf("Condiciones iniciales\n");
    leer_datos_iniciales();
    esc = 1e-3;//diam*1e-9;

    densidad = 2e19;
    volumen = pi*R*esc*R*esc*1e-6;
    nelectronesr = (densidad*volumen);
    nelectrones = npart[0] = nelectronesr;
    nhp = npart[1] = nelectrones;
    nh2p = npart[2] = 0;//0.2*nelectrones+1;
    nh20 = npart[3] = 9*nelectrones;
    printf("\ndensidad: %e volumen: %e nnr: %e nelectrones: %lld",densidad,volumen,nelectronesr,nelectrones);
    printf("\nnh20: %lld nhp: %lld nh2p: %lld cargatotal: %lld\n",nh20,nhp,nh2p,nh2p+nhp-nelectrones);
    //getchar();

    //tempe[0] = tempee; tempe[1] = tempei; tempe[2] = tempei; tempe[3] = tempei;
    for(int i=0; i<ncomp; i++)tempe[i]*=qe/kb;
    tempei = tempei*qe/kb;
    tempee = tempee*qe/kb;

    for(int i=0; i<ncomp; i++)printf("t_%d: %f ",i,tempe[i]); cout << endl;

    //printf("\ntempei: %f tempee: %f\n", tempei, tempee);
    //getchar();

    imprimir_datos_iniciales();
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_iniciales2(){
    leer_datos_iniciales2();
    esc = 1e-3;

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

    imprimir_datos_iniciales2();
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
                    vi=alea_f(-4*sqrt((kb*tempee)/masa[0]), 4*sqrt((kb*tempee)/masa[0]) );
                    fvi=exp(-vi*vi*masa[0]/(2*kb*tempee))*pow(masa[0]/(2*pi*kb*tempee), 1.0/2.0);
                }while(alea()*pow(masa[0]/(2*pi*kb*tempee), 1.0/2.0)>fvi);
                matriz_plasma[contadorm].vx[0] = vi;*/
            }
        }
    }
    nceldas=contadorm;
    n1i = nceldas+1;
    n2i = nceldas+2;
    for(i=1;i<=nceldas;i++){
        do{
            vi=alea_f(-4*sqrt((kb*tempee)/masa[0]), 4*sqrt((kb*tempee)/masa[0]) );
            fvi=exp(-vi*vi*masa[0]/(2*kb*tempee))*pow(masa[0]/(2*pi*kb*tempee), 1.0/2.0);
        }while(alea()*pow(masa[0]/(2*pi*kb*tempee), 1.0/2.0)>fvi);
        matriz_plasma[i].vx[0] = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempee)/masa[0]), 4*sqrt((kb*tempee)/masa[0]) );
            fvi=exp(-vi*vi*masa[0]/(2*kb*tempee))*pow(masa[0]/(2*pi*kb*tempee), 1.0/2.0);
        }while(alea()*pow(masa[0]/(2*pi*kb*tempee), 1.0/2.0)>fvi);
        matriz_plasma[i].vy[0] = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/masa[1]), 4*sqrt((kb*tempei)/masa[1]) );
            fvi=exp(-vi*vi*masa[1]/(2*kb*tempei))*pow(masa[1]/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(masa[1]/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vx[1] = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/masa[1]), 4*sqrt((kb*tempei)/masa[1]) );
            fvi=exp(-vi*vi*masa[1]/(2*kb*tempei))*pow(masa[1]/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(masa[1]/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vy[1] = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/masa[2]), 4*sqrt((kb*tempei)/masa[2]) );
            fvi=exp(-vi*vi*masa[2]/(2*kb*tempei))*pow(masa[2]/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(masa[2]/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vx[2] = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/masa[2]), 4*sqrt((kb*tempei)/masa[2]) );
            fvi=exp(-vi*vi*masa[2]/(2*kb*tempei))*pow(masa[2]/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(masa[2]/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vy[2] = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/masa[3]), 4*sqrt((kb*tempei)/masa[3]) );
            fvi=exp(-vi*vi*masa[3]/(2*kb*tempei))*pow(masa[3]/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(masa[3]/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vx[3] = vi;
        do{
            vi=alea_f(-4*sqrt((kb*tempei)/masa[3]), 4*sqrt((kb*tempei)/masa[3]) );
            fvi=exp(-vi*vi*masa[3]/(2*kb*tempei))*pow(masa[3]/(2*pi*kb*tempei), 1.0/2.0);
        }while(alea()*pow(masa[3]/(2*pi*kb*tempei), 1.0/2.0)>fvi);
        matriz_plasma[i].vy[3] = vi;

        matriz_plasma[i].v[0] = sqrt( matriz_plasma[i].vx[0]*matriz_plasma[i].vx[0] + matriz_plasma[i].vy[0]*matriz_plasma[i].vy[0] );
        matriz_plasma[i].v[1] = sqrt( matriz_plasma[i].vx[1]*matriz_plasma[i].vx[1] + matriz_plasma[i].vy[1]*matriz_plasma[i].vy[1] );
        matriz_plasma[i].v[2] = sqrt( matriz_plasma[i].vx[2]*matriz_plasma[i].vx[2] + matriz_plasma[i].vy[2]*matriz_plasma[i].vy[2] );
    }

    printf("nceldas: %i\tn1i: %i\tn2i: %i",nceldas,n1i,n2i);
    if( dat = fopen("celdas.dat","w") ){
        fprintf(dat,"#X Y Z\n");
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%e\t%e\t%e\n",matriz_plasma[i].x,matriz_plasma[i].y, sqrt( matriz_plasma[i].x*matriz_plasma[i].x + matriz_plasma[i].y*matriz_plasma[i].y ) );
        }
        fclose(dat);
    }
    if( dat = fopen("celdas2.dat","w") ){
        for(k=1;k<=2/fraccion;k++){
            fprintf(dat,"%d\t%d\t%d\n",1,k,nmatriz_plasma[1][k]);
            for(i=1;i<=R;i++){
                for(j=1;j<=(int)((pi*i)*0.25);j++){
                    fprintf(dat,"%d\t%d\t%d\n",i,j+(k-1)*(int)((pi*i)*0.25),nmatriz_plasma[i][j+(k-1)*(int)((pi*i)*0.25)]);
                }
            }
        }
        fclose(dat);
    }

    long int clases[201]={0};
    for(i=1;i<=nceldas;i++){
        clases[ int( 200*(matriz_plasma[i].vx[0]+4*sqrt((kb*tempee)/masa[0]))/(8*sqrt((kb*tempee)/masa[0])) )+1]++;
    }
    printf("Que show");
    //getchar();
    if( dat = fopen("histograma.dat","w") ){
        fprintf(dat,"#X\tY\n");
        for(i=1;i<=200;i++){
            fprintf(dat,"%e\t%i\n", -4*sqrt((kb*tempee)/masa[0])+ (i-1)*(8*sqrt((kb*tempee)/masa[0]))/(200) , clases[i] );
        }
        fclose(dat);
    }
}
////////////////////////////////////////////////////////////////////////////////////
void crear_matriz_plasma_rec(){
    //int i,j,k,contadorm=0,contadorm2=282697,contadorvec=0;
    //int i,j,k,contadorm=0,contadorm2=11289,contadorvec=0;
    int i,j,k,contadorm=0,contadorm2=2821,contadorvec=0;
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
    printf("\nnceldas: %i\tn1i: %i\tn2i: %i\n",nceldas,n1i,n2i);

    for(i=1;i<=nceldas;i++){
        vecinos(i);
    }

    /*for(i=1;i<=nceldas;i++){
        if(i%1000==0)printf("\rCelda actual: %i",i);
        matriz_plasma[i].vx[0] = aaadist_normal(0.0,1);
        matriz_plasma[i].vy[0] = aaadist_normal(0.0,1);
        matriz_plasma[i].vx[1] = aaadist_normal(0.0,1);
        matriz_plasma[i].vy[1] = aaadist_normal(0.0,1);
        matriz_plasma[i].vx[2] = aaadist_normal(0.0,1);
        matriz_plasma[i].vy[2] = aaadist_normal(0.0,1);
        matriz_plasma[i].vx[3] = aaadist_normal(0.0,1);
        matriz_plasma[i].vy[3] = aaadist_normal(0.0,1);

        matriz_plasma[i].v[0] = sqrt( matriz_plasma[i].vx[0]*matriz_plasma[i].vx[0] + matriz_plasma[i].vy[0]*matriz_plasma[i].vy[0] );
        matriz_plasma[i].v[1] = sqrt( matriz_plasma[i].vx[1]*matriz_plasma[i].vx[1] + matriz_plasma[i].vy[1]*matriz_plasma[i].vy[1] );
        matriz_plasma[i].v[2] = sqrt( matriz_plasma[i].vx[2]*matriz_plasma[i].vx[2] + matriz_plasma[i].vy[2]*matriz_plasma[i].vy[2] );
    }*/
    printf("\rImprimiendo archivos de celda\n");
    if( dat = fopen("celdas.dat","w") ){
        fprintf(dat,"#X Y Z\n");
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%e\t%e\t%e\n",matriz_plasma[i].x,matriz_plasma[i].y, matriz_plasma[i].rho );
        }
        fclose(dat);
    }
    if( dat = fopen("celdas2.dat","w") ){
        fprintf(dat,"i\tj\tx\ty\tncelda\n");
        for(i=1;i<=2*R*reso+1;i++){
            for(j=1;j<=2*R*reso+1;j++){
                fprintf(dat,"%d\t%d\t%f\t%f\t%d\t%f\t%f\n",i,j,matriz_plasma[nmatriz_plasma[i][j]].x,matriz_plasma[nmatriz_plasma[i][j]].y, nmatriz_plasma[i][j], matriz_plasma[nmatriz_plasma[i][j]].anchow, matriz_plasma[nmatriz_plasma[i][j]].anchol);
            }
        }
        fclose(dat);
    }
    ofstream ovecinos;
    ovecinos.open("vecinos.dat");
    for(int i1 = 1; i1<=nceldas; i1++){
        ovecinos << i1 << "\t";
        for(int j1 = 1; j1<= matriz_plasma[i1].nvecinos; j1++){
            ovecinos << matriz_plasma[i1].vecinos[j1] << "\t";
        }
        ovecinos << endl;
    }
    ovecinos.close();
    cout << "\nACABO VECINOS" << endl;
    //getchar(); getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void crear_matriz_plasma_rec2(){
    //int i,j,k,contadorm=0,contadorm2=11289,contadorvec=0;
    int i,j,k,contadorm=0,contadorm2=2821,contadorvec=0;
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
    n_grupos[0] = n_grupos[1] = n_grupos[2] = 2000*nceldas;
    //n_grupos[0] = n_grupos[1] = n_grupos[2] = 200*nceldas;
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
        matriz_plasma[i].vx[0] = aaadist_normal(0.0,1);
        matriz_plasma[i].vy[0] = aaadist_normal(0.0,1);
        matriz_plasma[i].vx[1] = aaadist_normal(0.0,1);
        matriz_plasma[i].vy[1] = aaadist_normal(0.0,1);
        matriz_plasma[i].vx[2] = aaadist_normal(0.0,1);
        matriz_plasma[i].vy[2] = aaadist_normal(0.0,1);
        matriz_plasma[i].vx[3] = aaadist_normal(0.0,1);
        matriz_plasma[i].vy[3] = aaadist_normal(0.0,1);

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
    matriz_plasma[a].carga = matriz_plasma[a].part[1] + matriz_plasma[a].part[2] - matriz_plasma[a].part[0];
    /*if(matriz_plasma[a].carga!=0&&p==pasoinicial){
        printf("\nQue show! a: %d carga_i: %lld 0: %lld 1: %lld 2: %lld 3: %lld",a,carga_i,matriz_plasma[a].part[0],matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].part[3]);
        getchar();
    }*/
}
////////////////////////////////////////////////////////////////////////////////////
void arreglo_inicial(void){
    int i,j,k,ni,dummy,res=10000;
    int nuevoono=0;
    long long int cargatotal=0,ee=0, h200=0, h2pp=0, hpp=0;

    printf("\rAsignando arreglo inicial\n");
    if( dat = fopen("posiciones_inicial.dat","r") ){
        fscanf(dat,"%i\t%i\n",&pasoinicial,&termaanterior);
        fclose(dat);
        printf("\npasoinicial: %i\tterma_anterior: %i",pasoinicial,termaanterior);
    }
    //getchar();
    if((pasoinicial>0)&&(ienergia==pasoinicial)){
        printf("\nContinuando posicion final de corrida anterior");
        //getchar();
        if( dat = fopen("posiciones_inicial.dat","r") ){
            fscanf(dat,"%i\t%i\t%f\t%f\t%f\t%f\t%d\t%d",&pasoinicial,&termaanterior,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
            /*if(termaanterior>pasoinicial)terma=termaanterior-pasoinicial;
            else terma = termaanterior;*/
            terma = termaanterior;
            printf("\nterma: %d\tp: %d",terma,p);
            for(i=1;i<=nceldas;i++){
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
        printf("\nArreglo nuevo");
        ienergia = 0;
        dat = fopen("datos/energia.dat","w");
        fclose(dat);
        if(solodinamica){
            /*int contador_inicial = 0;
            float i_radio = 5.0*reso;
            for(i=-i_radio;i<=i_radio;i++){
                for(j=-i_radio;j<=i_radio;j++){
                    if(i*i+j*j<=i_radio*i_radio)contador_inicial++;
                }
            }
            for(i=-i_radio;i<=i_radio;i++){
                for(j=-i_radio;j<=i_radio;j++){
                    for(k=0;k<4;k++){
                        if(i*i+j*j<=i_radio*i_radio)matriz_plasma[ nmatriz_plasma[i+int(R*reso)+1][j+int(R*reso)+1] ].part[k] = npart[k]/contador_inicial;
                    }
                }
            }
            for(k=0;k<4;k++){
                matriz_plasma[ nmatriz_plasma[int(R*reso)+1][int(R*reso)+1] ].part[k] += npart[k]%contador_inicial;
            }*/
            matriz_plasma[nmatriz_plasma[int(floor( -10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].part[0] = npart[0];
            matriz_plasma[nmatriz_plasma[int(floor( 10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].part[1] = npart[1];
            matriz_plasma[nmatriz_plasma[int(floor( -R*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].part[3] = npart[3];
        }else{
            for(i=1;i<=nceldas;i++){
                matriz_plasma[i].part[0] = nelectrones/nceldas;
                matriz_plasma[i].part[1] = nhp/nceldas;
                matriz_plasma[i].part[2] = nh2p/nceldas;
                matriz_plasma[i].part[3] = nh20/nceldas;
            }
            for(i=1;i<=nelectrones%nceldas;i++){
                matriz_plasma[i].part[0] += 1;
            }
            for(i=1;i<=nhp%nceldas;i++){
                matriz_plasma[i].part[1] += 1;
            }
            for(i=1;i<=nh2p%nceldas;i++){
                matriz_plasma[i].part[2] += 1;
            }
            for(i=1;i<=nh20%nceldas;i++){
                matriz_plasma[i].part[3] += 1;
            }
        }
            //matriz_plasma[ nmatriz_plasma[-15+int(R*reso)+1][int(R*reso)+1] ].part[0] = npart[0];
            //matriz_plasma[ nmatriz_plasma[15+int(R*reso)+1][int(R*reso)+1] ].part[1] = npart[1];
            //matriz_plasma[ nmatriz_plasma[15+int(R*reso)+1][int(R*reso)+1] ].part[2] = npart[2];
            //matriz_plasma[ nmatriz_plasma[int(R*reso)+1][int(R*reso)+1] ].part[3] = npart[3];
        /*for(i=1;i<=nceldas;i++){
            matriz_plasma[i].part[0] = 0;
            matriz_plasma[i].part[1] = 0;
            matriz_plasma[i].part[2]= 0;
            matriz_plasma[i].part[3] = 0;
        }
        float B = 0.2;
        int celdainicial1 = nmatriz_plasma[int(floor( 10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)];
        matriz_plasma[celdainicial1].part[0] = nelectrones;
        matriz_plasma[celdainicial1].part[1] = 0;
        matriz_plasma[celdainicial1].part[2] = 0;
        matriz_plasma[celdainicial1].part[3] = nh20;
        matriz_plasma[celdainicial1].vx[0] = 0;
        matriz_plasma[celdainicial1].vy[0] = 2*sqrt(2*kb*tempee/(pi*masa[0]));

        printf("\nradio de giro 0: %f\t3: %f",(masa[0]*matriz_plasma[celdainicial1].v[0])/(carga[0]*qe*B),(masa[3]*matriz_plasma[celdainicial1].v[3])/(carga[3]*qe*B));
        printf("\nv_0: %f\tv_3: %f",matriz_plasma[celdainicial1].v[0],matriz_plasma[celdainicial1].v[3]);
        int celdainicial2 = nmatriz_plasma[int(floor( -10.1*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)];
        matriz_plasma[celdainicial2].part[0] = 0;
        matriz_plasma[celdainicial2].part[1] = nhp;
        matriz_plasma[celdainicial2].part[2] = nh2p;
        matriz_plasma[celdainicial2].part[3] = 0;
        matriz_plasma[celdainicial1].vx[1] = 0;
        matriz_plasma[celdainicial1].vy[1] = -2*sqrt(2*kb*tempee/(pi*masa[1]));
        matriz_plasma[celdainicial1].vx[2] = 0;
        matriz_plasma[celdainicial1].vy[2] = -2*sqrt(2*kb*tempee/(pi*masa[2]));

        printf("\nradio de giro 1: %f\t2: %f",(masa[1]*matriz_plasma[celdainicial2].v[1])/(carga[1]*qe*B),(masa[2]*matriz_plasma[celdainicial2].v[2])/(carga[2]*qe*B));
        printf("\nv_1: %f\tv_2: %f",matriz_plasma[celdainicial2].v[1],matriz_plasma[celdainicial2].v[2]);
        getchar();*/
        printf("\nCalculando carga");
        for(i=1;i<=nceldas;i++){
            calc_carga(i);
            if(matriz_plasma[i].carga!=0){
                printf("\nQue show i: %d carga: %lld",i,matriz_plasma[i].carga );
                getchar();
            }
        }
    }
    printf("\nAsignando velocidades\n");
    for(i=1;i<=nceldas;i++){
        for(j=0;j<4;j++){
            matriz_plasma[i].vx[j]=0;
            matriz_plasma[i].vy[j]=0;
            matriz_plasma[i].v[j]=0;
        }
    }
    /*matriz_plasma[nmatriz_plasma[int(floor( -10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].vy[0] = -sqrt(2*kb*tempe[0]/(masa[0]));
    matriz_plasma[nmatriz_plasma[int(floor( -10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].v[0] = sqrt(2*kb*tempe[0]/(masa[0]));
    matriz_plasma[nmatriz_plasma[int(floor( -10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].vx[0] = 0;
    matriz_plasma[nmatriz_plasma[int(floor( 10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].vy[1] = -sqrt(2*kb*tempe[1]/(masa[1]));
    matriz_plasma[nmatriz_plasma[int(floor( 10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].v[1] = sqrt(2*kb*tempe[1]/(masa[1]));
    matriz_plasma[nmatriz_plasma[int(floor( 10.0*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].vx[1] = 0;
    matriz_plasma[nmatriz_plasma[int(floor( R*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].vy[3] = matriz_plasma[nmatriz_plasma[int(floor(R*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].v[3] = sqrt(2*kb*tempe[3]/(masa[3]));
    matriz_plasma[nmatriz_plasma[int(floor( R*reso) + floor(R*reso) + 1)][int(floor( 0.0*reso) + floor(R*reso) + 1)]].vx[3] = 0;*/
    for(i=1;i<=nceldas;i++){
        if(i%1000==0||i==nceldas)printf("\rCelda actual: %i",i);
        if(matriz_plasma[i].part[0]!=0){
            matriz_plasma[i].vx[0] = aaadist_normal(0.0,sqrt(kb*tempee/masa[0]));
            //matriz_plasma[i].vx[0] = 0;
            //if(matriz_plasma[i].vx[0]>=0) matriz_plasma[i].vx_1[0] = matriz_plasma[i].vx[0];
            //else matriz_plasma[i].vx_2[0] = matriz_plasma[i].vx[0];
            matriz_plasma[i].vy[0] = aaadist_normal(0.0,sqrt(kb*tempee/masa[0]));
            //matriz_plasma[i].vy[0] = -sqrt(kb*tempee/masa[0]);
            //if(matriz_plasma[i].vy[0]>=0) matriz_plasma[i].vy_1[0] = matriz_plasma[i].vy[0];
            //else matriz_plasma[i].vy_2[0] = matriz_plasma[i].vy[0];
            matriz_plasma[i].v[0] = norma(matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
        }
        if(matriz_plasma[i].part[1]!=0){
            matriz_plasma[i].vx[1] = aaadist_normal(0.0,sqrt(kb*tempei/masa[1]));
            //matriz_plasma[i].vx[1] = 0;
            //if(matriz_plasma[i].vx[1]>=0) matriz_plasma[i].vx_1[1] = matriz_plasma[i].vx[1];
            //else matriz_plasma[i].vx_2[1] = matriz_plasma[i].vx[1];
            matriz_plasma[i].vy[1] = aaadist_normal(0.0,sqrt(kb*tempei/masa[1]));
            //matriz_plasma[i].vy[1] = -sqrt(kb*tempei/masa[1]);
            //if(matriz_plasma[i].vy[1]>=0) matriz_plasma[i].vy_1[1] = matriz_plasma[i].vy[1];
            //else matriz_plasma[i].vy_2[1] = matriz_plasma[i].vy[1];
            matriz_plasma[i].v[1] = norma(matriz_plasma[i].vx[1],matriz_plasma[i].vy[1]);
        }
        if(matriz_plasma[i].part[2]!=0){
            matriz_plasma[i].vx[2] = aaadist_normal(0.0,sqrt(kb*tempei/masa[2]));
            //matriz_plasma[i].vx[2] = 0;
            //if(matriz_plasma[i].vx[2]>=0) matriz_plasma[i].vx_1[2] = matriz_plasma[i].vx[2];
            //else matriz_plasma[i].vx_2[2] = matriz_plasma[i].vx[2];
            matriz_plasma[i].vy[2] = aaadist_normal(0.0,sqrt(kb*tempei/masa[2]));
            //matriz_plasma[i].vy[2] = -sqrt(kb*tempei/masa[2]);
            //if(matriz_plasma[i].vy[2]>=0) matriz_plasma[i].vy_1[2] = matriz_plasma[i].vy[2];
            //else matriz_plasma[i].vy_2[2] = matriz_plasma[i].vy[2];
            matriz_plasma[i].v[2] = norma(matriz_plasma[i].vx[2],matriz_plasma[i].vy[2]);
        }
        if(matriz_plasma[i].part[3]!=0){
            matriz_plasma[i].vx[3] = aaadist_normal(0.0,sqrt(kb*tempei/masa[3]));
            //matriz_plasma[i].vx[3] = sqrt(kb*tempei/masa[3]);
            //if(matriz_plasma[i].vx[3]>=0) matriz_plasma[i].vx_1[3] = matriz_plasma[i].vx[3];
            //else matriz_plasma[i].vx_2[3] = matriz_plasma[i].vx[3];
            matriz_plasma[i].vy[3] = aaadist_normal(0.0,sqrt(kb*tempei/masa[3]));
            //matriz_plasma[i].vy[3] = 0;
            //if(matriz_plasma[i].vy[3]>=0) matriz_plasma[i].vy_1[3] = matriz_plasma[i].vy[3];
            //else matriz_plasma[i].vy_2[3] = matriz_plasma[i].vy[3];
            matriz_plasma[i].v[3] = norma(matriz_plasma[i].vx[3],matriz_plasma[i].vy[3]);
        }
    }
    printf("\nAntes de histograma de velocidades.");
    hacer_histograma(0,0);
    hacer_histograma(1,0);
    hacer_histograma(2,0);
    hacer_histograma(3,0);
    energia_cinetica(0);
    printf("\nDespues de histograma de velocidades.");
    retermalizacion();
    hacer_histograma(0,1);
    /*long int clases[201]={0};
    for(i=1;i<=nceldas;i++){
        clases[ int( 200*(matriz_plasma[i].vx[0]+5*sqrt((kb*tempee)/masa[0]))/(10*sqrt((kb*tempee)/masa[0])) )+1]++;
    }
    //getchar();
    if( dat = fopen("histograma.dat","w") ){
        printf("\nHaciendo histograma de velocidades");
        fprintf(dat,"#X\tY\n");
        for(i=1;i<=200;i++){
            fprintf(dat,"%e\t%i\n", -5*sqrt((kb*tempee)/masa[0])+ (i-1)*(10*sqrt((kb*tempee)/masa[0]))/(200) , clases[i] );
        }
        fclose(dat);
        printf("\nTermino de hacer histograma de velocidades");
    }*/

    /*if(matriz_plasma[1].part[0]==0){
        printf("\nAsignando electrones");
        for(i=1;i<=(int)(nelectrones/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].part[0]+=res;
        }
        printf("\nAsignando moleculas de hidrogeno");
        for(i=1;i<=(int)(nh20/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].part[3]+=res;
        }
        printf("\nAsignando hidrones");
        for(i=1;i<=(int)(nhp/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].part[1]+=res;
        }
        printf("\nAsignando deuterones");
        for(i=1;i<=(int)(nh2p/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].part[2]+=res;
        }
        ni = (int)((alea())*nceldas)+1;
        if(ni==nceldas+1)ni=nceldas;
        matriz_plasma[ni].part[0]+=nelectrones%res;
        ni = (int)((alea())*nceldas)+1;
        if(ni==nceldas+1)ni=nceldas;
        matriz_plasma[ni].part[3]+=nh20%res;
        ni = (int)((alea())*nceldas)+1;
        if(ni==nceldas+1)ni=nceldas;
        matriz_plasma[ni].part[1]+=nhp%res;
        ni = (int)((alea())*nceldas)+1;
        if(ni==nceldas+1)ni=nceldas;
        matriz_plasma[ni].part[2]+=nh2p%res;

        for(i=1;i<=nceldas;i++){
            calc_carga(i);
        }
    }*/
    /*for(i=1;i<=nceldas;i++){
        if(i==nmatriz_plasma[ int((R+20)*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].part[0]=nelectrones/2.0;
            matriz_plasma[i].part[2]=0;
            matriz_plasma[i].part[1]=0;
            matriz_plasma[i].part[3]=0;
        }
        else if(i==nmatriz_plasma[ int((R+15)*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].part[0]=nelectrones/2.0;
            matriz_plasma[i].part[2]=0;
            matriz_plasma[i].part[1]=0;
            matriz_plasma[i].part[3]=0;
        }else if(i==nmatriz_plasma[ int((R-10)*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].part[0]=0;
            matriz_plasma[i].part[2]=0;
            matriz_plasma[i].part[1]=nhp/2.0;
            matriz_plasma[i].part[3]=nh20/2.0;
        }else if(i==nmatriz_plasma[ int((R-5)*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].part[0]=0;
            matriz_plasma[i].part[2]=0;
            matriz_plasma[i].part[1]=nhp/2.0;
            matriz_plasma[i].part[3]=nh20/2.0;
        }else if(i==nmatriz_plasma[ int(R*reso) + 1 ][ int(R*reso) + 1 ]){
            matriz_plasma[i].part[0]=0;
            matriz_plasma[i].part[2]=nh2p;
            matriz_plasma[i].part[1]=0;
            matriz_plasma[i].part[3]=0;
        }else{
            matriz_plasma[i].part[0]=0;
            matriz_plasma[i].part[2]=0;
            matriz_plasma[i].part[1]=0;
            matriz_plasma[i].part[3]=0;
        }
        calc_carga(i);
    }*/
    /*if((matriz_plasma[1].part[0]==0)&&(matriz_plasma[1].part[3]==0)&&(matriz_plasma[1].part[1]==0)&&(matriz_plasma[1].part[2]==0)&&(matriz_plasma[1].carga==0)){
        for(i=1;i<=5;i++){
            for(j=1;j<=2;j++){
                matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[0] = nelectrones/10;
                matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[2] = nh2p/10;
                matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[1] = nhp/10;
                matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[3] = nh20/10;
                cargatotal += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[2] + matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[1] - matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[0];
                printf("\ni: %i j: %i electrones: %I64d h20: %I64d hp: %I64d h2p: %I64d carga: %I64d",i,j,nelectrones/10,nh20/10,nhp/10,nh2p/10,-nelectrones/10+nh2p/10+nhp/10);
                ee += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[0];
                h200 += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[3];
                h2pp += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[2];
                hpp += matriz_plasma[ nmatriz_plasma[30+i][15+j] ].part[1];
                calc_carga( nmatriz_plasma[30+i][15+j] );
            }
        }
        matriz_plasma[500].part[0] = 4;
        matriz_plasma[500].part[3] = 6;
        matriz_plasma[500].part[1] = 5;
        matriz_plasma[500].part[2] = 9;
        calc_carga(500);
        cargatotal += matriz_plasma[500].part[1] + matriz_plasma[500].part[2] - matriz_plasma[500].part[0];
        printf("\nCargatotal: %I64d ee: %I64d h200: %I64d hpp: %I64d h2pp: %I64d",cargatotal,ee-nelectrones,h200-nh20,hpp-nhp,h2pp-nh2p);
    }*/
    sprintf(salidac,"datos/posiciones_2.dat");
    if( dat=fopen(salidac,"w") ){
        for(i=1;i<=nceldas;i++){
            printposiciones(i);
        }
        fclose(dat);
    }
    //printf("\nYa termino posiciones2");
    //getchar();
    p = terma;
    actu_salida();
    p = 0;
}
////////////////////////////////////////////////////////////////////////////////////
void arreglo_inicial_part(void){
    float xrand, yrand, dummy, R_inicial = R;
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
            part_cont=0;
            do{
                xrand = alea_f(-R_inicial,R_inicial);
                yrand = alea_f(-R_inicial,R_inicial);
                if(xrand*xrand+yrand*yrand<R_inicial*R_inicial){
                    if((10*part_cont)%n_grupos[j]==0)printf("\rj: %d part_cont: %d",j,part_cont);
                    part_cont++;
                    part_plasma[part_cont][j].x = xrand;
                    part_plasma[part_cont][j].y = yrand;
                    part_plasma[part_cont][j].vx = aaadist_normal(0.0,sqrt(kb*tempe[j]/masa[j]));
                    part_plasma[part_cont][j].vy = aaadist_normal(0.0,sqrt(kb*tempe[j]/masa[j]));
                    //part_plasma[part_cont][j].vx = alea_i(-1,1)*2*sqrt(kb*tempe[j]/masa[j]);
                    //part_plasma[part_cont][j].vy = alea_i(-1,1)*2*sqrt(kb*tempe[j]/masa[j]);
                    part_plasma[part_cont][j].v = norma(part_plasma[part_cont][j].vx,part_plasma[part_cont][j].vy);
                    part_plasma[part_cont][j].part = npart3[j];
                    if(j==3){
                        printf("\nnpart_%d: %lld part_%d.part: %lld",j,npart3[j],j,part_plasma[part_cont][j].part);
                        getchar();
                    }
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
        cout << "Despues de part" << endl;
        importar_part_a_mallado();
        cout << "Despues de importar" << endl;
        p = 0;
        actu_salida2();
        for(int j=0; j<ncomp; j++)hacer_histograma2(j,0);
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
            icelda = part_a_mallado(i,j);
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
int part_a_mallado(int a, int b){
    int ix, iy;
    ix = floor( part_plasma[a][b].x*reso ) + floor(R*reso) + 1;
    iy = floor( part_plasma[a][b].y*reso ) + floor(R*reso) + 1;
    return(nmatriz_plasma[ix][iy]);
}
////////////////////////////////////////////////////////////////////////////////////
void hacer_histograma(int a, int b){
    const int int_clases = 200;
    long long int clasesx[int_clases+1]={0},clasesy[int_clases+1]={0};
    int i;
    char nombre[50];
    float vxmin = matriz_plasma[1].vx[a], vxmax = matriz_plasma[1].vx[a], vymin = matriz_plasma[1].vy[a], vymax = matriz_plasma[1].vy[a];
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
        if(matriz_plasma[i].vx[a]<vxmin)vxmin=matriz_plasma[i].vx[a];
        if(matriz_plasma[i].vx[a]>vxmax)vxmax=matriz_plasma[i].vx[a];
        if(matriz_plasma[i].vy[a]<vymin)vymin=matriz_plasma[i].vy[a];
        if(matriz_plasma[i].vy[a]>vymax)vymax=matriz_plasma[i].vy[a];
        //fprintf(dat,"%d\t%e\n",i,matriz_plasma[i].vx[a]);
        //fprintf(dat2,"%d\t%e\n",i,matriz_plasma[i].vy[a]);
    }
    vxmin-=1;vxmax+=1;
    vymin-=1;vymax+=1;
    //fclose(dat);fclose(dat2);
    for(i=1;i<=nceldas;i++){
            if(int( int_clases*(matriz_plasma[i].vx[a]-vxmin)/(vxmax-vxmin) )+1<0){
                printf("\nVXX: %e\tvxmin: %e\tvxmax: %e",matriz_plasma[i].vx[a],vxmin,vxmax);
                getchar();
            }
            if(int( int_clases*(matriz_plasma[i].vy[a]-vymin)/(vymax-vymin) )+1<0){
                printf("\nVYY: %e\tvymin: %e\tvymax: %e",matriz_plasma[i].vy[a],vymin,vymax);
                getchar();
            }
        clasesx[ int( int_clases*(matriz_plasma[i].vx[a]-vxmin)/(vxmax-vxmin) )+1]+=matriz_plasma[i].part[a];
        clasesy[ int( int_clases*(matriz_plasma[i].vy[a]-vymin)/(vymax-vymin) )+1]+=matriz_plasma[i].part[a];
    }
    sprintf(nombre,"histogramas/histv_%d_%d.dat",a,b);
    if( dat = fopen(nombre,"w") ){
        fprintf(dat,"#X\tY\n");
        for(i=1;i<=int_clases;i++){
            /*if(clases[i]>300){
                printf("\nQue show! clases[%d]: %d",i,clases[i]);
                getchar();
            }*/
            fprintf(dat,"%e\t%e\t%e\t%e\n", vxmin + (i-1)*(vxmax-vxmin)/int_clases , 1.0*clasesx[i]/npart[a], vymin + (i-1)*(vymax-vymin)/int_clases , 1.0*clasesy[i]/npart[a] );
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
void hacer_histograma2(int a, int b){
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
void vel_txt(char a[10],int b){
    FILE *daaat;
    sprintf(a,"vel%d.txt",b);
    printf("\nImprimiendo: %s",a);
    daaat = fopen(a,"w");
    for(int i=1;i<=nceldas;i++){
        if(matriz_plasma[i].part[0]!=0||matriz_plasma[i].part[1]!=0||matriz_plasma[i].part[2]!=0||matriz_plasma[i].part[3]!=0)fprintf(daaat,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i,matriz_plasma[i].vx[0],matriz_plasma[i].vy[0],matriz_plasma[i].vx[1],matriz_plasma[i].vy[1],matriz_plasma[i].vx[2],matriz_plasma[i].vy[2],matriz_plasma[i].vx[3],matriz_plasma[i].vy[3]);
    }
    fclose(daaat);
    getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void vel_txt2(char a[10],int b){
    FILE *daaat;
    sprintf(a,"vel%d.txt",b);
    printf("\nImprimiendo: %s",a);
    daaat = fopen(a,"w");
    for(int i=1;i<=nceldas;i++){
        if(matriz_plasma2[i].part[0]!=0||matriz_plasma2[i].part[1]!=0||matriz_plasma2[i].part[2]!=0||matriz_plasma2[i].part[3]!=0)fprintf(daaat,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i,matriz_plasma2[i].vx[0],matriz_plasma2[i].vy[0],matriz_plasma2[i].vx[1],matriz_plasma2[i].vy[1],matriz_plasma2[i].vx[2],matriz_plasma2[i].vy[2],matriz_plasma2[i].vx[3],matriz_plasma2[i].vy[3]);
    }
    fclose(daaat);
    getchar();
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
////////////////////////////////////////////////////////////////////////////////////ELIGE EN TODO EL ESPACIO
void mover_particulas(void){
    float nalea;
    do{
        tale = alea();
        if(tale<0.25&&nelectrones!=0){
            tipo = 0;
        }else if(tale<0.5&&tale>=0.25&&nhp!=0){
            tipo = 1;
        }else if(tale<0.75&&tale>=0.5&&nh2p!=0){
            tipo = 2;
        }else if(tale>=0.75&&nh20!=0){
            tipo = 3;
        }else{
            tipo = 4;
        }
    }while(tipo==4);
    n1 = int(nceldas*alea())+1;
    if(n1==nceldas+1)n1=nceldas;
    n2 = int(nceldas*alea())+1;
    if(n2==nceldas+1)n2=nceldas;

    while(n2==n1||matriz_plasma[n2].part[tipo]==0){
        n2 = int(nceldas*alea())+1;
        if(n2==nceldas+1)n2=nceldas;
    }
    matriz_plasma[n1i] = matriz_plasma[n1];
    matriz_plasma[n2i] = matriz_plasma[n2];

    //if(matriz_plasma[n1i].carga!=0||matriz_plasma[n2i].carga!=0)printf("\npaso: %d n1: %d n2: %d n1i.q: %lld n2i.q: %lld",p,n1,n2,matriz_plasma[n1i].carga,matriz_plasma[n2i].carga);

    contador[tipo]++;

    /*int inalea;
    printf("\n");
    for(int _i=1;_i<=1000000;_i++){
        inalea = rand();
        nalea = 1.0*inalea/RAND_MAX;//alea();
        nale = nalea*matriz_plasma[n2].part[tipo];
        nale += 1;
        nale = inalea*matriz_plasma[n2].part[tipo]/RAND_MAX;
        printf("\ri: %d",_i);
        if(nale>matriz_plasma[n2].part[tipo]){
            printf("\nQue show! :v");
            getchar();
        }
        //if(nale!=inalea*matriz_plasma[n2].part[tipo]/RAND_MAX+1){
            //printf("\nnale: %I64d\tnalea: %f\tn2_part[%d]: %I64d\totronumero[%d]: %I64d",nale,nalea,tipo,matriz_plasma[n2].part[tipo],tipo,inalea*matriz_plasma[n2].part[tipo]/RAND_MAX+1);
            //getchar();
        //}
    }*/
    rasdnd = rand();
    nale = (rasdnd/(1.0*RAND_MAX))*matriz_plasma[n2].part[tipo];
    if(nale>matriz_plasma[n2].part[tipo])nale=matriz_plasma[n2].part[tipo];
    //if(matriz_plasma[n1i].part[tipo]+nale<1){
    //    printf("\ndividiendo entre cero!\nn1_part[tipo]: %I64d\tnale: %I64d",matriz_plasma[n1i].part[tipo],nale);
    //    getchar();
    //}
    //matriz_plasma[n1].vx_1[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vx_1[tipo]+nale*matriz_plasma[n2].vx_1[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
    //matriz_plasma[n1].vy_1[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vy_1[tipo]+nale*matriz_plasma[n2].vy_1[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
    //matriz_plasma[n1].vx_2[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vx_2[tipo]+nale*matriz_plasma[n2].vx_2[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
    //matriz_plasma[n1].vy_2[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vy_2[tipo]+nale*matriz_plasma[n2].vy_2[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
    //matriz_plasma[n1].vx[tipo] = matriz_plasma[n1].vx_1[tipo]+matriz_plasma[n1].vx_2[tipo];
    //matriz_plasma[n1].vy[tipo] = matriz_plasma[n1].vy_1[tipo]+matriz_plasma[n1].vy_2[tipo];
    if(nale!=0||matriz_plasma[n1i].part[tipo]!=0){
        matriz_plasma[n1].vx[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vx[tipo]+nale*matriz_plasma[n2].vx[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
        matriz_plasma[n1].vy[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vy[tipo]+nale*matriz_plasma[n2].vy[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
        matriz_plasma[n1].v[tipo] = norma(matriz_plasma[n1].vx[tipo],matriz_plasma[n1].vy[tipo]);
    }else{
        matriz_plasma[n1].vx[tipo] = 0;
        matriz_plasma[n1].vy[tipo] = 0;
        matriz_plasma[n1].v[tipo] = norma(matriz_plasma[n1].vx[tipo],matriz_plasma[n1].vy[tipo]);
    }
    matriz_plasma[n1].part[tipo] +=  nale;
    matriz_plasma[n2].part[tipo] -=  nale;

    if((matriz_plasma[n2].part[tipo]<0)||(matriz_plasma[n1].part[tipo]<0)){
        rechazo = 1;
        contador_rechazo++;
    }
    calc_carga(n1);
    calc_carga(n2);
    //if(matriz_plasma[n1i].carga!=0||matriz_plasma[n2i].carga!=0)printf(" n1.q: %lld n2.q: %lld",matriz_plasma[n1].carga,matriz_plasma[n2].carga);
    //if(matriz_plasma[n1i].carga!=0||matriz_plasma[n2i].carga!=0)getchar();
}
////////////////////////////////////////////////////////////////////////////////////ELIGE ALREDEDOR DE LA CELDA N2
void mover_particulas_new(void){
    float nalea;
    int i, j, int_x, int_y;
    do{
        tale = alea();
        if(tale<0.25&&nelectrones!=0){
            tipo = 0;
        }else if(tale<0.5&&tale>=0.25&&nhp!=0){
            tipo = 1;
        }else if(tale<0.75&&tale>=0.5&&nh2p!=0){
            tipo = 2;
        }else if(tale>=0.75&&nh20!=0){
            tipo = 3;
        }else{
            tipo = 4;
        }
    }while(tipo==4);
    //n1 = int(nceldas*alea())+1;
    //if(n1==nceldas+1)n1=nceldas;
    do{
        n2 = int(nceldas*alea())+1;
        if(n2==nceldas+1)n2=nceldas;
    }while(matriz_plasma[n2].part[tipo]==0);
    //n2 = 1;
    do{
        int_x = matriz_plasma[n2].xi + int(R*reso) + 1 + (alea_i(-1,1));
        int_y = matriz_plasma[n2].yi + int(R*reso) + 1 + (alea_i(-1,1));
    }while(nmatriz_plasma[int_x][int_y]>nceldas||int_x<1||int_x>R_reso||int_y<1||int_y>R_reso||nmatriz_plasma[int_x][int_y]==n2);
    n1 = nmatriz_plasma[int_x][int_y];
    //printf("\nn1: %d\tn2: %d",n1,n2);
    //getchar();
    //elegir_n1();

    //if(p<1){
        matriz_plasma[n1i] = matriz_plasma[n1];
        matriz_plasma[n2i] = matriz_plasma[n2];
        //if(tipo!=3)printf("\npart[%d]: %lld",tipo,matriz_plasma[n2i].part[tipo]);

        /*printf("\nn1: %i n2: %i nale: %i",n1,n2,nale);
        imprimir_celda_plasma(n1);
        imprimir_celda_plasma(n2);*/

        contador[tipo]++;
        /*nalea = alea();
        nale = nalea*matriz_plasma[n2].part[tipo];
        nale += 1;*/
        rasdnd = rand();
        nale = (1.0*rasdnd/RAND_MAX)*matriz_plasma[n2].part[tipo];
        if(nale>matriz_plasma[n2].part[tipo])nale=matriz_plasma[n2].part[tipo];
        /*if(matriz_plasma[n1i].part[tipo]+nale<1){
            printf("\ndividiendo entre cero!\nn1_part[tipo]: %I64d\tnale: %I64d",matriz_plasma[n1i].part[tipo],nale);
            getchar();
        }*/
        /*matriz_plasma[n1].vx_1[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vx_1[tipo]+nale*matriz_plasma[n2].vx_1[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
        matriz_plasma[n1].vy_1[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vy_1[tipo]+nale*matriz_plasma[n2].vy_1[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
        matriz_plasma[n1].vx_2[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vx_2[tipo]+nale*matriz_plasma[n2].vx_2[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
        matriz_plasma[n1].vy_2[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vy_2[tipo]+nale*matriz_plasma[n2].vy_2[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
        matriz_plasma[n1].vx[tipo] = matriz_plasma[n1].vx_1[tipo]+matriz_plasma[n1].vx_2[tipo];
        matriz_plasma[n1].vy[tipo] = matriz_plasma[n1].vy_1[tipo]+matriz_plasma[n1].vy_2[tipo];*/
        /*matriz_plasma[n1].vx[tipo] = ((matriz_plasma[n1i].part[tipo]-nale)*matriz_plasma[n1i].vx[tipo]+2*nale*matriz_plasma[n2].vx[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
        matriz_plasma[n1].vy[tipo] = ((matriz_plasma[n1i].part[tipo]-nale)*matriz_plasma[n1i].vy[tipo]+2*nale*matriz_plasma[n2].vy[tipo])/(matriz_plasma[n1i].part[tipo]+nale);*/
        if(nale!=0||matriz_plasma[n1i].part[tipo]!=0){
            matriz_plasma[n1].vx[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vx[tipo]+nale*matriz_plasma[n2].vx[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
            matriz_plasma[n1].vy[tipo] = (matriz_plasma[n1i].part[tipo]*matriz_plasma[n1i].vy[tipo]+nale*matriz_plasma[n2].vy[tipo])/(matriz_plasma[n1i].part[tipo]+nale);
            matriz_plasma[n1].v[tipo] = norma(matriz_plasma[n1].vx[tipo],matriz_plasma[n1].vy[tipo]);
        }
        else{
            matriz_plasma[n1].vx[tipo] = 0;
            matriz_plasma[n1].vy[tipo] = 0;
            matriz_plasma[n1].v[tipo] = norma(matriz_plasma[n1].vx[tipo],matriz_plasma[n1].vy[tipo]);
        }
        matriz_plasma[n1].part[tipo] +=  nale;
        matriz_plasma[n2].part[tipo] -=  nale;

        if((matriz_plasma[n2].part[tipo]<0)||(matriz_plasma[n1].part[tipo]<0)){
            rechazo = 1;
            contador_rechazo++;
        }

        calc_carga(n1);
        calc_carga(n2);

        /*if(tale<0.25){
            tipo = 1;
            contador[0]++;
            nalea = alea();
            nale = nalea*matriz_plasma[n2].part[0];
            nale += 1;
            if(nale==matriz_plasma[n2].part[0]+1)nale=matriz_plasma[n2].part[0];
            //nale = alea_i(1,matriz_plasma[n2].part[0]);
            if(matriz_plasma[n1i].part[0]+nale<1){
                printf("\ndividiendo entre cero!");
                getchar();
            }
            matriz_plasma[n1].vx[0] = (matriz_plasma[n1i].part[0]*matriz_plasma[n1i].vx[0]+nale*matriz_plasma[n2].vx[0])/(matriz_plasma[n1i].part[0]+nale);
            matriz_plasma[n1].vy[0] = (matriz_plasma[n1i].part[0]*matriz_plasma[n1i].vy[0]+nale*matriz_plasma[n2].vy[0])/(matriz_plasma[n1i].part[0]+nale);
            matriz_plasma[n1].v[0] = norma(matriz_plasma[n1].vx[0],matriz_plasma[n1].vy[0]);
            matriz_plasma[n1].part[0] +=  nale;
            matriz_plasma[n2].part[0] -=  nale;
            if((matriz_plasma[n2].part[0]<0)||(matriz_plasma[n1].part[0]<0)) rechazo=1;
        }
        else if(tale<0.5){
            tipo = 2;
            contador[1]++;
            nalea = alea();
            nale = nalea*matriz_plasma[n2].part[1];
            nale += 1;
            if(nale==matriz_plasma[n2].part[1]+1)nale=matriz_plasma[n2].part[1];
            //nale = alea_i(1,matriz_plasma[n2].part[1]);
            if(matriz_plasma[n1i].part[1]+nale<1){
                printf("\ntipo: %i  dividiendo entre cero!",tipo);
                getchar();
            }
            matriz_plasma[n1].vx[1] = (matriz_plasma[n1i].part[1]*matriz_plasma[n1i].vx[1]+nale*matriz_plasma[n2].vx[1])/(matriz_plasma[n1i].part[1]+nale);
            matriz_plasma[n1].vy[1] = (matriz_plasma[n1i].part[1]*matriz_plasma[n1i].vy[1]+nale*matriz_plasma[n2].vy[1])/(matriz_plasma[n1i].part[1]+nale);
            matriz_plasma[n1].v[1] = norma(matriz_plasma[n1].vx[1],matriz_plasma[n1].vy[1]);
            matriz_plasma[n1].part[1] +=  nale;
            matriz_plasma[n2].part[1] -=  nale;
            if((matriz_plasma[n2].part[1]<0)||(matriz_plasma[n1].part[1]<0)) rechazo=1;
        }
        else if(tale<0.75){
            tipo = 2;
            contador[1]++;
            nalea = alea();
            nale = nalea*matriz_plasma[n2].part[2];
            nale += 1;
            if(nale==matriz_plasma[n2].part[2]+1)nale=matriz_plasma[n2].part[2];
            //nale = alea_i(1,matriz_plasma[n2].part[2]);
            if(matriz_plasma[n1i].part[2]+nale<1){
                printf("\ntipo: %i  dividiendo entre cero!",tipo);
                getchar();
            }
            matriz_plasma[n1].vx[2] = (matriz_plasma[n1i].part[2]*matriz_plasma[n1i].vx[2]+nale*matriz_plasma[n2].vx[2])/(matriz_plasma[n1i].part[2]+nale);
            matriz_plasma[n1].vy[2] = (matriz_plasma[n1i].part[2]*matriz_plasma[n1i].vy[2]+nale*matriz_plasma[n2].vy[2])/(matriz_plasma[n1i].part[2]+nale);
            matriz_plasma[n1].v[2] = norma(matriz_plasma[n1].vx[2],matriz_plasma[n1].vy[2]);
            matriz_plasma[n1].part[2] +=  nale;
            matriz_plasma[n2].part[2] -=  nale;
            if((matriz_plasma[n2].part[2]<0)||(matriz_plasma[n1].part[2]<0)) rechazo=1;
        }
        else{
            tipo = 3;
            contador[2]++;
            nalea = alea();
            nale = nalea*matriz_plasma[n2].part[3];
            nale += 1;
            if(nale==matriz_plasma[n2].part[3]+1)nale=matriz_plasma[n2].part[3];
            //nale = alea_i(1,matriz_plasma[n2].part[3]);
            if(matriz_plasma[n1i].part[3]+nale<1){
                printf("\ntipo: %i  dividiendo entre cero!",tipo);
                getchar();
            }
            matriz_plasma[n1].vx[3] = (matriz_plasma[n1i].part[3]*matriz_plasma[n1i].vx[3]+nale*matriz_plasma[n2].vx[3])/(matriz_plasma[n1i].part[3]+nale);
            matriz_plasma[n1].vy[3] = (matriz_plasma[n1i].part[3]*matriz_plasma[n1i].vy[3]+nale*matriz_plasma[n2].vy[3])/(matriz_plasma[n1i].part[3]+nale);
            matriz_plasma[n1].v[3] = norma(matriz_plasma[n1].vx[3],matriz_plasma[n1].vy[3]);
            matriz_plasma[n1].part[3] +=  nale;
            matriz_plasma[n2].part[3] -=  nale;
            if((matriz_plasma[n2].part[3]<0)||(matriz_plasma[n1].part[3]<0)) rechazo=1;
        }
        calc_carga(n1);
        calc_carga(n2);*/
        /*if(rechazo==1){
            printf("tipo: %i\tnalea: %f\tnale: %I64d",tipo,nalea,nale);
            getchar();
        }*/
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
void mover_particulas_part(void){
    int idx, idy;
    int c_pared = 0;
    tipo = alea_i( 0, ncomp-2 );
    //printf("\ntipo: %d contador_tipo: %d",tipo,contador[tipo]);
    contador[tipo]++;
    //printf("\ncontador_tipo: %d",contador[tipo]);
    alea_part = alea_i( 1, n_grupos[tipo] );

    n2 = part_a_mallado(alea_part,tipo);
    matriz_plasma[n2i]=matriz_plasma[n2];
    part_plasma[0][tipo] = part_plasma[alea_part][tipo];
    do{
        c_pared++;
        if(rechazo==1)part_plasma[alea_part][tipo]=part_plasma[0][tipo];
        do{
            idx = alea_i(-1,1);
            idy = alea_i(-1,1);
            if(c_pared>100){
                printf("\nidx: %f idy: %f",idx,idy);
            }
        }while(idx==0&&idy==0);
        //printf("\nANTES tipo: %d a_part: %d px_0: %f py_0: %f mx: %f my: %f",tipo,alea_part,part_plasma[0][tipo].x,part_plasma[0][tipo].y,matriz_plasma[n2].x,matriz_plasma[n2].y);
        part_plasma[alea_part][tipo].x += idx/reso;
        part_plasma[alea_part][tipo].y += idy/reso;
        //part_plasma[alea_part][tipo].x = alea_f(-R,R);
        //part_plasma[alea_part][tipo].y = alea_f(-R,R);
        n1 = part_a_mallado(alea_part,tipo);
        matriz_plasma[n1i]=matriz_plasma[n1];
        rechazo = 0;
        if(norma(part_plasma[alea_part][tipo].x, part_plasma[alea_part][tipo].y)>R)rechazo=1;
        //pared();
        if(c_pared>100){
            imprimir_celda_plasma_part(0,tipo);
            imprimir_celda_plasma_part(alea_part,tipo);
            printf("\ncontador_pared: %d",c_pared);
            getchar();
        }
        if(c_pared>110)rechazo=0;
    }while((rechazo==1));
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
    if(rho>=R){
        rechazo=1;
        //printf("\nQue show!");
        //getchar();
    }
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda_plasma(int a){
    printf("\nCelda: %i rho: %f phi: %f\nelectrones: %lld h20: %lld hp: %lld h2p: %lld carga: %lld",a,matriz_plasma[a].rho,matriz_plasma[a].phi,matriz_plasma[a].part[0],matriz_plasma[a].part[3],matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].carga);
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda_plasma_vec(int a, int b){
    printf("\nCelda: %i x: %f y: %f\n vx[%i]: %f vy[%i]: %f\nelectrones: %lld h20: %lld hp: %lld h2p: %lld carga: %lld",a,matriz_plasma[a].x,matriz_plasma[a].y,b,matriz_plasma[a].vx[b],matriz_plasma[a].vy[b],matriz_plasma[a].part[0],matriz_plasma[a].part[3],matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].carga);
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda_plasma_vec2(int a, int b){
    printf("\nCelda: %i x: %f y: %f\n vx[%i]: %f vy[%i]: %f\nelectrones: %lld h20: %lld hp: %lld h2p: %lld carga: %lld",a,matriz_plasma[a].x,matriz_plasma[a].y,b,matriz_plasma2[a].vx[b],b,matriz_plasma2[a].vy[b],matriz_plasma2[a].part[0],matriz_plasma2[a].part[3],matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].carga);
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
    de_coulomb = ef - ei;
    if(p>=terma){
        /*if(fabs(autoenergia(n1i))>fabs(ei)||fabs(autoenergia(n2i))>fabs(ei)||fabs(autoenergia(n1))>fabs(ef)||fabs(autoenergia(n2))>fabs(ef)){
            printf("\ntipo: %i ef: %e ei: %e\naE_n1i: %e aE_n2i: %e aE_n1: %e aE_n2: %e",tipo,ef,ei,autoenergia(n1i),autoenergia(n2i),autoenergia(n1),autoenergia(n2));
            bool1 = true;
        }*/
        /*printf("\npaso: %d tipo: %d n1: %d n2: %d n1.q: %lld n2.q: %lld E_n1n2: %e AE_n1: %e AE_n2: %e",p,tipo,n1,n2,matriz_plasma[n1].carga,matriz_plasma[n2].carga,(qe*qe*matriz_plasma[n1].carga*matriz_plasma[n2].carga)/(4*pi*epce*epsi*distancianormal(n1,n2)*esc),autoenergia(n1),autoenergia(n2));
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
    float aenergia, k_auto;
    //aenergia = k_auto*(6*matriz_plasma[a].carga*matriz_plasma[a].carga*qe*qe)/(5*4*pi*epce*epsi*matriz_plasma[a].anchow);
    k_auto = 1.4866047991;//2.0*log(1+sqrt(2))+2.0*(1-sqrt(2))/3.0;
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
        dt_termico+=(esc*distancianormal(n1,n2)*part_plasma[alea_part][tipo].part)/part_plasma[alea_part][tipo].v;
        //printf("\np: %d d: %f part_desp: %lld",p,distancianormal(n1,n2),part_plasma[alea_part][tipo].part);
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
void dinamica(void){
    int i,j;
    int celdaobj,nmov;
    //long long int suma1[4]={0},suma2[4]={0},suma3[4]={0},suma4[4]={0},suma5[4]={0},particulas[2822][4]={0};
    float vxx, vyy, dd, tt, giro_radio;
    /*FILE *daat;
    daat = fopen("velocidades_dina.txt","w");
    for(i=1;i<=nceldas;i++){
        fprintf(daat,"%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i,matriz_plasma[i].vx[0],matriz_plasma[i].vy[0],matriz_plasma[i].vx[1],matriz_plasma[i].vy[1],matriz_plasma[i].vx[2],matriz_plasma[i].vy[2],matriz_plasma[i].vx[3],matriz_plasma[i].vy[3]);
    }
    fclose(daat);
    getchar();/*
    //retermalizacion();
    /*for(i=1;i<=nceldas;i++){
        suma1[0] += matriz_plasma[i].part[0];
        suma1[1] += matriz_plasma[i].part[1];
        suma1[2] += matriz_plasma[i].part[2];
        suma1[3] += matriz_plasma[i].part[3];
    }
    for(i=1;i<=nceldas;i++){
        suma2[0] += matriz_plasma2[i].part[0];
        suma2[1] += matriz_plasma2[i].part[1];
        suma2[2] += matriz_plasma2[i].part[2];
        suma2[3] += matriz_plasma2[i].part[3];
    }*/
    //printf("\nelec1: %lld\thp1: %lld\th2p1: %lld\th201: %lld\n",suma1[0],suma1[1],suma1[2],suma1[3]);
    //printf("elec2: %lld\thp2: %lld\th2p2: %lld\th202: %lld\n",suma2[0],suma2[1],suma2[2],suma2[3]);
    for(i=1;i<=nceldas;i++){
        for(j=0;j<=3;j++){
            matriz_plasma2[i].part[j] = 0;
            matriz_plasma2[i].vx[j] = 0;
            matriz_plasma2[i].vy[j] = 0;
            //matriz_plasma2[i].vx_1[j] = 0;
            //matriz_plasma2[i].vy_2[j] = 0;
        }
    }
    //float normaa=0;
    for(i=1;i<=nceldas;i++){//mover particulas debido a la dinamica(8
        /*suma3[0]=suma3[1]=suma3[2]=suma3[3]=0;
        suma4[0]=suma4[1]=suma4[2]=suma4[3]=0;
        suma5[0]=suma5[1]=suma5[2]=suma5[3]=0;
        if(matriz_plasma[celda_dir_rec_tipob(i,0)].rho>R){
            printf("ASDASD,0\n");
            celda_dir_rec_tipob_error(i,0);
        }
        if(matriz_plasma[celda_dir_rec_tipob(i,1)].rho>R){
            printf("ASDASD,1\n");
            celda_dir_rec_tipob_error(i,1);
        }
        if(matriz_plasma[celda_dir_rec_tipob(i,2)].rho>R){
            printf("ASDASD,2\n");
            celda_dir_rec_tipob_error(i,2);
        }
        if(matriz_plasma[celda_dir_rec_tipob(i,3)].rho>R){
            printf("ASDASD,3\n");
            celda_dir_rec_tipob_error(i,3);
        }*/
        /*for(j=0;j<4;j++){
            if(matriz_plasma[i].part[j]!=0){
                if((masa[j]*matriz_plasma[i].v[j])/(qe*B)>0){
                    giro_radio = (carga[j]!=0)?(masa[j]*matriz_plasma[i].v[j])/(qe*B):1e10;
                    celdaobj = celda_dir_rec_tipob(i,j);
                    //dd = 2*giro_radio*asin(esc*distancianormal(i,celdaobj)/(2*giro_radio));
                    //dd = giro_radio*atan(esc*distancianormal(i,celdaobj)/giro_radio);
                    dd = esc*distancianormal(i,celdaobj)/reso;
                    tt = dd/matriz_plasma[i].v[j];
                    ndt[i] += tt;

                    vxx = matriz_plasma[i].vx[j];
                    vyy = matriz_plasma[i].vy[j];
                    matriz_plasma[i].vx[j] = vxx*cos(carga[j]*tt*qe*B/masa[j]) + vyy*sin(carga[j]*tt*qe*B/masa[j]);
                    matriz_plasma[i].vy[j] = -vxx*sin(carga[j]*tt*qe*B/masa[j]) + vyy*cos(carga[j]*tt*qe*B/masa[j]);
                    matriz_plasma2[celdaobj].vx[j] += matriz_plasma[i].part[j]*matriz_plasma[i].vx[j];
                    matriz_plasma2[celdaobj].vy[j] += matriz_plasma[i].part[j]*matriz_plasma[i].vy[j];
                    matriz_plasma2[celdaobj].part[j] += matriz_plasma[i].part[j];
                    matriz_plasma[i].part[j] = 0;
                    matriz_plasma[i].vx[j] = 0;
                    matriz_plasma[i].vy[j] = 0;
                    matriz_plasma[i].v[j] = 0;
                }
                else{
                    vxx = matriz_plasma[i].vx[j];
                    vyy = matriz_plasma[i].vy[j];
                    matriz_plasma2[i].vx[j] += matriz_plasma[i].part[j]*vxx;
                    matriz_plasma2[i].vy[j] += matriz_plasma[i].part[j]*vyy;
                    matriz_plasma2[i].part[j] += matriz_plasma[i].part[j];
                    matriz_plasma[i].part[j] = 0;
                    matriz_plasma[i].vx[j] = 0;
                    matriz_plasma[i].vy[j] = 0;
                    matriz_plasma[i].v[j] = 0;
                }
            }
        }*/
        if(matriz_plasma[i].part[0]!=0||npart[0]!=0){//&&(masa[0]*matriz_plasma[i].v[0])/(carga[0]*qe*B)>1.0*esc/reso){
            if((masa[0]*matriz_plasma[i].v[0])/(qe*B)>1.0*esc/reso){
            //if((masa[0]*matriz_plasma[i].v[0])/(qe*B)>0){
                giro_radio = (masa[0]*matriz_plasma[i].v[0])/(qe*B);
                celdaobj = celda_dir_rec_tipob(i,0);
                //particulas[i][0]+=matriz_plasma[i].part[0];
                /*if(matriz_plasma[celdaobj].rho > R){
                    printf("\nTIPO: 0\trho: %f\n",matriz_plasma[celdaobj].rho);
                    getchar();
                }*/
                /*printf("celda_dir_rec: %i\tcelda_dir_rec_tipob: %i\n",celda_dir_rec(i),celda_dir_rec_tipob(i,0));
                getchar();*/
                //celdaobj = celda_dir_rec_neg(i);
                //dd = 2*giro_radio*asin(esc*distancianormal(i,celdaobj)/(2*giro_radio));
                //dd = giro_radio*atan(esc*distancianormal(i,celdaobj)/giro_radio);
                dd = esc*distancianormal(i,celdaobj)/reso;
                tt = dd/matriz_plasma[i].v[0];
                ndt[i] += tt;
                //celdaobj = celda_dir(i);
                /*CORRECCION AHORITAprintf("\n0 radio: %e\tdd: %e\ttt: %e\tv[0]: %e",giro_radio,dd,tt,matriz_plasma[i].v[0]);
                printf("\ndistancia: %e\tesc: %e\targ_asin: %e",distancianormal(i,celdaobj),esc,esc*distancianormal(i,celdaobj)/(2*giro_radio));
                getchar();*/

                vxx = matriz_plasma[i].vx[0];
                vyy = matriz_plasma[i].vy[0];
                //matriz_plasma[i].vx[0] = vxx + ndt[i]*acel[0];
                //matriz_plasma[i].vy[0] = vyy + ndt[i]*acel[1];
                matriz_plasma[i].vx[0] = vxx*cos(carga[0]*tt*qe*B/masa[0]) + vyy*sin(carga[0]*tt*qe*B/masa[0]);
                matriz_plasma[i].vy[0] = -vxx*sin(carga[0]*tt*qe*B/masa[0]) + vyy*cos(carga[0]*tt*qe*B/masa[0]);
                /*printf("\nELECTRON\nCelda inicial: %d\tCelda objetivo: %d\tdd: %e\ttt: %e",i,celdaobj,dd,tt);
                printf("\nANTES\nvxx_0: %e\tvyy_0: %e\tv_0: %e",vxx,vyy,sqrt(vxx*vxx+vyy*vyy));
                printf("\nDESPUES\nvx_0: %e\tvvy_0: %e\tv_0: %e\n",matriz_plasma[i].vx[0],matriz_plasma[i].vy[0],sqrt(matriz_plasma[i].vx[0]*matriz_plasma[i].vx[0]+matriz_plasma[i].vy[0]*matriz_plasma[i].vy[0]));*/
                //getchar();
                matriz_plasma2[celdaobj].vx[0] += matriz_plasma[i].part[0]*matriz_plasma[i].vx[0];
                matriz_plasma2[celdaobj].vy[0] += matriz_plasma[i].part[0]*matriz_plasma[i].vy[0];
                matriz_plasma2[celdaobj].part[0] += matriz_plasma[i].part[0];
                /*if(matriz_plasma[i].vy[0]!=1){
                    printf("\nVALORES RAROS! vx: %f vy: %f\n",matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
                    getchar();
                }*/
                //getchar();
                matriz_plasma[i].part[0] = 0;
                matriz_plasma[i].vx[0] = 0;
                matriz_plasma[i].vy[0] = 0;
                matriz_plasma[i].v[0] = 0;
                /*if(i==nmatriz_plasma[ 62 ][ 62 ]){
                    printf("\nceldaobj: %i\n",celdaobj);
                }*/
                /*celdaobj = int(alea()*nceldas)+1;
                if(celdaobj==nceldas+1)celdaobj=nceldas;
                nmov = int( alea_f(1, matriz_plasma[i].part[0] ) );
                matriz_plasma[i].part[0]-= nmov;
                matriz_plasma[celdaobj].part[0]+= nmov;
                matriz_plasma[i].t[0]=distancianormal(i,celdaobj)/matriz_plasma[i].v[0];

                celdaobj = int(alea()*nceldas)+1;
                if(celdaobj==nceldas+1)celdaobj=nceldas;
                nmov = int( alea_f(1, matriz_plasma[i].part[1] ) );
                matriz_plasma[i].part[1]-= nmov;
                matriz_plasma[celdaobj].part[1]+= nmov;
                matriz_plasma[i].t[1]=distancianormal(i,celdaobj)/matriz_plasma[i].v[1];

                celdaobj = int(alea()*nceldas)+1;
                if(celdaobj==nceldas+1)celdaobj=nceldas;
                nmov = int( alea_f(1, matriz_plasma[i].part[2] ) );
                matriz_plasma[i].part[2]-= nmov;
                matriz_plasma[celdaobj].part[2]+= nmov;
                matriz_plasma[i].t[2]=distancianormal(i,celdaobj)/matriz_plasma[i].v[2];*/
            }
            else{
                //printf("\n0 giroradio: %f",(masa[0]*matriz_plasma[i].v[0])/(qe*B));
                //printf("\n0 masa: %e\tv: %f\tqe: %e\tB: %f",masa[0],matriz_plasma[i].v[0],qe,B);
                //particulas[i][0]+=matriz_plasma[i].part[0];
                tt = 0;
                vxx = matriz_plasma[i].vx[0];
                vyy = matriz_plasma[i].vy[0];
                //matriz_plasma[i].vx[0] = vxx + ndt[i]*acel[0];
                //matriz_plasma[i].vy[0] = vyy + ndt[i]*acel[1];
                matriz_plasma[i].vx[0] = vxx*cos(carga[0]*tt*qe*B/masa[0]) + vyy*sin(carga[0]*tt*qe*B/masa[0]);
                matriz_plasma[i].vy[0] = -vxx*sin(carga[0]*tt*qe*B/masa[0]) + vyy*cos(carga[0]*tt*qe*B/masa[0]) ;
                matriz_plasma2[i].vx[0] += matriz_plasma[i].part[0]*matriz_plasma[i].vx[0];
                matriz_plasma2[i].vy[0] += matriz_plasma[i].part[0]*matriz_plasma[i].vy[0];
                matriz_plasma2[i].part[0] += matriz_plasma[i].part[0];
                matriz_plasma[i].part[0] = 0;
                matriz_plasma[i].vx[0] = 0;
                matriz_plasma[i].vy[0] = 0;
                matriz_plasma[i].v[0] = 0;
            }
        }
        if(matriz_plasma[i].part[1]!=0||npart[1]!=0){//&&(masa[1]*matriz_plasma[i].v[1])/(carga[1]*qe*B)>1.0*esc/reso){
            if((masa[1]*matriz_plasma[i].v[1])/(qe*B)>1.0*esc/reso){
            //if((masa[1]*matriz_plasma[i].v[1])/(qe*B)>0){
                giro_radio = (masa[1]*matriz_plasma[i].v[1])/(qe*B);
                celdaobj = celda_dir_rec_tipob(i,1);
                //particulas[i][1]+=matriz_plasma[i].part[1];
                if(matriz_plasma[celdaobj].rho > R){
                    printf("\nTIPO: 1\trho: %f\n",matriz_plasma[celdaobj].rho);
                    getchar();
                }
                //dd = 2*giro_radio*asin(esc*distancianormal(i,celdaobj)/(2*giro_radio));
                //dd = giro_radio*atan(esc*distancianormal(i,celdaobj)/giro_radio);
                dd = esc*distancianormal(i,celdaobj)/reso;
                tt = dd/matriz_plasma[i].v[1];
                ndt[i] += tt;

                vxx = matriz_plasma[i].vx[1];
                vyy = matriz_plasma[i].vy[1];
                matriz_plasma[i].vx[1] = vxx*cos(carga[1]*tt*qe*B/masa[1]) + vyy*sin(carga[1]*tt*qe*B/masa[1]);
                matriz_plasma[i].vy[1] = vyy*cos(carga[1]*tt*qe*B/masa[1]) - vxx*sin(carga[1]*tt*qe*B/masa[1]);
                matriz_plasma2[celdaobj].vx[1] += matriz_plasma[i].part[1]*matriz_plasma[i].vx[1];
                matriz_plasma2[celdaobj].vy[1] += matriz_plasma[i].part[1]*matriz_plasma[i].vy[1];
                matriz_plasma2[celdaobj].part[1] += matriz_plasma[i].part[1];
                /*if(matriz_plasma[i].vy[0]!=1){
                    printf("\nVALORES RAROS! vx: %f vy: %f\n",matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
                    getchar();
                }*/
                //getchar();
                /*printf("\nHP\nCelda inicial: %d\tCelda objetivo: %d\tdd: %e\ttt: %e",i,celdaobj,dd,tt);
                printf("\nANTES\nvxx_1: %e\tvyy_1: %e\tv_1: %e\n",vxx,vyy,sqrt(vxx*vxx+vyy*vyy));
                printf("\nDESPUES\nvx_1: %e\tvvy_1: %e\tv_1: %e\n",matriz_plasma[i].vx[1],matriz_plasma[i].vy[1],sqrt(matriz_plasma[i].vx[1]*matriz_plasma[i].vx[1]+matriz_plasma[i].vy[1]*matriz_plasma[i].vy[1]));*/
                //getchar();
                /*printf("Entro aqui v;\ni: %d celdabj: %d",i,celdaobj);
                printf("\nx: %f\ty: %f",matriz_plasma[i].x,matriz_plasma[i].y);
                printf("\nvx_1: %f\tvy_1: %f",matriz_plasma[i].vx[1],matriz_plasma[i].vy[1]);
                printf("\nCELDA OBJETIVO\nx: %f\ty: %f",matriz_plasma[celdaobj].x,matriz_plasma[celdaobj].y);
                printf("\nvx_1: %f\tvy_1: %f",matriz_plasma[celdaobj].vx[1],matriz_plasma[celdaobj].vy[1]);
                getchar();*/
                matriz_plasma[i].part[1] = 0;
                matriz_plasma[i].vx[1] = 0;
                matriz_plasma[i].vy[1] = 0;
                matriz_plasma[i].v[1] = 0;
            }
            else{
                //printf("\n1 giroradio: %f",(masa[1]*matriz_plasma[i].v[1])/(qe*B));
                //printf("\n1 masa: %f\tv: %f\tqe: %e\tB: %f",masa[1],matriz_plasma[i].v[1],qe,B);
                //particulas[i][1]+=matriz_plasma[i].part[1];
                tt = 0;
                vxx = matriz_plasma[i].vx[1];
                vyy = matriz_plasma[i].vy[1];
                matriz_plasma[i].vx[1] = vxx*cos(carga[1]*tt*qe*B/masa[1]) + vyy*sin(carga[1]*tt*qe*B/masa[1]);
                matriz_plasma[i].vy[1] = vyy*cos(carga[1]*tt*qe*B/masa[1]) - vxx*sin(carga[1]*tt*qe*B/masa[1]);
                matriz_plasma2[i].vx[1] += matriz_plasma[i].part[1]*matriz_plasma[i].vx[1];
                matriz_plasma2[i].vy[1] += matriz_plasma[i].part[1]*matriz_plasma[i].vy[1];
                matriz_plasma2[i].part[1] += matriz_plasma[i].part[1];
                matriz_plasma[i].part[1] = 0;
                matriz_plasma[i].vx[1] = 0;
                matriz_plasma[i].vy[1] = 0;
                matriz_plasma[i].v[1] = 0;
            }
        }
        if(matriz_plasma[i].part[2]!=0||npart[2]!=0){//&&(masa[2]*matriz_plasma[i].v[2])/(carga[2]*qe*B)>1.0*esc/reso){
            if((masa[2]*matriz_plasma[i].v[2])/(qe*B)>1.0*esc/reso){
            //if((masa[2]*matriz_plasma[i].v[2])/(qe*B)>0){
                giro_radio = (masa[2]*matriz_plasma[i].v[2])/(qe*B);
                celdaobj = celda_dir_rec_tipob(i,2);
                //particulas[i][2]+=matriz_plasma[i].part[2];
                if(matriz_plasma[celdaobj].rho > R){
                    printf("\nTIPO: 2\trho: %f\n",matriz_plasma[celdaobj].rho);
                    getchar();
                }
                //dd = 2*giro_radio*asin(esc*distancianormal(i,celdaobj)/(2*giro_radio));
                //dd = giro_radio*atan(esc*distancianormal(i,celdaobj)/giro_radio);
                dd = esc*distancianormal(i,celdaobj)/reso;
                tt = dd/matriz_plasma[i].v[2];
                ndt[i] += tt;

                vxx = matriz_plasma[i].vx[2];
                vyy = matriz_plasma[i].vy[2];
                matriz_plasma[i].vx[2] = vxx*cos(carga[2]*tt*qe*B/masa[2]) + vyy*sin(carga[2]*tt*qe*B/masa[2]);
                matriz_plasma[i].vy[2] = vyy*cos(carga[2]*tt*qe*B/masa[2]) - vxx*sin(carga[2]*tt*qe*B/masa[2]);
                matriz_plasma2[celdaobj].vx[2] += matriz_plasma[i].part[2]*matriz_plasma[i].vx[2];
                matriz_plasma2[celdaobj].vy[2] += matriz_plasma[i].part[2]*matriz_plasma[i].vy[2];
                matriz_plasma2[celdaobj].part[2] += matriz_plasma[i].part[2];
                /*if(matriz_plasma[i].vy[0]!=1){
                    printf("\nVALORES RAROS! vx: %f vy: %f\n",matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
                    getchar();
                }*/
                //getchar();
                /*printf("\nH2P\nCelda inicial: %d\tCelda objetivo: %d\tdd: %e\ttt: %e",i,celdaobj,dd,tt);
                printf("\nANTES\nvxx_2: %e\tvyy_2: %e\tv_2: %e\n",vxx,vyy,sqrt(vxx*vxx+vyy*vyy));
                printf("\nDESPUES\nvx_2: %e\tvvy_2: %e\tv_2: %e\n",matriz_plasma[i].vx[2],matriz_plasma[i].vy[2],sqrt(matriz_plasma[i].vx[2]*matriz_plasma[i].vx[2]+matriz_plasma[i].vy[2]*matriz_plasma[i].vy[2]));*/
                //getchar();
                matriz_plasma[i].part[2] = 0;
                matriz_plasma[i].vx[2] = 0;
                matriz_plasma[i].vy[2] = 0;
                matriz_plasma[i].v[2] = 0;
            }
            else{
                //printf("\n2 giroradio: %f",(masa[2]*matriz_plasma[i].v[2])/(qe*B));
                //printf("\n2 masa: %f\tv: %f\tqe: %e\tB: %f",masa[2],matriz_plasma[i].v[2],qe,B);
                //particulas[i][2]+=matriz_plasma[i].part[2];
                tt = 0;
                vxx = matriz_plasma[i].vx[2];
                vyy = matriz_plasma[i].vy[2];
                matriz_plasma[i].vx[2] = vxx*cos(carga[2]*tt*qe*B/masa[2]) + vyy*sin(carga[2]*tt*qe*B/masa[2]);
                matriz_plasma[i].vy[2] = vyy*cos(carga[2]*tt*qe*B/masa[2]) - vxx*sin(carga[2]*tt*qe*B/masa[2]);
                matriz_plasma2[i].vx[2] += matriz_plasma[i].part[2]*matriz_plasma[i].vx[2];
                matriz_plasma2[i].vy[2] += matriz_plasma[i].part[2]*matriz_plasma[i].vy[2];
                matriz_plasma2[i].part[2] += matriz_plasma[i].part[2];
                matriz_plasma[i].part[2] = 0;
                matriz_plasma[i].vx[2] = 0;
                matriz_plasma[i].vy[2] = 0;
                matriz_plasma[i].v[2] = 0;
            }
        }
        if(matriz_plasma[i].part[3]!=0||npart[3]!=0){//&&(masa[3]*matriz_plasma[i].v[3])/(carga[3]*qe*B)>1.0*esc/reso){
            //if((masa[3]*matriz_plasma[i].v[3])/(carga[3]*qe*B)>1.0*esc/reso){
                celdaobj = celda_dir_rec_tipob(i,3);
                /*if(celdaobj==i){
                    printf("\ni: %d\tceldaobj: %d\nn %lld\tvx[%d]: %f\tvy[%d]: %f",i,celdaobj,matriz_plasma[i].part[3],3,matriz_plasma[i].vx[3],3,matriz_plasma[i].vy[3]);
                    getchar();
                }*/
                //particulas[i][3]+=matriz_plasma[i].part[3];
                if(matriz_plasma[celdaobj].rho > R){
                    printf("\nTIPO: 3\trho: %f\n",matriz_plasma[celdaobj].rho);
                    getchar();
                }
                dd = esc*distancianormal(i,celdaobj)/reso;
                tt = dd/matriz_plasma[i].v[3];
                ndt[i] += tt;

                vxx = matriz_plasma[i].vx[3];
                vyy = matriz_plasma[i].vy[3];
               // matriz_plasma[i].vx[3] = vxx*cos(carga[3]*tt*qe*B/masa[3]) + vyy*sin(carga[3]*tt*qe*B/masa[3]);
                //matriz_plasma[i].vy[3] = vyy*cos(carga[3]*tt*qe*B/masa[3]) - vxx*sin(carga[3]*tt*qe*B/masa[3]);
                matriz_plasma2[celdaobj].vx[3] += matriz_plasma[i].part[3]*matriz_plasma[i].vx[3];
                matriz_plasma2[celdaobj].vy[3] += matriz_plasma[i].part[3]*matriz_plasma[i].vy[3];
                matriz_plasma2[celdaobj].part[3] += matriz_plasma[i].part[3];
                /*if(matriz_plasma[i].vy[0]!=1){
                    printf("\nVALORES RAROS! vx: %f vy: %f\n",matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
                    getchar();
                }*/
                //getchar();
                matriz_plasma[i].part[3] = 0;
                matriz_plasma[i].vx[3] = 0;
                matriz_plasma[i].vy[3] = 0;
                matriz_plasma[i].v[3] = 0;
            //}
            /*else{
                particulas[i][3]+=matriz_plasma[i].part[3];
                vxx = matriz_plasma[i].vx[3];
                vyy = matriz_plasma[i].vy[3];
                matriz_plasma[i].vx[3] = vxx*cos(carga[3]*tt*qe*B/masa[3]) + vyy*sin(carga[3]*tt*qe*B/masa[3]);
                matriz_plasma[i].vy[3] = vyy*cos(carga[3]*tt*qe*B/masa[3]) - vxx*sin(carga[3]*tt*qe*B/masa[3]);
                matriz_plasma2[i].vx[3] += matriz_plasma[i].part[3]*matriz_plasma[i].vx[3];
                matriz_plasma2[i].vy[3] += matriz_plasma[i].part[3]*matriz_plasma[i].vy[3];
                matriz_plasma2[i].part[3] += matriz_plasma[i].part[3];
                matriz_plasma[i].part[3] = 0;
                matriz_plasma[i].vx[3] = 0;
                matriz_plasma[i].vy[3] = 0;
                matriz_plasma[i].v[3] = 0;

            }*/
        }
        /*for(j=1;j<=nceldas;j++){
            suma3[0]+=matriz_plasma2[j].part[0];
            suma3[1]+=matriz_plasma2[j].part[1];
            suma3[2]+=matriz_plasma2[j].part[2];
            suma3[3]+=matriz_plasma2[j].part[3];
        }
        for(j=1;j<=nceldas;j++){
            suma4[0]+=particulas[j][0];
            suma4[1]+=particulas[j][1];
            suma4[2]+=particulas[j][2];
            suma4[3]+=particulas[j][3];
        }
        for(j=1;j<=nceldas;j++){
            suma5[0]+=matriz_plasma[j].part[0];
            suma5[1]+=matriz_plasma[j].part[1];
            suma5[2]+=matriz_plasma[j].part[2];
            suma5[3]+=matriz_plasma[j].part[3];
        }
        if(suma3[0]!=suma4[0]||suma3[1]!=suma4[1]||suma3[2]!=suma4[2]||suma3[3]!=suma4[3]){
            printf("\nQUE SHOW i: %i\n",i);
            for(j=0;j<=3;j++){
                printf("suma3%i: %lld\t",j,suma3[j]);
            }
            printf("\n");
            for(j=0;j<=3;j++){
                printf("suma4%i: %lld\t",j,suma4[j]);
            }
            printf("\n");
            for(j=0;j<=3;j++){
                printf("suma5%i: %lld\t",j,suma5[j]);
            }
            printf("\nQue show");
            imprimir_celda_plasma_vec(i,0);
            imprimir_celda_plasma_vec2(celdaobj,0);
            getchar();
        }*/
    }
    /*for(i=1;i<=nceldas;i++){
        if(matriz_plasma[i].part[0]!=0||matriz_plasma[i].part[1]!=0||matriz_plasma[i].part[2]!=0||matriz_plasma[i].part[3]!=0){
            printf("\nQUE SHOW");
            imprimir_celda_plasma_rec(i);
            getchar();
        }
    }*/
    for(i=1;i<=nceldas;i++){
        for(j=0;j<4;j++){
            if(matriz_plasma2[i].part[j]!=0){
                vxx = (matriz_plasma[i].vx[j]*matriz_plasma[i].part[j]+matriz_plasma2[i].vx[j])/(matriz_plasma[i].part[j] + matriz_plasma2[i].part[j]);
                vyy = (matriz_plasma[i].vy[j]*matriz_plasma[i].part[j]+matriz_plasma2[i].vy[j])/(matriz_plasma[i].part[j] + matriz_plasma2[i].part[j]);
                if(matriz_plasma[i].part[j]&&matriz_plasma2[i].part[j]){
                    printf("\nQue show!\nvxx: %e\tvyy: %e\tm[%d].part[%d]: %lld\tm2[%d].part[%d]: %lld\n",vxx, vyy,i,j,matriz_plasma[i].part[j],i,j,matriz_plasma2[i].part[j]);
                    getchar();
                }
                matriz_plasma[i].vx[j] = vxx;
                matriz_plasma[i].vy[j] = vyy;
                matriz_plasma[i].v[j] = norma(matriz_plasma[i].vx[j],matriz_plasma[i].vy[j]);
                matriz_plasma[i].part[j] = matriz_plasma2[i].part[j];
            }
        }
        /*if(matriz_plasma2[i].part[0]!=0){
            vxx = (matriz_plasma[i].vx[0]*matriz_plasma[i].part[0]+matriz_plasma2[i].vx[0])/(matriz_plasma[i].part[0] + matriz_plasma2[i].part[0]);
            vyy = (matriz_plasma[i].vy[0]*matriz_plasma[i].part[0]+matriz_plasma2[i].vy[0])/(matriz_plasma[i].part[0] + matriz_plasma2[i].part[0]);
            matriz_plasma[i].vx[0] = vxx;
            matriz_plasma[i].vy[0] = vyy;
            matriz_plasma[i].v[0] = norma(matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
            matriz_plasma[i].part[0] = matriz_plasma2[i].part[0];
        }
        if(matriz_plasma2[i].part[1]!=0){
            vxx = (matriz_plasma[i].vx[1]*matriz_plasma[i].part[1]+matriz_plasma2[i].vx[1])/(matriz_plasma[i].part[1] + matriz_plasma2[i].part[1]);
            vyy = (matriz_plasma[i].vy[1]*matriz_plasma[i].part[1]+matriz_plasma2[i].vy[1])/(matriz_plasma[i].part[1] + matriz_plasma2[i].part[1]);
            matriz_plasma[i].vx[1] = vxx;
            matriz_plasma[i].vy[1] = vyy;
            matriz_plasma[i].v[1] = norma(matriz_plasma[i].vx[1],matriz_plasma[i].vy[1]);
            matriz_plasma[i].part[1] = matriz_plasma2[i].part[1];
        }
        if(matriz_plasma2[i].part[2]!=0){
            vxx = (matriz_plasma[i].vx[2]*matriz_plasma[i].part[2]+matriz_plasma2[i].vx[2])/(matriz_plasma[i].part[2] + matriz_plasma2[i].part[2]);
            vyy = (matriz_plasma[i].vy[2]*matriz_plasma[i].part[2]+matriz_plasma2[i].vy[2])/(matriz_plasma[i].part[2] + matriz_plasma2[i].part[2]);
            matriz_plasma[i].vx[2] = vxx;
            matriz_plasma[i].vy[2] = vyy;
            matriz_plasma[i].v[2] = norma(matriz_plasma[i].vx[2],matriz_plasma[i].vy[2]);
            matriz_plasma[i].part[2] = matriz_plasma2[i].part[2];
        }
        if(matriz_plasma2[i].part[3]!=0){
            vxx = (matriz_plasma[i].vx[3]*matriz_plasma[i].part[3]+matriz_plasma2[i].vx[3])/(matriz_plasma[i].part[3] + matriz_plasma2[i].part[3]);
            vyy = (matriz_plasma[i].vy[3]*matriz_plasma[i].part[3]+matriz_plasma2[i].vy[3])/(matriz_plasma[i].part[3] + matriz_plasma2[i].part[3]);
            matriz_plasma[i].vx[3] = vxx;
            matriz_plasma[i].vy[3] = vyy;
            matriz_plasma[i].v[3] = norma(matriz_plasma[i].vx[3],matriz_plasma[i].vy[3]);
            matriz_plasma[i].part[3] = matriz_plasma2[i].part[3];
        }*/
        calc_carga(i);
    }
    /*suma3[0]=suma3[1]=suma3[2]=suma3[3]=0;
    suma4[0]=suma4[1]=suma4[2]=suma4[3]=0;
    suma5[0]=suma5[1]=suma5[2]=suma5[3]=0;*/
    /*for(i=1;i<=nceldas;i++){
        suma3[0] += matriz_plasma[i].part[0];
        suma3[1] += matriz_plasma[i].part[1];
        suma3[2] += matriz_plasma[i].part[2];
        suma3[3] += matriz_plasma[i].part[3];
    }
    //printf("elec3: %lld\thp3: %lld\th2p3: %lld\th203: %lld\n",suma3[0],suma3[1],suma3[2],suma3[3]);
    for(i=1;i<=nceldas;i++){
        suma4[0] += matriz_plasma2[i].part[0];
        suma4[1] += matriz_plasma2[i].part[1];
        suma4[2] += matriz_plasma2[i].part[2];
        suma4[3] += matriz_plasma2[i].part[3];
    }
    for(i=1;i<=nceldas;i++){
        suma5[0] += particulas[i][0];
        suma5[1] += particulas[i][1];
        suma5[2] += particulas[i][2];
        suma5[3] += particulas[i][3];
    }*/
    //printf("elec4: %lld\thp4: %lld\th2p4: %lld\th204: %lld\n",suma4[0],suma4[1],suma4[2],suma4[3]);
    //printf("elec5: %lld\thp5: %lld\th2p5: %lld\th205: %lld\n",suma5[0],suma5[1],suma5[2],suma5[3]);
    //printf("delec: %lld\tdhp: %lld\tdh2p: %lld\tdh20: %lld\n",suma1[0]-suma4[0],suma1[1]-suma4[1],suma1[2]-suma4[2],suma1[3]-suma4[3]);
    //getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void retermalizacion(void){
    /*int i,j,xi,yi;
    for(i=0;i<R_reso;i++){
        for(j=0;j<4;j++){
            xi = nmatriz_plasma[ int(R*reso) + 1 ][ i + 1 ];
            yi = nmatriz_plasma[ i + 1 ][ int(R*reso) + 1 ];
            //matriz_plasma[ yi ].vy[j] = 2;
            //matriz_plasma[ xi ].vx[j] = -2;
            //printf("\nANTES vx: %f\tvy: %f",matriz_plasma[ xi ].vx[j],matriz_plasma[ yi ].vy[j]);
            matriz_plasma[ yi ].vy[j] = signo(matriz_plasma[ yi ].vy[j])*fabs( sqrt(kb*tempe[j]/masa[j]));//aaadist_normal(0.0,sqrt(kb*tempe[j]/masa[j])) );
            matriz_plasma[ xi ].vx[j] = signo(matriz_plasma[ xi ].vx[j])*fabs( sqrt(kb*tempe[j]/masa[j]));//aaadist_normal(0.0,sqrt(kb*tempe[j]/masa[j])) );
            //printf("\nDESPUES vx: %f\tvy: %f",matriz_plasma[ xi ].vx[j],matriz_plasma[ yi ].vy[j]);
            //getchar();
        }
        //printf("\ni: %d\tj: %d\tx_1: %f\ty_1: %f\tx_2: %f\ty_2: %f",i,j,matriz_plasma[ nmatriz_plasma[ i + 1 ][ int(R*reso)+1 ] ].x,matriz_plasma[ nmatriz_plasma[ i + 1 ][ int(R*reso)+1 ] ].y,matriz_plasma[ nmatriz_plasma[ int(R*reso)+1 ][ i + 1 ] ].x,matriz_plasma[ nmatriz_plasma[ int(R*reso)+1 ][ i + 1 ] ].y);
        //getchar();
    }*/
    double sumatemp[4] = {0}, retempe[4];
    for(int i=1; i<=nceldas; i++){
        sumatemp[0] += 0.5*masa[0]*matriz_plasma[i].part[0]*matriz_plasma[i].v[0]*matriz_plasma[i].v[0];
        sumatemp[1] += 0.5*masa[1]*matriz_plasma[i].part[1]*matriz_plasma[i].v[1]*matriz_plasma[i].v[1];
        sumatemp[3] += 0.5*masa[3]*matriz_plasma[i].part[3]*matriz_plasma[i].v[3]*matriz_plasma[i].v[3];
    }
    sumatemp[0] = sumatemp[0]/npart[0];
    sumatemp[1] = sumatemp[1]/npart[1];
    sumatemp[3] = sumatemp[3]/npart[3];
    retempe[0] = sumatemp[0]/qe;
    retempe[1] = sumatemp[1]/qe;
    retempe[3] = sumatemp[3]/qe;
    //printf("\nANTES\nretemp 0: %e 1: %e 3: %e",retempe[0],retempe[1],retempe[3]);
    //printf("\ntemp 0: %e 1: %e 3: %e",tempe[0],tempe[1],tempe[3]);
    for(int i=1; i<=nceldas; i++){
            matriz_plasma[i].vx[0] *= sqrt(tempe[0]/retempe[0]);
            matriz_plasma[i].vx[1] *= sqrt(tempe[1]/retempe[1]);
            matriz_plasma[i].vx[3] *= sqrt(tempe[3]/retempe[3]);
            matriz_plasma[i].vy[0] *= sqrt(tempe[0]/retempe[0]);
            matriz_plasma[i].vy[1] *= sqrt(tempe[1]/retempe[1]);
            matriz_plasma[i].vy[3] *= sqrt(tempe[3]/retempe[3]);
            matriz_plasma[i].v[0] = norma(matriz_plasma[i].vx[0],matriz_plasma[i].vy[0]);
            matriz_plasma[i].v[1] = norma(matriz_plasma[i].vx[1],matriz_plasma[i].vy[1]);
            matriz_plasma[i].v[3] = norma(matriz_plasma[i].vx[3],matriz_plasma[i].vy[3]);
    }
    /*sumatemp[0]=sumatemp[1]=sumatemp[3]=0;
    for(int i=1; i<=nceldas; i++){
        sumatemp[0] += 0.5*masa[0]*matriz_plasma[i].part[0]*matriz_plasma[i].v[0]*matriz_plasma[i].v[0];
        sumatemp[1] += 0.5*masa[1]*matriz_plasma[i].part[1]*matriz_plasma[i].v[1]*matriz_plasma[i].v[1];
        sumatemp[3] += 0.5*masa[3]*matriz_plasma[i].part[3]*matriz_plasma[i].v[3]*matriz_plasma[i].v[3];
    }
    sumatemp[0] = sumatemp[0]/npart[0];
    sumatemp[1] = sumatemp[1]/npart[1];
    sumatemp[3] = sumatemp[3]/npart[3];
    retempe[0] = sumatemp[0]/qe;
    retempe[1] = sumatemp[1]/qe;
    retempe[3] = sumatemp[3]/qe;
    printf("\nDESPUES\nretemp 0: %e 1: %e 3: %e",retempe[0],retempe[1],retempe[3]);
    printf("\ntemp 0: %e 1: %e 3: %e",tempe[0],tempe[1],tempe[3]);*/
    //getchar();

}
////////////////////////////////////////////////////////////////////////////////////
void actu_salida(void){
	int i, j;
	long long int cargatotaaal=0;
	float xx, yy;
    sprintf(salidac,"datos/posiciones%i.dat",(p-terma)/actu);
    if( dat=fopen(salidac,"w") ){
        fprintf(dat,"%i\t%i\t%f\t%f\t%f\t%f\t%d\t%d",p,terma,1.0*contador_a[0]/contador[0],((contador[1]>0)?1.0*contador_a[1]/contador[1]:0.0),1.0*contador_a[2]/contador[2],1.0*contador_a[3]/contador[3],contador_rechazo,rechazo_met_mov);
        for(i=1;i<=nceldas;i++){
            printposiciones(i);
        }
        fclose(dat);
    }

    sprintf(salidac,"datos/electron%i.dat",(p-terma)/actu);
    if( iprint[0]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(i=1;i<=nceldas;i++){
            //for(j=1;j<=8;j++){
                //xx = (int)(matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
                //yy = (int)(matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
                //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
                //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
                //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[0]);
                //fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[0]);
                //fprintf(dat,"%f\t%f\t%lld\n", xx, yy, matriz_plasma[i].part[0]);
                fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[0]);
            //}
        }
        fclose(dat);
    }

    sprintf(salidac,"datos/hp%i.dat",(p-terma)/actu);
    if( iprint[1]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"hp\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[1]);
        }
        fclose(dat);
    }

    sprintf(salidac,"datos/h2+%i.dat",(p-terma)/actu);
    if( iprint[2]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h2+\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[2]);
        }
        fclose(dat);
    }

    sprintf(salidac,"datos/h20%i.dat",(p-terma)/actu);
    if( iprint[3]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h20\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[3]);
        }
        fclose(dat);
    }

    sprintf(salidac,"datos/todas%i.dat",(p-terma)/actu);
    if( iprint[4]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[2]+matriz_plasma[i].part[1]+matriz_plasma[i].part[3]+matriz_plasma[i].part[0]);
        }
        fclose(dat);
    }

    sprintf(salidac,"datos/carga%i.dat",(p-terma)/actu);
    if( iprint[5]==1 ){
        dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"carga\"\n");
        fprintf(dat, "#X\tY\tZ\n");
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].carga);
            cargatotaaal += matriz_plasma[i].carga;
        }
        fprintf(dat,"carga total: %lld\n",cargatotaaal);
        fclose(dat);
    }
}
////////////////////////////////////////////////////////////////////////////////////
void actu_salida2(void){
	long long int cargatotaaal=0;
	int max_n_grupos = n_grupos[0];
	float xx, yy;
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
        fprintf(dat, "#X\tY\tZ\n");
        for(int i=1;i<=nceldas;i++){
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].carga);
            cargatotaaal += matriz_plasma[i].carga;
        }
        fprintf(dat,"carga total: %lld\n",cargatotaaal);
        fclose(dat);
    }
}
////////////////////////////////////////////////////////////////////////////////////
void salida_prom(void){
	int i, j;
	//long long int cargatotaaal=0;
	//float xx, yy;

    sprintf(salidac,"datos/promedios/gelectron.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"electrones\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = (int)(matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
            //yy = (int)(matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[0]);
            //fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].part[0]);
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, gelectron[i]/(p-terma));
        //}
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/h2+.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h2+\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, gh2p[i]/(p-terma));
        //}
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/h20.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h20\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, gh20[i]/(p-terma));
        //}
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/hp.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"hp\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, ghp[i]/(p-terma));
        //}
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/todas.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, (gelectron[i]+gh2p[i]+gh20[i]+ghp[i])/(p-terma));
        //}
    }
    fclose(dat);

    sprintf(salidac,"datos/promedios/carga.dat");
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"carga\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        //for(j=1;j<=8;j++){
            //xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            //yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%lld\n", matriz_plasma[i].x, matriz_plasma[i].y, (gh2p[i]+ghp[i]-gelectron[i])/(p-terma));
            //cargatotaaal += matriz_plasma[i].carga;
        //}
    }
    //fprintf(dat,"carga total: %lld\n",cargatotaaal);
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void salida_prom2(void){
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
                    d = distancianormal(i,j);
                    energia_i += (qe*qe*matriz_plasma[i].carga*matriz_plasma[j].carga)/(4*pi*epce*epsi*d*esc);
                }
            }
        }
        aenergia_total += autoenergia(i);
        ecoulomb_total += energia_i;
    }
    energia_total = aenergia_total + 0.5*ecoulomb_total;
    if(p>pasoinicial){
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
        //printf("\ncelda: %d ec_t: %e ae_t: %e e_total: %e d_%d_%d: %f",i,ecoulomb_total,aenergia_total,energia_total,dummy_i,dummy_j,distancianormal(dummy_i,dummy_j));
        imprimir_celda_plasma_rec(dummy_i);
        imprimir_celda_plasma_rec(dummy_j);
        printf("\nae_%d: %e ae_%d: %e ec: %e",dummy_i,autoenergia(dummy_i),dummy_j,autoenergia(dummy_j),(qe*qe*matriz_plasma[dummy_i].carga*matriz_plasma[dummy_j].carga)/(4*pi*epce*epsi*distancianormal(dummy_i,dummy_j)*esc));
        getchar();
    }*/
    return(energia_total);
}
////////////////////////////////////////////////////////////////////////////////////
float energia_cinetica(int a){
    int i,j;
    float suma=0;
    char nombre[50];
    for(i=1;i<=nceldas;i++){
        for(j=0;j<=3;j++){
            if(matriz_plasma[i].part[j]!=0){
                suma += 0.5*matriz_plasma[i].part[j]*masa[j]*matriz_plasma[i].v[j]*matriz_plasma[i].v[j];
                /*if(i==20&&j==0){
                    printf("\nEnergia[%d][%d]: %e",i,j,0.5*matriz_plasma[i].part[j]*masa[j]*matriz_plasma[i].v[j]*matriz_plasma[i].v[j]);
                    getchar();
                }*/
            }
        }
    }
    sprintf(nombre,"histogramas/energia_cinetica%d.dat",a);
    if( dat = fopen(nombre,"w") ){
        for(i=1;i<=nceldas;i++){
            fprintf(dat,"%d",i);
            for(j=0;j<=3;j++){
                fprintf(dat,"\t%e",0.5*matriz_plasma[i].part[j]*masa[j]*matriz_plasma[i].v[j]*matriz_plasma[i].v[j]);
            }
            fprintf(dat,"\n",i);
        }
    }
    fclose(dat);
    if( dat = fopen("energia_cinetica.dat","a") ){
        fprintf(dat,"%d\t%e\n",a,suma);
        fclose(dat);
    }
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
void printposiciones(int a){
//fprintf(dat,"\n%4.0i\t%9.0I64d\t%9.0I64d\t%9.0I64d\t%9.0I64d\t%9.0i",i,matriz_plasma[i].part[0],matriz_plasma[i].part[3],matriz_plasma[i].part[1],matriz_plasma[i].part[2],matriz_plasma[i].carga);
    if(matriz_plasma[a].part[0]==0){
        if(matriz_plasma[a].part[3]==0 ){
            if(matriz_plasma[a].part[1]==0){
                if(matriz_plasma[a].part[2]==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%s\t%s\t%s\t%s",a,cero,cero,cero,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%s\t%s\t%s\t%10.0lld",a,cero,cero,cero,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%s\t%s\t%10.0lld\t%s",a,cero,cero,cero,matriz_plasma[a].part[2],cero);

                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%s\t%s\t%10.0lld\t%10.0lld",a,cero,cero,cero,matriz_plasma[a].part[2],matriz_plasma[a].carga);
                    }
                }
            }
            else{
                if(matriz_plasma[a].part[2]==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%s\t%10.0lld\t%s\t%s",a,cero,cero,matriz_plasma[a].part[1],cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%s\t%10.0lld\t%s\t%10.0lld",a,cero,cero,matriz_plasma[a].part[1],cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%s\t%10.0lld\t%10.0lld\t%s",a,cero,cero,matriz_plasma[a].part[1],matriz_plasma[a].part[2],cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%s\t%10.0lld\t%10.0lld\t%10.0lld",a,cero,cero,matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].carga);
                    }
                }
            }
        }
        else{
            if(matriz_plasma[a].part[1]==0){
                if(matriz_plasma[a].part[2]==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%10.0lld\t%s\t%s\t%s",a,cero,matriz_plasma[a].part[3],cero,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%10.0lld\t%s\t%s\t%10.0lld",a,cero,matriz_plasma[a].part[3],cero,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%10.0lld\t%s\t%10.0lld\t%s",a,cero,matriz_plasma[a].part[3],cero,matriz_plasma[a].part[2],cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%10.0lld\t%s\t%10.0lld\t%10.0lld",a,cero,matriz_plasma[a].part[3],cero,matriz_plasma[a].part[2],matriz_plasma[a].carga);
                    }
                }
            }
            else{
                if(matriz_plasma[a].part[2]==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%10.0lld\t%10.0lld\t%s\t%s",a,cero,matriz_plasma[a].part[3],matriz_plasma[a].part[1],cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%10.0lld\t%10.0lld\t%s\t%10.0lld",a,cero,matriz_plasma[a].part[3],matriz_plasma[a].part[1],cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%s\t%10.0lld\t%10.0lld\t%10.0lld\t%s",a,cero,matriz_plasma[a].part[3],matriz_plasma[a].part[1],matriz_plasma[a].part[2],cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%s\t%10.0lld\t%10.0lld\t%10.0lld\t%10.0lld",a,cero,matriz_plasma[a].part[3],matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].carga);
                    }
                }
            }
        }
    }
    else{
        if(matriz_plasma[a].part[3]==0 ){
            if(matriz_plasma[a].part[1]==0){
                if(matriz_plasma[a].part[2]==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0lld\t%s\t%s\t%s\t%s",a,matriz_plasma[a].part[0],cero,cero,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0lld\t%s\t%s\t%s\t%10.0lld",a,matriz_plasma[a].part[0],cero,cero,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0lld\t%s\t%s\t%10.0lld\t%s",a,matriz_plasma[a].part[0],cero,cero,matriz_plasma[a].part[2],cero);

                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0lld\t%s\t%s\t%10.0lld\t%10.0lld",a,matriz_plasma[a].part[0],cero,cero,matriz_plasma[a].part[2],matriz_plasma[a].carga);
                    }
                }
            }
            else{
                if(matriz_plasma[a].part[2]==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0lld\t%s\t%10.0lld\t%s\t%s",a,matriz_plasma[a].part[0],cero,matriz_plasma[a].part[1],cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0lld\t%s\t%10.0lld\t%s\t%10.0lld",a,matriz_plasma[a].part[0],cero,matriz_plasma[a].part[1],cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0lld\t%s\t%10.0lld\t%10.0lld\t%s",a,matriz_plasma[a].part[0],cero,matriz_plasma[a].part[1],matriz_plasma[a].part[2],cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0lld\t%s\t%10.0lld\t%10.0lld\t%10.0lld",a,matriz_plasma[a].part[0],cero,matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].carga);
                    }
                }
            }
        }
        else{
            if(matriz_plasma[a].part[1]==0){
                if(matriz_plasma[a].part[2]==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0lld\t%10.0lld\t%s\t%s\t%s",a,matriz_plasma[a].part[0],matriz_plasma[a].part[3],cero,cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0lld\t%10.0lld\t%s\t%s\t%10.0lld",a,matriz_plasma[a].part[0],matriz_plasma[a].part[3],cero,cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0lld\t%10.0lld\t%s\t%10.0lld\t%s",a,matriz_plasma[a].part[0],matriz_plasma[a].part[3],cero,matriz_plasma[a].part[2],cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0lld\t%10.0lld\t%s\t%10.0lld\t%10.0lld",a,matriz_plasma[a].part[0],matriz_plasma[a].part[3],cero,matriz_plasma[a].part[2],matriz_plasma[a].carga);
                    }
                }
            }
            else{
                if(matriz_plasma[a].part[2]==0){
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0lld\t%10.0lld\t%10.0lld\t%s\t%s",a,matriz_plasma[a].part[0],matriz_plasma[a].part[3],matriz_plasma[a].part[1],cero,cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0lld\t%10.0lld\t%10.0lld\t%s\t%10.0lld",a,matriz_plasma[a].part[0],matriz_plasma[a].part[3],matriz_plasma[a].part[1],cero,matriz_plasma[a].carga);
                    }
                }
                else{
                    if(matriz_plasma[a].carga==0){
                        fprintf(dat,"\n%i\t%10.0lld\t%10.0lld\t%10.0lld\t%10.0lld\t%s",a,matriz_plasma[a].part[0],matriz_plasma[a].part[3],matriz_plasma[a].part[1],matriz_plasma[a].part[2],cero);
                    }
                    else{
                        fprintf(dat,"\n%i\t%10.0lld\t%10.0lld\t%10.0lld\t%10.0lld\t%10.0lld",a,matriz_plasma[a].part[0],matriz_plasma[a].part[3],matriz_plasma[a].part[1],matriz_plasma[a].part[2],matriz_plasma[a].carga);
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
        gelectron[i] += matriz_plasma[i].part[0];
        gh20[i] += matriz_plasma[i].part[3];
        ghp[i] += matriz_plasma[i].part[1];
        gh2p[i] += matriz_plasma[i].part[2];
    }
}
////////////////////////////////////////////////////////////////////////////////////
void gdr_plasma2(void){
    for(int i=1;i<=nceldas;i++){
        for(int j=0;j<ncomp;j++)gpart[i][j] += matriz_plasma[i].part[j];
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir(int a){
    float normav, normar, rpv, phipv;
    int int_rho, int_phi, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vx[0]*matriz_plasma[a].vx[0] + matriz_plasma[a].vy[0]*matriz_plasma[a].vy[0] );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vx[0];  v1[1] = matriz_plasma[a].vy[0]; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/

    rpv = (matriz_plasma[a].vx[0]*matriz_plasma[a].x)/(normar*normav) + (matriz_plasma[a].vy[0]*matriz_plasma[a].y)/(normar*normav);
    phipv = -(matriz_plasma[a].vx[0]*matriz_plasma[a].y)/(normar*normav) + (matriz_plasma[a].vy[0]*matriz_plasma[a].x)/(normar*normav);

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
    float normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vx[0]*matriz_plasma[a].vx[0] + matriz_plasma[a].vy[0]*matriz_plasma[a].vy[0] );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    /*if(matriz_plasma[a].part[0]!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vx[0];  v1[1] = matriz_plasma[a].vy[0]; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    //pphi = (atan2(matriz_plasma[a].vy[0],matriz_plasma[a].vx[0])+((atan2(matriz_plasma[a].vy[0],matriz_plasma[a].vx[0]))<0?2*pi:0.0));
    //printf("pphi: %f",pphi);
    //getchar();
    if(matriz_plasma[a].vx[0]/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vx[0]/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vy[0]/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vy[0]/normav<1.0/3.0){
        int_y = 0;
    }else{
        int_y = 1;
    }
    int_i = floor( matriz_plasma[a].x*reso) + floor(R*reso) + 1 + int_x;
    int_j = floor( matriz_plasma[a].y*reso) + floor(R*reso) + 1 + int_y;
    /*if(a==10000000000000000){
        printf("\nANTES DE CELDA_ERROR\nCELDA INICIAL: %i i: %i j: %i", a, (int)( matriz_plasma[a].x*reso + R*reso + 1 ), (int)( matriz_plasma[a].y*reso + R*reso + 1 ));
        printf("\nCELDA OBJETIVO: %i i: %i j: %i", nmatriz_plasma[ int_i ][ int_j ], int_i, int_j );
        printf("\nint_x: %i int_y: %i",int_x, int_y);
        printf("\ny*reso: %f int(y*reso): %i\n",matriz_plasma[a].y*reso, (((int) (matriz_plasma[a].y*reso + 32768.5)) - 32768));
        printf("\nR*reso: %f int(R*reso): %i\n",R*reso, int(R*reso));
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vx[0],matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vy[0]);
        imprimir_celda_plasma_rec(a);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        getchar();
    }*/
    /*if(nmatriz_plasma[ int_i ][ int_j ]>=12345||nmatriz_plasma[ int_i ][ int_j ]==0){
        printf("int_i: %i\tint_j: %i",int_i,int_j);
        return( nmatriz_plasma[ int_i ][ int_j ] );*/
    //rrho = sqrt(matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x+matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y);
    if(matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho>R){
        if(a==221099000){
            printf("Que show r: %f",sqrt(int_i*int_i+int_j*int_j));
            getchar();
        }
        return( nmatriz_plasma[ int_i-int_x ][ int_j-int_y ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec_tipob(int a, int b){
    float normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = norma(matriz_plasma[a].vx[b],matriz_plasma[a].vy[b]);//sqrt( matriz_plasma[a].vx[b]*matriz_plasma[a].vx[b] + matriz_plasma[a].vy[b]*matriz_plasma[a].vy[b] );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );

    if(matriz_plasma[a].vx[b]/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vx[b]/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vy[b]/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vy[b]/normav<1.0/3.0){
        int_y = 0;
    }else{
        int_y = 1;
    }
    //printf("\nQue show 1! a: %d x: %f y: %f flr: %f flx: %f fly: %f ix: %d iy: %d",a,matriz_plasma[a].x,matriz_plasma[a].y,floor( R*reso ),floor( matriz_plasma[a].x*reso ),floor( matriz_plasma[a].y*reso ),int_x,int_y);
    int_i = floor( matriz_plasma[a].x*reso ) + floor(R*reso) + 1 + int_x;
    int_j = floor( matriz_plasma[a].y*reso ) + floor(R*reso) + 1 + int_y;
    if(int_i>R_reso)int_i=R_reso;
    if(int_j>R_reso)int_j=R_reso;
    //printf("\tint_x: %d\tint_y: %d\tint_i: %d\tint_j: %d",int_x,int_y,int_i,int_j);
    //printf("ASD");
    //printf("\nint_x: %d\tint_y: %d\tint_i: %d\tint_j: %d\tnmatriz: %d\trho: %f",int_x,int_y,int_i,int_j,nmatriz_plasma[int_i][int_j],matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho);
    /*if(a==10000000000000000){
        printf("\nANTES DE CELDA_ERROR\nCELDA INICIAL: %i i: %i j: %i", a, (int)( matriz_plasma[a].x*reso + R*reso + 1 ), (int)( matriz_plasma[a].y*reso + R*reso + 1 ));
        printf("\nCELDA OBJETIVO: %i i: %i j: %i", nmatriz_plasma[ int_i ][ int_j ], int_i, int_j );
        printf("\nint_x: %i int_y: %i",int_x, int_y);
        printf("\ny*reso: %f int(y*reso): %i\n",matriz_plasma[a].y*reso, (((int) (matriz_plasma[a].y*reso + 32768.5)) - 32768));
        printf("\nR*reso: %f int(R*reso): %i\n",R*reso, int(R*reso));
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vx[0],matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vy[0]);
        imprimir_celda_plasma_rec(a);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        getchar();
    }*/
    /*if(nmatriz_plasma[int_i][int_j]==0){
        printf("a: %i\tint_i: %i\tint_j: %i\tint_x: %i\tint_y: %i",a,int_i,int_j,int_x,int_y);
        getchar();
    }*/
    /*if(nmatriz_plasma[int_i][int_j]==0){
        printf("QUE SHOWWWWWW\ncelda: %i celdaobj: %i",a,nmatriz_plasma[int_i][int_j]);
        printf("\nint_x: %i\tint_y: %i\tint_i: %i\tint_j: %i",int_x,int_y,int_i,int_j);
        getchar();
    }*/
    if(matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho>R||int_i==0||int_j==0||int_i>R_reso||int_j>R_reso||(int_x==0&&int_y==0)){
        /*if(a==221099000){
            printf("Que show r: %f",sqrt(int_i*int_i+int_j*int_j));
            getchar();
        }*/
        return( nmatriz_plasma[ int_i-int_x ][ int_j-int_y ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
void celda_dir_rec_tipob_error(int a, int b){
    float normav, normar;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vx[b]*matriz_plasma[a].vx[b] + matriz_plasma[a].vy[b]*matriz_plasma[a].vy[b] );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );

    if(matriz_plasma[a].vx[b]/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vx[b]/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vy[b]/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vy[b]/normav<1.0/3.0){
        int_y = 0;
    }else{
        int_y = 1;
    }
    int_i = floor( matriz_plasma[a].x*reso) + floor(R*reso) + 1 + int_x;
    int_j = floor( matriz_plasma[a].y*reso) + floor(R*reso) + 1 + int_y;
    printf("a: %i\tx: %f\ty: %f\n",a,matriz_plasma[a].x,matriz_plasma[a].y);
    printf("vx: %f\tvy: %f\n",matriz_plasma[a].vx[b]/normav,matriz_plasma[a].vy[b]/normav);
    printf("int_i: %i\tint_j: %i\tint_x: %i\tint_y: %i\n",int_i,int_j,int_x,int_y);
    getchar();
    if(matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho>R||int_i==0||int_j==0){
        //return( nmatriz_plasma[ int_i-int_x ][ int_j-int_y ] );
    }else{
        //return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec_error(int a){
    float normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vx[0]*matriz_plasma[a].vx[0] + matriz_plasma[a].vy[0]*matriz_plasma[a].vy[0] );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    /*if(matriz_plasma[a].part[0]!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vx[0];  v1[1] = matriz_plasma[a].vy[0]; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    //pphi = (atan2(matriz_plasma[a].vy[0],matriz_plasma[a].vx[0])+((atan2(matriz_plasma[a].vy[0],matriz_plasma[a].vx[0]))<0?2*pi:0.0));
    //printf("pphi: %f",pphi);
    //getchar();
    if(matriz_plasma[a].vx[0]/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vx[0]/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vy[0]/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vy[0]/normav<1.0/3.0){
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
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vx[0],matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vy[0]);
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
        /*if(a==221099000){
            printf("Que show r: %f",sqrt(int_i*int_i+int_j*int_j));
            getchar();
        }*/
        return( nmatriz_plasma[ int_i-int_x ][ int_j-int_y ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec_hp(int a){
    float normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vx[1]*matriz_plasma[a].vx[1] + matriz_plasma[a].vy[1]*matriz_plasma[a].vy[1] );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    /*if(matriz_plasma[a].part[0]!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vx[0];  v1[1] = matriz_plasma[a].vy[0]; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    //pphi = (atan2(matriz_plasma[a].vy[0],matriz_plasma[a].vx[0])+((atan2(matriz_plasma[a].vy[0],matriz_plasma[a].vx[0]))<0?2*pi:0.0));
    //printf("pphi: %f",pphi);
    //getchar();
    if(matriz_plasma[a].vx[1]/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vx[1]/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vy[1]/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vy[1]/normav<1.0/3.0){
        int_y = 0;
    }else{
        int_y = 1;
    }
    int_i = floor( matriz_plasma[a].x*reso) + floor(R*reso) + 1 + int_x;
    int_j = floor( matriz_plasma[a].y*reso) + floor(R*reso) + 1 + int_y;
    if(matriz_plasma[a].part[0]>=1e30){
        printf("celda inicial: %i celda objetivo: %i int_i: %i int_j:%i\n",a,nmatriz_plasma[ int_i ][ int_j ],int_i,int_j);
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vx[0],matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vy[0]);
        imprimir_celda_plasma_rec(a);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        getchar();
    }
    //rrho = sqrt(matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].x+matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y*matriz_plasma[ nmatriz_plasma[ int_i ][ int_j ] ].y);
    if(matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho>R){
            /*printf("Que hongo, que show\ni: %i r: %f\n",nmatriz_plasma[int_i][int_j],matriz_plasma[ nmatriz_plasma[int_i][int_j] ].rho);
            printf("i: %i r: %f",nmatriz_plasma[int_i-int_x][int_j-int_y],matriz_plasma[ nmatriz_plasma[int_i-int_x][int_j-int_y] ].rho);
            getchar();*/
        return( nmatriz_plasma[ int_i-int_x ][ int_j-int_y ] );
    }else{
        return( nmatriz_plasma[ int_i ][ int_j ] );
    }
}
////////////////////////////////////////////////////////////////////////////////////
int celda_dir_rec_error_hp(int a){
    float normav, normar;//, rrho, pphi;
    int int_x, int_y, int_i, int_j;
    normav = sqrt( matriz_plasma[a].vx[1]*matriz_plasma[a].vx[1] + matriz_plasma[a].vy[1]*matriz_plasma[a].vy[1] );
    normar = sqrt( matriz_plasma[a].x*matriz_plasma[a].x + matriz_plasma[a].y*matriz_plasma[a].y );
    /*if(matriz_plasma[a].part[0]!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vx[0];  v1[1] = matriz_plasma[a].vy[0]; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    //pphi = (atan2(matriz_plasma[a].vy[0],matriz_plasma[a].vx[0])+((atan2(matriz_plasma[a].vy[0],matriz_plasma[a].vx[0]))<0?2*pi:0.0));
    //printf("pphi: %f",pphi);
    //getchar();
    if(matriz_plasma[a].vx[1]/normav<-(1.0/3.0)){
        int_x = -1;
    }else if(matriz_plasma[a].vx[1]/normav<1.0/3.0){
        int_x = 0;
    }else{
        int_x = 1;
    }
    if(matriz_plasma[a].vy[1]/normav<-(1.0/3.0)){
        int_y = -1;
    }else if(matriz_plasma[a].vy[1]/normav<1.0/3.0){
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
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vx[1],matriz_plasma[a].vy[1]);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vx[1],matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vy[1]);
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
    float normav, normar, rrho, xx, yy, dt = 2e0;//campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    float acelneg[2];
    int int_i, int_j, nceldaobj;
    /*if(matriz_plasma[a].part[0]!=0){
        printf("a: %i vx: %f vy: %f",a,matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        getchar();
    }*/
    //if(fabs(normar - matriz_plasma[a].rho)>1e-5)printf("celda: %i normar: %f rho_a: %f resta: %e\n", a, normar, matriz_plasma[a].rho, normar - matriz_plasma[a].rho);
    /*float v1[3], v2[3];
    vector<float> d1, d2, d3;
    v1[0] = matriz_plasma[a].vx[0];  v1[1] = matriz_plasma[a].vy[0]; v1[2] = 0.0;
    v2[0] = 0.0;  v2[1] = 0.0; v2[2] = campo_magnetico[int(matriz_plasma[a].rho*10)].bz;
    d1 = vect(v1); d2 = vect(v2); d3 = cruz(d1,d2);*/
    acelneg[0] = -qe*matriz_plasma[a].vy[0]*B/masa[0];
    acelneg[1] = qe*matriz_plasma[a].vx[0]*B/masa[0];
    xx = matriz_plasma[a].x + matriz_plasma[a].vx[0]*dt + 0.5*acelneg[0]*dt*dt;
    yy = matriz_plasma[a].y + matriz_plasma[a].vy[0]*dt + 0.5*acelneg[1]*dt*dt;
    printf("x_a_0: %f y_a_0: %f\tx_a_t: %f y_a_t: %f",matriz_plasma[a].x,matriz_plasma[a].y,xx,yy);
    printf("\ndato1: %f dato2: %f dato3: %f dato4: %f", matriz_plasma[a].x,matriz_plasma[a].vx[0],dt,acelneg[0] );
    printf("\ndato1: %f dato2: %f dato3: %f dato4: %f", matriz_plasma[a].y,matriz_plasma[a].vy[0],dt,acelneg[1] );
    getchar();
    int_i = (int)( xx + 61 + 1 );
    int_j = (int)( yy + 61 + 1 );
    nceldaobj = nmatriz_plasma[ int_i ][ int_j ];
    /*if(matriz_plasma[a].part[0]>=1e30){
        printf("celda inicial: %i celda objetivo: %i int_i: %i int_j:%i\n",a,nmatriz_plasma[ int_i ][ int_j ],int_i,int_j);
        printf("vx_a: %f vy_a: %f\n",matriz_plasma[a].vx[0],matriz_plasma[a].vy[0]);
        printf("vx_obj: %f vy_onj: %f\nAAAA",matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vx[0],matriz_plasma[nmatriz_plasma[ int_i ][ int_j ]].vy[0]);
        imprimir_celda_plasma_rec(a);
        printf("\nOBJ");
        imprimir_celda_plasma_rec(nmatriz_plasma[ int_i ][ int_j ]);
        getchar();
    }*/
    rrho = sqrt(matriz_plasma[ nceldaobj ].x*matriz_plasma[ nceldaobj ].x+matriz_plasma[ nceldaobj ].y*matriz_plasma[ nceldaobj ].y);
    if(rrho>R){
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
void dinamica2(void){
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
    for(int j=0; j<ncomp; j++){
        //for(int i=0; i<ncomp; i++){
            //printf("\nnpart2_%d: %d",i,n_grupos[i]);
            //getchar();
        //}
        if(n_grupos[j]>0){
            for(int i=1; i<=n_grupos[j]; i++){
                //printf("\nANTES DE DINAMICA i: %d j: %d\n",i,j);
                //if((10*i)%n_grupos[j]==0)printf("\ri: %d\n",i);
                //dinamica2();
                double x,y,vx,vy,ax,ay;
                float mass = masa[j], dt;
                //dt = dtau/100.0;
                if(solodinamica){
                    dt = 0.01*( B>1e-5?((2*pi*min_mass)/(100.0*max_charge*qe*B)):(esc/(100.0*max_v*reso)) );
                }else{
                    dt = dt_termico/(100.0*(npart[0]+npart[1]));
                }
                int charge = carga[j];
                double Ex=0.0,Ey=0.0,Bz=B;
                double x0=part_plasma[i][j].x*esc, v0x=part_plasma[i][j].vx, y0 = part_plasma[i][j].y*esc, v0y = part_plasma[i][j].vy;
                double periodo = (2*pi*mass)/(abs(charge)*qe*Bz);
                int ktau = 100;
                if(j>=ncomp-1)dt=0;
                //printf("\nq: %e m: %e q/m: %e",charge*qe,mass*m_u,charge*qe/mass/m_u);
                //printf("\nipart: %d p_x: %f p_y: %f x0: %f y0: %f",ipart,part_plasma[1][0].x,part_plasma[1][0].y,x0,y0);
                //printf("\ncarga: %d masa: %e charge: %d mass: %e",carga[0],masa[0],charge,mass);

                //dat=fopen("dinamica2.dat","w");
                //fprintf(dat,"#t\tx\tvx\ty\tvy\n");
                //fprintf(dat,"%f\t%f\t%f\t%f\t%f\n",0.0,x0,v0x,y0,v0y);
                //for(int i=1;i<=int(periodo/dt)+1;i++){
                long long int k;
                for(k=1;k<=ktau;k++){
                    ax = ((charge*qe)/(mass))*(Ex+v0y*Bz);
                    ay = ((charge*qe)/(mass))*(Ey-v0x*Bz);
                    x = x0 + v0x*dt + 0.5*ax*dt*dt;
                    vx = v0x + ax*dt;
                    y = y0 + v0y*dt + 0.5*ay*dt*dt;
                    vy = v0y + ay*dt;
                    /*printf("ax: %f ay: %f\n%f\t%f\t%f\t%f\t%f\n",ax,ay,i*dt,x,vx,y,vy);
                    getchar();*/
                    /*if(i%itau==0){
                        fprintf(dat,"%e\t%f\t%f\t%f\t%f\n",i*dt,x,vx,y,vy);
                        //if(jcomp==0)printf("\rp: %d j: %d itau: %d i: %lld",p,jcomp, itau, i);
                    }*/
                    //printf("\nx0: %e y0: %e vx0: %e vy0: %e",x0,y0,v0x,v0y);
                    //printf("\nx: %e y: %e vx: %e vy: %e ax: %e ay: %e",x,y,vx,vy,ax,ay);
                    //getchar();
                    x0 = x; v0x = vx;
                    y0 = y; v0y = vy;
                }
                part_plasma[i][j].x = x/esc;
                part_plasma[i][j].y = y/esc;
                part_plasma[i][j].vx = vx;
                part_plasma[i][j].vy = vy;
                //dinamica2(i,j,B>1e-5?((2*pi*min_mass)/(100.0*max_charge*qe*B)):(esc/(100.0*max_v*reso)) );
                //dinamica2(i,j,1e-7 );

                float rho2 = norma(part_plasma[i][j].x,part_plasma[i][j].y);
                if(rho2>R){
                    intercambiar(i,n_grupos[j],j);
                    n_grupos[j]--;
                }
            }
            for(int i=1; i<=n_grupos[j]; i++){
                //printf("\nANTES DE DINAMICA i: %d j: %d\n",i,j);
                float rho2 = norma(part_plasma[i][j].x,part_plasma[i][j].y);
                if(rho2>R){
                    printf("\ni: %d n_grupos[%d]: %d",i,j,n_grupos[j]);
                    getchar();
                }
            }
        }else{
            printf("\nn_grupos[%d]: %d",j,n_grupos[j]);
            getchar();
        }
    }

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
}
