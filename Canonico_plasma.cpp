#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <tgmath.h>
#include <windows.h>

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
void leer_datos_iniciales(void);
void leer_campo_magnetico(void);
void imprimir_datos_iniciales(void);
void imprimir_celda_plasma(int a);
void condiciones_iniciales(void);
void arreglo_inicial(void);
void crear_matriz_plasma(void);
void calc_carga(int a);

float distanciacelda(int a, int b);
float dist_ima(int a, int b, int c);
int signo(float a);
void mover_particulas(void);
void metropolis_plasma(int a);

float de_plasma(void);
float calcular_de_mov(void);
void gdr_plasma(void);
void actu_salida(void);
void salida_prom(void);
float energia(void);
float ecoulomb(int a, int b, int c);
void beta(void);

void printposiciones(int a);

//////////////////////////////Constantes
int pasos, actu, terma;
float esc, densidad, volumen, nelectronesr, tempee, tempei;
long long int nh20, nhp, nh2p, nelectrones;
float R;
float m1, m2, m3, m4;
///////////////////////////////////////////////////////////////////////Variables globales
int p=0, tipo, ienergia;
int rechazo;
int nceldas, n1, n2, n1i, n2i;                                //numero de celda Nueva/Vieja Fnicial/Final
char cero[] = "         0";
///////////////////////////////////////////////////////////////////////Contadores
int c_mov=0,c_mova=1;
int c_uno=0,c_unoa=0,c_dos=0,c_dosa=0,c_tres=0,c_tresa=0;
int rechazo_neg=0, rechazo_met_mov=0;

long long int gelectron[1410], gh20[1410], gh2p[1410], ghp[1410];

struct smatriz_plasma{
	float rho,phi,x,y;
	int carga;
	long long int h20,hp,h2p,electrones;
    float k1x, k1y, k1z, k2x, k2y, k2z, volumen;
}matriz_plasma[61*50+1];

struct sarreglo{                        //pongo el +1 porque me gusta empezar los ciclos desde i,j = 1
    float z;
    float r[36+1];                      //36 valores posibles para r [0.0,3.5] de 0.1 en 0.1
    float br[36+1];
    float bz[36+1];
    float absb[36+1+24];
    float A[36+1];
    float dbz[36+1];
    float dbr[36+1];
    float fieldindex[36+1];
}arreglo[371+1];                        //371 valores posibles para z [0.0,37.0] de 0.1 en 0.1

int nmatriz_plasma[60+1][50];
////////////////////////////////////////////////////////////////////////Variables para salida
FILE *dat;
char salidac[10];

main(){
    int i,j,k;
    srand((unsigned)time(NULL));
    /*system("mkdir datos");
    system("mkdir datos/datitos");*/
    CreateDirectory ("datos", NULL);
    CreateDirectory ("datos/promedios", NULL);
    for(i=1;i<=1407;i++){
        matriz_plasma[i].electrones=matriz_plasma[i].h20=matriz_plasma[i].h2p=matriz_plasma[i].hp=matriz_plasma[i].carga = 0;
    }

	condiciones_iniciales();
	crear_matriz_plasma();
	arreglo_inicial();
	leer_campo_magnetico();
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
        if(p%1000==0){
            if(p%100000==0){
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
void leer_datos_iniciales(){
    float dummy;
    dat=fopen("entrada.txt","r");
    fscanf(dat,"Numero de pasos: %i\n", &pasos);
    fscanf(dat,"Actualizacion: %i\n", &actu);
    fscanf(dat,"Termalizacion: %i\n", &terma);
    fscanf(dat,"Masas (kg): %f, %f, %f, %f\n",&m1,&m2,&m3,&m4);
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
    printf("Masas (kg): %e, %e, %e, %e\n",m1,m2,m3,m4);
    printf("R (mm): %f\n",R);
    printf("Temperatura electrones (K): %f\n",tempee);
    printf("Temperatura iones (K): %f\n\n",tempei);
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_iniciales(){
    leer_datos_iniciales();
    esc = 1e-3;//diam*1e-9;

    densidad = 2e19;
    volumen = pi*R*esc*R*esc*0.125*1e-6;
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

    printf("te: %f ti: %f",tempee,tempei);
    getchar();
    getchar();
    getchar();
    getchar();
    getchar();
    getchar();

    //printf("\ntempei: %f tempee: %f\n", tempei, tempee);
    //getchar();

    imprimir_datos_iniciales();
}

////////////////////////////////////////////////////////////////////////////////////
void crear_matriz_plasma(){
    int i,j,contadorm,contadorvec=0;

    contadorm=1;
    nmatriz_plasma[1][1] = 1;
    //matriz_plasma[1].rho = 1/sqrt(2);
    matriz_plasma[1].rho = 1-0.5;
    matriz_plasma[1].phi = pi/8;
    matriz_plasma[1].x = matriz_plasma[1].rho*cos(matriz_plasma[1].phi);
    matriz_plasma[1].y = matriz_plasma[1].rho*sin(matriz_plasma[1].phi);
    matriz_plasma[1].volumen = pi/8;

    for(i=2;i<=R;i++){
        for(j=1;j<=(int)((pi*i)*0.25);j++){
            contadorm++;
            nmatriz_plasma[i][j]=contadorm;
            matriz_plasma[contadorm].rho = sqrt(pow(i,2)+pow(i-1,2))/sqrt(2);
            //matriz_plasma[contadorm].rho = i-0.5;
            matriz_plasma[contadorm].phi = pi/(8*((int)((pi*i)*0.25)))+(j-1)*(pi/(4*((int)((pi*i)*0.25))));
            matriz_plasma[contadorm].x = matriz_plasma[contadorm].rho*cos(matriz_plasma[contadorm].phi);
            matriz_plasma[contadorm].y = matriz_plasma[contadorm].rho*sin(matriz_plasma[contadorm].phi);
            matriz_plasma[contadorm].volumen =  (pi*(i*i-(i-1)*(i-1))*0.125*esc*esc*1e-6)/(int)((pi*i)*0.25);
        }
    }
    nceldas=contadorm;
    printf("nceldas: %i",nceldas);
    getchar();
    getchar();
    getchar();
    getchar();
    getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void calc_carga(int a){
    matriz_plasma[a].carga = matriz_plasma[a].h2p + matriz_plasma[a].hp - matriz_plasma[a].electrones;
}
////////////////////////////////////////////////////////////////////////////////////
void arreglo_inicial(){
    int i,j,ni,dummy,res=1000;
    long long int cargatotal=0,ee=0, h200=0, h2pp=0, hpp=0;

    dat = fopen("datos/posiciones.dat","r");
    for(i=1;i<=nceldas;i++){
        fscanf(dat,"\n%i\t%i\t%i\t%i\t%i\t%i",&dummy,&matriz_plasma[i].electrones,&matriz_plasma[i].h20,&matriz_plasma[i].hp,&matriz_plasma[i].h2p,&matriz_plasma[i].carga);
        if((matriz_plasma[i].electrones==0)||(matriz_plasma[i].h20==0)||(matriz_plasma[i].hp==0)||(matriz_plasma[i].h2p==0)||(matriz_plasma[i].carga==0)){
            imprimir_celda_plasma(i);
        }
    }
    fclose(dat);
    if(matriz_plasma[1].electrones==0){
        for(i=1;i<=(int)(nelectrones/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].electrones+=res;
        }
        for(i=1;i<=(int)(nh20/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].h20+=res;
        }
        for(i=1;i<=(int)(nhp/res);i++){
            ni = (int)((alea())*nceldas)+1;
            if(ni==nceldas+1)ni=nceldas;
            matriz_plasma[ni].hp+=res;
        }
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
    }
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

    actu_salida();
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
                d = distanciacelda(n1i,i);
                ei += (qe*qe*matriz_plasma[n1i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                d = distanciacelda(n1,i);
                ef += (qe*qe*matriz_plasma[n1].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                /*ei += 8*ecoulomb(n1i,i,1);
                ef += 8*ecoulomb(n1,i,1);*/
                d = distanciacelda(n2i,i);
                ei += (qe*qe*matriz_plasma[n2i].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                d = distanciacelda(n2,i);
                ef += (qe*qe*matriz_plasma[n2].carga*matriz_plasma[i].carga)/(4*pi*epce*epsi*d*esc);
                /*ei += 8*ecoulomb(n2i,i,1);
                ef += 8*ecoulomb(n2,i,1);*/
            }
            //getchar();
        }
    }
    d = distanciacelda(n1i,n2i);
    ei += (qe*qe*matriz_plasma[n1i].carga*matriz_plasma[n2i].carga)/(4*pi*epce*epsi*d*esc);
    //ei += 8*ecoulomb(n1i,n2i,1);
    d = distanciacelda(n1,n2);
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
        for(j=1;j<=8;j++){
            //xx = (int)(matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
            //yy = (int)(matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4)*10)/10.0;
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            //fprintf(dat,"%2.5f, %2.5f, %9.0i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
            //fprintf(dat,"%f\t%f\t%i\n", matriz_plasma[i].x, matriz_plasma[i].y, matriz_plasma[i].electrones);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].electrones);
        }
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
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].h2p);
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/h20%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"h20\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].h20);
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/hp%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"hp\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].hp);
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/todas%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"n_particulas\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%I64d\n", xx, yy, matriz_plasma[i].h2p+matriz_plasma[i].hp+matriz_plasma[i].h20+matriz_plasma[i].electrones);
        }
    }
    fclose(dat);

    sprintf(salidac,"datos/carga%i.dat",(p-terma)/actu);
    dat=fopen(salidac,"w");
    //fprintf(dat, "\"x\", \"y\", \"carga\"\n");
    fprintf(dat, "#X\tY\tZ\n");
    for(i=1;i<=nceldas;i++){
        for(j=1;j<=8;j++){
            xx = matriz_plasma[i].rho*cos(matriz_plasma[i].phi+(j-1)*pi/4);
            yy = matriz_plasma[i].rho*sin(matriz_plasma[i].phi+(j-1)*pi/4);
            fprintf(dat,"%f\t%f\t%i\n", xx, yy, matriz_plasma[i].carga);
            cargatotaaal += matriz_plasma[i].carga;
        }
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
                    d = distanciacelda(i,j);
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
