#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#define MAXPART 1000
#define pi 3.1415926535897932384
#define epce 8.8541878176e-12
#define epsi 78.5
#define qe 1.60219e-19
#define kb 1.3806488e-23
#define na 6.02217e23
const int maxl = 20;
const int maxw = 20;

float alea(void);
void posiciones_iniciales(void);
void leer_datos_iniciales(void);
void imprimir_datos_iniciales(void);
void imprimir_celda(int a);
void imprimir_particula(int a);
void imprimir_celda_f(int a);
void imprimir_particula_f(int a);
void condiciones_iniciales(void);
void crear_matriz(void);
void crear_matriz2(void);
void posicion_promedio_celda(int a);
void condiciones_periodicas(void);
void condiciones_periodicas2(void);
void esferas_duras(void);
void esferas_duras_part(int a);

float distancia(int a, int b);
float distanciaz(int a,int b);
float distanciacelda(int a, int b);
float distanciazcelda(int a, int b);
float distanciapartcel(int a, int b);
float distanciazpartcel(int a, int b);
void intercambiar(int a, int b);
int signo(float a);
int maximo(int a, int b);
int numcelda(int a);
int numcelda2(int a);
void elegir_particula(void);
void mover(void);
void insertar_particulas(void);
void sustraer_particulas(void);
void metropolis_mov(void);
void metropolis_insertar(void);
void metropolis_sustraer(void);

float de_insus(int a, int b);
float de_mov(void);
float calcular_de_mov(void);

float encuentrab(void);
void sumas(void);
void gdr(void);
void salida(void);
void salidageo(void);

//////////////////////////////Constantes
int pasos, actu, terma;
float winicial, w, l, diam, esc, csal, densidad, volumen, sigmainicial, sigma, ner, nnr, tempe;
int ne, nn, np, npart;
float r;
float dx, dy, dz;
int v1, v2;
float B, lngamma;

const int clases = 200;
float distrmp, distrmn, dw, dl;

///////////////////////////////////////////////////////////////////////Variables globales
int p=0;
int n1,n2;
int rechazo;
int nceldas, ncvi, ncvf, ncni, ncnf;                                //numero de celda Nueva/Vieja Fnicial/Final

long sumapx, sumapy, sumapz;
long sumanx, sumany, sumanz;
///////////////////////////////////////////////////////////////////////Variables para verificar
int mis;
float dxx,dyy,dzz;
int siestaenlosvecinos=0;
int particulainmovil=0;
int particulasobrepuesta=0;
int ni=0;
///////////////////////////////////////////////////////////////////////Contadores
int c_mov=1,c_mova=1,c_ins=1, c_insa=1, c_sus=1, c_susa=1;
int rechazo_met=0, rechazo_tras=0, rechazo_per=0, rechazo_met_mov=0, rechazo_esf_mov=0;
////////////////////////////////////////////////////////////////////////Arreglos

long dgdrpx[clases+1], dgdrnx[clases+1];
long dgdrpy[clases+1], dgdrny[clases+1];
long dgdrpz[clases+1], dgdrnz[clases+1];

struct particulas{
    float x,y,z;
    int carga,inmovilidad;
}part[MAXPART+1];
//int matriz[15][15][30][1][27][1][15][1];//x,y,z,#matriz,vecinos,#part,listapart,carga
struct smatriz{
	float x,y,z;
	int carga,nparticulas,nvecinos;
	int vecinos[50];
	int particulas[150];
}matriz[50*50*50+1];

int nmatriz[50][50][50];

struct smatriz2{
	float x,y,z;
	int carga,nparticulas,nvecinos;
	int vecinos[28];
	int particulas[1500];
}matriz2[61*50+1];

int nmatriz2[60+1][50];
////////////////////////////////////////////////////////////////////////Variables para salida
FILE *dat;
char salidac[10];

main(){
    int i,j,k;
    int cargatotal=0;
    for(i=0;i<=clases;i++){
        dgdrpx[i]=dgdrpy[i]=dgdrpz[i]=dgdrnx[i]=dgdrny[i]=dgdrnz[i]=0;
    }
    system("mkdir datos");
    system("mkdir temp");
    srand((unsigned)time(NULL));

	condiciones_iniciales();
	crear_matriz2();
	crear_matriz();
	posiciones_iniciales();
	salidageo();

	part[1000].x = 0.8;
    part[1000].y = 0.4;
    part[1000].x = part[1000].x*1000000;
    part[1000].y = part[1000].y*1000000;

    printf("\nparticula 1000\n x: %f y: %f ncelda2: %i",part[1000].x,part[1000].y,numcelda2(1000));
    getchar();

	/*part[0].x = 0.5;
	for(i=1;i<=winicial;i++){
		part[0].y = 0.5;
		for(j=1;j<=winicial;j++){
			part[0].z = 1.0;
			for(k=1;k<=(l-1);k++){
				printf("\nCelda %i:: x:%f y:%f z:%f\nParticula:: x:%f y:%f z:%f",numcelda(0),matriz[ numcelda(0) ].x, matriz[ numcelda(0) ].y, matriz[ numcelda(0) ].z, part[0].x, part[0].y, part[0].z);
        		getchar();
				if(part[0].z + 1<= (l-0.5) ){
				part[0].z++;
				}
			}
			if(part[0].y+1<=winicial)part[0].y++;
		}
		if(part[0].x+1<=winicial)part[0].x++;
	}*/

    printf("\nWr: %f Sigma: %e npart:%i, nn:%i, ne:%i, np: %i ",w,sigma,npart,nn,ne,np);
    printf("\nner: %f nnr: %f Volumen: %e  B: %f ",ner, nnr, volumen, B);
    printf("\ndistrmp: %f distrmn: %f", distrmp, distrmn);

    if(lngamma==0) B = encuentrab();
    for(i=0;i<=clases;i++){
        dgdrpx[i]=dgdrpy[i]=dgdrpz[i]=dgdrnx[i]=dgdrny[i]=dgdrnz[i]=0;
    }
    c_mov=c_mova=c_ins=c_insa=c_sus=c_susa=1;
    rechazo_met=rechazo_per=rechazo_tras=rechazo_met_mov=0;
    nn = nnr; np = nn+ne; npart = nn + np;

    for(i=1;i<=npart;i++){
        part[i].inmovilidad = 0;
    }

    /*for(i=1;i<=npart;i++){
        printf("\ni:%i,x:%f,y:%f,z:%f,q:%f",i,part[i].x,part[i].y,part[i].z,part[i].carga);
    }
    getchar();*/
    for(p=1;p<=pasos;p++){

        rechazo = 0;
        //if(p%1==0)esferas_duras();

        if(alea()<1.8&&npart>0){
            mis = 1;
            c_mov++;
            elegir_particula();
            mover();
            condiciones_periodicas();

            if((part[n1].inmovilidad>=15000)&&(p>=-1000000)){
                //printf("\nLa particula %i está inmovil",n1);
                //printf("\nx: %f y: %f z: %f carga: %i",part[0].x, part[0].y, part[0].z, part[0].carga);
                //printf("\ndxx: %f dyy: %f dzz: %f",dxx,dyy,dzz);
                //imprimir_particula(0);
                //imprimir_particula(n1);
                //esferas_duras_part(n1);

                particulainmovil=1;
            }
            else{
                particulainmovil=0;
            }

            if(rechazo == 0){
                metropolis_mov();
            }
            else{
                rechazo_per++;
                part[n1] = part[0];
                //printf("\n\nSe rechazo el movimiento por condiciones periodicas.");
            }
            //getchar();
        }
        else{
            if(alea()<0.5&&nn>=1){
                mis = 2;
                c_sus++;
                sustraer_particulas();
                metropolis_sustraer();
            }
            else{
                mis = 3;
                c_ins++;
                insertar_particulas();
                metropolis_insertar();
            }
        }
        cargatotal=0;
        for(i=1;i<=nceldas;i++){
            cargatotal += matriz[i].carga;
        }
        if(cargatotal!=ne){
            //printf("\nCARGA TOTAL: %i n1: %i n1celda: %i",cargatotal,n1,numcelda(n1));
            /*printf("\nPARTICULAS\n i    x      y      z     carga");
            for(i=1;i<=npart;i++){
                printf("\n%3.0i %2.3f %2.3f %2.3f %2.0i",i,part[i].x,part[i].y,part[i].z,part[i].carga);
            }*/
            salidageo();
            //getchar();
        }
        gdr();
        sumas();
        salida();
        if((p%actu==0)&&(p>terma))salidageo();
        /*if(p%1000000==0){
            sprintf(salidac,"posicion_celdas_%i.dat",p/(1000000));
            dat=fopen(salidac,"w");
            for(i=1;i<=nceldas;i++){
                imprimir_celda_f(i);
            }
            fclose(dat);
        }*/
        if(p%1000==0) printf("\rPaso: %i Acept. Mov: %1.5f Rechazo per: %1.5f Rechazo met: %1.5f Rechazo esf: %1.5f",p,c_mova/(c_mov*1.0),rechazo_per/(c_mov*1.0), rechazo_met_mov/(c_mov*1.0),rechazo_esf_mov/(c_mov*1.0));
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
    fscanf(dat,"W, L (en diametros): %f, %f\n",&winicial,&l);
    fscanf(dat,"Diametro (nm): %f\n",&diam);
    fscanf(dat,"Sigma (microC/m^2): %f\n",&sigmainicial);
    fscanf(dat,"Temperatura (K): %f\n",&tempe);
    fscanf(dat,"Concentracion salina (mol/l): %f\n",&csal);
    fscanf(dat,"dx, dy, dz: %f, %f, %f\n",&dx, &dy, &dz);
    fscanf(dat,"Log(gamma): %f\n",&lngamma);
    fclose(dat);

}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_datos_iniciales(){
    printf("Numero de pasos: %i\n", pasos);
    printf("Actualizacion: %i\n", actu);
    printf("Termalizacion: %i\n", terma);
    printf("Valencias: %i, %i\n",v1,v2);
    printf("W, L (en diametros): %f, %f\n",winicial,l);
    printf("Diametro (nm): %f\n",diam);
    printf("Sigma (microC/m^2): %f\n",sigmainicial);
    printf("Temperatura (K): %f\n",tempe);
    printf("Concentracion salina (mol/l): %f\n",csal);
    printf("dx, dy, dz: %f, %f, %f\n",dx, dy, dz);
    printf("Log(gamma): %f\n",lngamma);
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_iniciales(){
leer_datos_iniciales();
esc = diam*1e-9;
densidad = csal*1000*na;
sigma = sigmainicial/100;
ner = (sigma*winicial*winicial*esc*esc)/(qe);
ne = (int)(ner);
if(sigmainicial!=0){
    //w = sqrt((100*qe*ne)/(sigmainicial*esc*esc));
    w = winicial;
}else{
    w = winicial;
}
volumen = w*esc*w*esc*(l-1)*esc;
nnr = densidad*volumen;
nn = (int)(nnr);
np = nn + ne;
npart = nn + np;

nceldas = winicial*winicial*(l-1);
ncvi = 0;
ncvf = nceldas + 1;
ncni = nceldas + 2;
ncnf = nceldas + 3;

distrmp = (nnr+ne)/(1.0*clases);
distrmn = (nnr)/(1.0*clases);

distrmp = (nn+ne)/(1.0*clases);
distrmn = (nn)/(1.0*clases);

dw=w/(1.0*clases);
dl=(l-1)/(1.0*clases);

B = log(nnr*nnr)+2*lngamma;
imprimir_datos_iniciales();
}
////////////////////////////////////////////////////////////////////////////////////
void crear_matriz(){
    int i,j,k,i1,j1,k1,contadorm=0,contadorvec=0;
    float xx, yy, zz;

    //printf("\nHASTA AQUI SI LLEGA");
    //getchar();

    for(i=1;i<=winicial;i++){
        for(j=1;j<=winicial;j++){
            for(k=1;k<=(l-1);k++){
                contadorm++;
                nmatriz[i][j][k]=contadorm;
                matriz[contadorm].x = i;
                matriz[contadorm].y = j;
                matriz[contadorm].z = k;
            }
        }
    }
    nceldas = contadorm;
    for(i=1;i<=nceldas;i++){
        contadorvec=0;
        for(j=1;j<=nceldas;j++){
            xx = matriz[j].x; yy = matriz[j].y; zz = matriz[j].z;
            if(fabs(matriz[i].x - matriz[j].x)==(int)(winicial-1)){
            //if(fabs(matriz[i].x - matriz[j].x)>(w*0.5)){
                if(matriz[i].x - matriz[j].x>0){
                    xx += (int)(winicial);
                    //xx += w;
                }
                else{
                    xx -= (int)(winicial);
                    //xx -= w;
                }
            }
            if(fabs(matriz[i].y - matriz[j].y)==(int)(winicial-1)){
            //if(fabs(matriz[i].y - matriz[j].y)>(w*0.5)){
                if(matriz[i].y - matriz[j].y>0){
                    yy += (int)(winicial);
                    //yy += w;
                }
                else{
                    yy -= (int)(winicial);
                    //yy -= w;
                }
            }
            if( (fabs(matriz[i].x-xx)<2)&&(fabs(matriz[i].y-yy)<2)&&(fabs(matriz[i].z-zz)<2) ){
                contadorvec++;
                matriz[i].nvecinos++;
                matriz[i].vecinos[contadorvec] = j;
            }
        }
    }
    for(i=1;i<=nceldas;i++){
    	if(matriz[i].x==(int)(winicial)){
		    matriz[i].x -= (1+(winicial-w))/2;
		}else{
			matriz[i].x -= 0.5;
		}
		if(matriz[i].y==(int)(winicial)){
		    matriz[i].y -= (1+(winicial-w))/2;
		}else{
			matriz[i].y -= 0.5;
		}
    }
    /*printf("\n%f %f",w,l);
    for(i=1;i<=nceldas;i++){
        printf("\n%4.0i %4.3f %4.3f %4.3f",i,matriz[i].x,matriz[i].y,matriz[i].z);
    }
    getchar();*/
    sprintf(salidac,"vecinos.dat");
    dat=fopen(salidac,"w");
    for(i=1;i<=nceldas;i++){
        fprintf(dat,"\n%3.0i    ",i);
        for(j=1;j<=matriz[i].nvecinos;j++){
        	fprintf(dat,"%4.0i",matriz[i].vecinos[j]);
		}
    }
    fclose(dat);
}
////////////////////////////////////////////////////////////////////////////////////
void crear_matriz2(){
    int i,j,contadorm=0,contadorvec=0;
    float xx, yy;
    r = 60000000;                                                                   //Radio en cm, todo lo demas en nm
    contadorm++;
    nmatriz2[1][1] = 1;
    for(i=2;i<=(int)(r/(1000000));i++){
        for(j=1;j<=(int)((pi*i)*0.25);j++){
            contadorm++;
            nmatriz2[i][j] = contadorm;
            matriz2[contadorm].x = i;
            matriz2[contadorm].y = j;
        }
    }
    nceldas=contadorm;
    for(i=1;i<=nceldas;i++){
        contadorvec = 0;
        for(j=1;j<=nceldas;j++){
            xx = matriz2[j].x; yy = matriz2[j].y;
            if(matriz2[i].y == 1){
                if(yy==(int)((pi*matriz2[j].x)*0.25)){

                }
            }

        }
        contadorvec++;
        matriz2[i].nvecinos++;
        matriz2[i].vecinos[contadorvec] = i;
    }
    sprintf(salidac,"vecinos2.dat");
    dat=fopen(salidac,"w");
    for(i=1;i<=nceldas;i++){
        fprintf(dat,"\n%3.0i    ",i);
        for(j=1;j<=27;j++){
        	fprintf(dat,"%4.0i",matriz2[i].vecinos[j]);
		}
    }
    fclose(dat);
    /*struct smatriz2{
        float x,y,z;
        int carga,nparticulas,nvecinos;
        int vecinos[28];
        int particulas[1500];
    }matriz2[61*50+1];

    int nmatriz2[60+1][50];*/
}
////////////////////////////////////////////////////////////////////////////////////
int numcelda(int a){
	return(nmatriz[(int)(part[a].x)+1][(int)(part[a].y)+1][(int)(part[a].z-0.5)+1]);
}
////////////////////////////////////////////////////////////////////////////////////
int numcelda2(int a){
    float rn;
    int fraccion;
    rn = sqrt( pow(part[a].x,2) + pow(part[a].y,2) );
    fraccion =  (int)(((int)(rn/(1000000))+1)*pi*0.25);
    printf("\nrn: %f fraccion: %i entrada1: %i entrada2: %i",rn,fraccion,(int)(rn/(1000000))+1,(int)( fraccion*atan(part[a].y/part[a].x)/(0.25*pi) )+1);
    if( fraccion >= 1 ){
        return( nmatriz2[ (int)(rn/(1000000))+1 ][ (int)( fraccion*atan(part[a].y/part[a].x)/(0.25*pi) )+1 ] );
    }
    else{
        return( nmatriz2[1][1] );
    }

}
////////////////////////////////////////////////////////////////////////////////////
void posicion_promedio_celda(int a){
    int i;
    float xcm,ycm,zcm,xsuma=0,ysuma=0,zsuma=0;

    /*if(matriz[a].nparticulas!=0){
        for(i=1;i<=matriz[a].nparticulas;i++){
            xsuma += part[ matriz[a].particulas[i] ].x;
            ysuma += part[ matriz[a].particulas[i] ].y;
            zsuma += part[ matriz[a].particulas[i] ].z;
        }

        xcm = xsuma/(1.0*matriz[a].nparticulas);
        ycm = ysuma/(1.0*matriz[a].nparticulas);
        zcm = zsuma/(1.0*matriz[a].nparticulas);

        matriz[a].x = xcm;
        matriz[a].y = ycm;
        matriz[a].z = zcm;
    }*/
    if(matriz[a].carga!=0){
        for(i=1;i<=matriz[a].nparticulas;i++){
            xsuma += part[ matriz[a].particulas[i] ].carga*part[ matriz[a].particulas[i] ].x;
            ysuma += part[ matriz[a].particulas[i] ].carga*part[ matriz[a].particulas[i] ].y;
            zsuma += part[ matriz[a].particulas[i] ].carga*part[ matriz[a].particulas[i] ].z;
        }

        xcm = xsuma/(1.0*matriz[a].carga);
        ycm = ysuma/(1.0*matriz[a].carga);
        zcm = zsuma/(1.0*matriz[a].carga);

        matriz[a].x = xcm;
        matriz[a].y = ycm;
        matriz[a].z = zcm;
    }
}
////////////////////////////////////////////////////////////////////////////////////
void posiciones_iniciales(){
    int i,j,k,i1;
	int n=0;
    for(i=1;i<w;i++){
        if(n>=npart)break;
        for(j=1;j<w;j++){
            if(n>=npart)break;
            for(k=1;k<=(l-1);k++){
                if(n>=npart)break;
                n++;
                part[n].x = i;
                part[n].y = j;
                part[n].z = k;
                for(i1=1;i1<=150;i1++){
					if(matriz[ numcelda(n) ].particulas[i1]==0){
						matriz[ numcelda(n) ].particulas[i1] = n;
						//printf("\n %i	%i",numcelda(n),matriz[ numcelda(n) ].particulas[i1]);
						//getchar();
						break;
					}
				}

				//printf("\n%i	%i	%i	%i %i",n,(int)(part[n].x),(int)(part[n].y),(int)(part[n].z-0.5)+1,nmatriz[(int)(part[n].x)][(int)(part[n].y)][(int)(part[n].z-0.5)+1]);
				//getchar();

				if(n<=np){
                    part[n].carga = v1;
                }
                else{
                    part[n].carga = v2;
                }

                matriz[ numcelda(n) ].nparticulas++;
                matriz[ numcelda(n) ].carga+=part[n].carga;

            }
        }
    }
    for(i=1;i<=nceldas;i++){
        posicion_promedio_celda(i);
    }
    /*for(i=1;i<=nceldas;i++){
	    if(matriz[i].nparticulas!=0){
			printf("\nCelda: %3.0i Numero de particulas: %i Carga: %2.0i nvecinos : %2.0i Particulas: ",i,matriz[i].nparticulas,matriz[i].carga, matriz[i].nvecinos);
			for(j=1;j<=14;j++){
				if(matriz[i].particulas[j]!=0){
					printf("%i, ",matriz[i].particulas[j]);
				}
			}
		}
	}
	getchar();*/
}
////////////////////////////////////////////////////////////////////////////////////
void elegir_particula(void){
    //if(particulasobrepuesta==0){
    n1 = int(npart*alea())+1;
    if(n1==npart+1)n1=npart;
    part[n1].inmovilidad++;
    //}
    /*printf("\n\n\nSe eligio la particula n1: %i",n1);
    imprimir_particula(n1);
    imprimir_celda(numcelda(n1));*/
}
////////////////////////////////////////////////////////////////////////////////////
void mover(void){
    //if(particulasobrepuesta==0){
    dxx = (alea()-0.5)*dx;
    dyy = (alea()-0.5)*dy;
    dzz = (alea()-0.5)*dz;

    part[0] = part[n1];

    part[n1].x += dxx;//(alea()-0.5)*dx;
    part[n1].y += dyy;//(alea()-0.5)*dy;
    part[n1].z += dzz;//(alea()-0.5)*dz;
    //}
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_periodicas(void){
    int i;

    if(part[n1].x>=w) part[n1].x -= w;
    if(part[n1].x<0) part[n1].x += w;
    if(part[n1].y>=w) part[n1].y -= w;
    if(part[n1].y<0) part[n1].y += w;
    if((part[n1].z<0.5)||(part[n1].z>=(l-0.5))){
        rechazo = 1;
    }
    else{
        matriz[ncvi] = matriz[ numcelda(0) ];
        matriz[ncvf] = matriz[ numcelda(0) ];
        matriz[ncni] = matriz[ numcelda(n1) ];
        matriz[ncnf] = matriz[ numcelda(n1) ];

        if(numcelda(n1)!=numcelda(0)){
            matriz[ncvf].carga-=part[0].carga;
            for(i=1;i<=matriz[ncvf].nparticulas;i++){
                if(matriz[ncvf].particulas[i]==n1){
                    matriz[ncvf].particulas[0] = matriz[ncvf].particulas[i];
                    matriz[ncvf].particulas[i] = matriz[ncvf].particulas[ matriz[ncvf].nparticulas ];
                    matriz[ncvf].particulas[ matriz[ncvf].nparticulas ] = matriz[ncvf].particulas[0];
                }
            }
            matriz[ncvf].nparticulas-=1;

            matriz[ncnf].carga+=part[n1].carga;
            matriz[ncnf].nparticulas+=1;
            matriz[ncnf].particulas[ matriz[ncnf].nparticulas ] = n1;
        }

        posicion_promedio_celda(ncvf);

        matriz[ncvi].x = part[0].x;
        matriz[ncvi].y = part[0].y;
        matriz[ncvi].z = part[0].z;
        matriz[ncnf].x = part[n1].x;
        matriz[ncnf].y = part[n1].y;
        matriz[ncnf].z = part[n1].z;
    }
    /*printf("\n\nSe movio la particula %i a la nueva posicion",n1);
    imprimir_particula(n1);
    imprimir_celda(ncnf);*/
    //getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void condiciones_periodicas2(void){
    int i;
    float r_part;
    r_part = sqrt(pow(part[n1].x,2)+pow(part[n1].y,2));

    if(part[n1].y<0){
        part[n1].x = (part[n1].x-part[n1].y)/sqrt(2);
        part[n1].y = (part[n1].x+part[n1].y)/sqrt(2);
    }
    if(part[n1].y>=part[n1].x){
        part[n1].x = (part[n1].x+part[n1].y)/sqrt(2);
        part[n1].y = (-part[n1].x+part[n1].y)/sqrt(2);
    }

    //if((part[n1].z<0.5)||(part[n1].z>=(l-0.5))){
    if( (r_part > 6000000 )||(part[n1].x<0)||(part[n1].y<-part[n1].x) ){
        rechazo = 1;
    }
    else{
        matriz[ncvi] = matriz[ numcelda(0) ];
        matriz[ncvf] = matriz[ numcelda(0) ];
        matriz[ncni] = matriz[ numcelda(n1) ];
        matriz[ncnf] = matriz[ numcelda(n1) ];

        if(numcelda(n1)!=numcelda(0)){
            matriz[ncvf].carga-=part[0].carga;
            for(i=1;i<=matriz[ncvf].nparticulas;i++){
                if(matriz[ncvf].particulas[i]==n1){
                    matriz[ncvf].particulas[0] = matriz[ncvf].particulas[i];
                    matriz[ncvf].particulas[i] = matriz[ncvf].particulas[ matriz[ncvf].nparticulas ];
                    matriz[ncvf].particulas[ matriz[ncvf].nparticulas ] = matriz[ncvf].particulas[0];
                }
            }
            matriz[ncvf].nparticulas-=1;

            matriz[ncnf].carga+=part[n1].carga;
            matriz[ncnf].nparticulas+=1;
            matriz[ncnf].particulas[ matriz[ncnf].nparticulas ] = n1;
        }

        posicion_promedio_celda(ncvf);

        matriz[ncvi].x = part[0].x;
        matriz[ncvi].y = part[0].y;
        matriz[ncvi].z = part[0].z;
        matriz[ncnf].x = part[n1].x;
        matriz[ncnf].y = part[n1].y;
        matriz[ncnf].z = part[n1].z;
    }
    /*printf("\n\nSe movio la particula %i a la nueva posicion",n1);
    imprimir_particula(n1);
    imprimir_celda(ncnf);*/
    //getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void esferas_duras(void){
    int i,j;
    for(i=1;i<=npart;i++){
        for(j=1;j<=npart;j++){
            if(i>j){
                if(distancia(i,j)<1){

                    /*printf("\nparticulas %i y %i pegadas d: %f paso: %i",i,j,distancia(i,j),p);
                    printf("\ni x: %f y: %f z: %f celda: %i",part[i].x,part[i].y,part[i].z,numcelda(i));
                    printf(" x: %f y: %f z: %f",matriz[numcelda(i)].x,matriz[numcelda(i)].y,matriz[numcelda(i)].z);
                    printf("\ni x: %f y: %f z: %f celda: %i",part[j].x,part[j].y,part[j].z,numcelda(j));
                    printf(" x: %f y: %f z: %f",matriz[numcelda(j)].x,matriz[numcelda(j)].y,matriz[numcelda(j)].z);
                    printf("\n dxx: %f dyy: %f dzz: %f",dxx,dyy,dzz);
                    */
                    printf("\nn1: %i i: %i j: %i d: %f",n1,i,j,distancia(i,j));
                    imprimir_celda( numcelda(i) );
                    imprimir_celda( numcelda(j) );
                    rechazo = 1;
                    //rechazo_esf_mov++;
                    particulasobrepuesta=1;
                }
            }
        }
    }
    if(rechazo==1)getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void esferas_duras_part(int a){
    int i;
    for(i=1;i<=npart;i++){
        if(i!=a){
            if(distancia(a,i)<1){
                printf("\nn1: %i a: %i i: %i d: %f",n1,a,i,distancia(a,i));
                imprimir_celda( numcelda(a) );
                imprimir_celda( numcelda(i) );
            }
        }
    }
    //getchar();
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda(int a){
    int i;

    printf("\n\nCelda: %i x: %f y: %f z: %f carga: %i",a,matriz[a].x,matriz[a].y,matriz[a].z,matriz[a].carga);
    for(i=1;i<=matriz[a].nparticulas;i++){
        imprimir_particula( matriz[a].particulas[i] );
    }
    //printf("\n");

}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_particula(int a){
    printf("\nParticula: %i x: %f y: %f z: %f carga: %i",a,part[a].x,part[a].y,part[a].z,part[a].carga);
}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_celda_f(int a){
    int i;

    fprintf(dat,"\n\nCelda: %i x: %f y: %f z: %f carga: %i npart: %i",a,matriz[a].x,matriz[a].y,matriz[a].z,matriz[a].carga,matriz[a].nparticulas);
    for(i=1;i<=matriz[a].nparticulas;i++){
        imprimir_particula_f( matriz[a].particulas[i] );
    }
    //printf("\n");

}
////////////////////////////////////////////////////////////////////////////////////
void imprimir_particula_f(int a){
    fprintf(dat,"\nParticula: %i x: %f y: %f z: %f carga: %i",a,part[a].x,part[a].y,part[a].z,part[a].carga);
}
////////////////////////////////////////////////////////////////////////////////////
float distancia(int a,int b){
    float dist2, xx, yy;//, zz;

    if(fabs(part[a].x-part[b].x)>(w/2.0)){
        if((part[a].x-part[b].x)>0){
            xx = part[b].x + w;
        }
        else{
            xx = part[b].x - w;
        }
    }
    else
    {
        xx = part[b].x;
    }

    if(fabs(part[a].y-part[b].y)>(w/2.0)){
        if((part[a].y-part[b].y)>0){
            yy = part[b].y + w;
        }
        else{
            yy = part[b].y - w;
        }
    }
    else
    {
        yy = part[b].y;
    }
    /*if(fabs(part[a].z-part[b].z)>((l-1)/2.0)){
        if((part[a].z-part[b].z)>0){
            zz = part[b].z + (l-1);
        }
        else
        {
            zz = part[b].z - (l-1);
        }
    }
    else
    {
        zz = part[b].z;
    }*/

    //dist2 = sqrt(pow(part[a].x-xx,2)+pow(part[a].y-yy,2)+pow(part[a].z-zz,2));
    dist2 = sqrt(pow(part[a].x-xx,2)+pow(part[a].y-yy,2)+pow(part[a].z-part[b].z,2));
    //printf("\na:%i, b:%i, xa:%f, ya:%f, za:%f\nxb:%f, yb:%f, zb:%f, dist:%f",a,b,part[a].x,part[a].y,part[a].z,part[b].x,part[b].y,part[b].z,dist2);
    //getchar();
    return(dist2);
}
////////////////////////////////////////////////////////////////////////////////////
float distanciacelda(int a,int b){
    float dist2, xx, yy;//, zz;

    if(fabs(matriz[a].x-matriz[b].x)>(w/2.0)){
        if((matriz[a].x-matriz[b].x)>0){
            xx = matriz[b].x + w;
        }
        else{
            xx = matriz[b].x - w;
        }
    }
    else
    {
        xx = matriz[b].x;
    }

    if(fabs(matriz[a].y-matriz[b].y)>(w/2.0)){
        if((matriz[a].y-matriz[b].y)>0){
            yy = matriz[b].y + w;
        }
        else{
            yy = matriz[b].y - w;
        }
    }
    else
    {
        yy = matriz[b].y;
    }
    /*if(fabs(matriz[a].z-matriz[b].z)>((l-1)/2.0)){
        if((matriz[a].z-matriz[b].z)>0){
            zz = matriz[b].z + (l-1);
        }
        else
        {
            zz = matriz[b].z - (l-1);
        }
    }
    else
    {
        zz = matriz[b].z;
    }*/

    //dist2 = sqrt(pow(matriz[a].x-xx,2)+pow(matriz[a].y-yy,2)+pow(matriz[a].z-zz,2));
    dist2 = sqrt(pow(matriz[a].x-xx,2)+pow(matriz[a].y-yy,2)+pow(matriz[a].z-matriz[b].z,2));
    //printf("\na:%i, b:%i, xa:%f, ya:%f, za:%f\nxb:%f, yb:%f, zb:%f, dist:%f",a,b,matriz[a].x,matriz[a].y,matriz[a].z,matriz[b].x,matriz[b].y,matriz[b].z,dist2);
    //getchar();
    return(dist2);
}
////////////////////////////////////////////////////////////////////////////////////
float distanciapartcel(int a,int b){
    float dist2, xx, yy;//, zz;

    if(fabs(part[a].x-matriz[b].x)>(w/2.0)){
        if((part[a].x-matriz[b].x)>0){
            xx = matriz[b].x + w;
        }
        else{
            xx = matriz[b].x - w;
        }
    }
    else
    {
        xx = matriz[b].x;
    }

    if(fabs(part[a].y-matriz[b].y)>(w/2.0)){
        if((part[a].y-matriz[b].y)>0){
            yy = matriz[b].y + w;
        }
        else{
            yy = matriz[b].y - w;
        }
    }
    else
    {
        yy = matriz[b].y;
    }
    /*if(fabs(matriz[a].z-matriz[b].z)>((l-1)/2.0)){
        if((matriz[a].z-matriz[b].z)>0){
            zz = matriz[b].z + (l-1);
        }
        else
        {
            zz = matriz[b].z - (l-1);
        }
    }
    else
    {
        zz = matriz[b].z;
    }*/

    //dist2 = sqrt(pow(matriz[a].x-xx,2)+pow(matriz[a].y-yy,2)+pow(matriz[a].z-zz,2));
    dist2 = sqrt(pow(part[a].x-xx,2)+pow(part[a].y-yy,2)+pow(part[a].z-matriz[b].z,2));
    return(dist2);
}
////////////////////////////////////////////////////////////////////////////////////
float distanciaz(int a,int b){
    float distz;

    distz = fabs(part[a].z - part[b].z);

    return(distz);
}
////////////////////////////////////////////////////////////////////////////////////
float distanciazcelda(int a,int b){
    float distz;

    distz = fabs(matriz[a].z - matriz[b].z);

    return(distz);
}
////////////////////////////////////////////////////////////////////////////////////
float distanciazpartcel(int a,int b){
    float distz;

    distz = fabs(part[a].z - matriz[b].z);

    return(distz);
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
int maximo(int a,int b){
    if(a>=b){
        return(a);
    }
    else{
        return(b);
    }
}
////////////////////////////////////////////////////////////////////////////////////
void intercambiar(int a, int b){
    part[0]=part[a];
    part[a]=part[b];
    part[b]=part[0];
}
////////////////////////////////////////////////////////////////////////////////////
float de_mov(void){
    int i,j,k,k1,sonvecinos=0,n1celda,n0celda,vecmax;
    int contadorvec_i=1, contadorvec_f=1;
    float eicoul = 0, efcoul = 0, eisig = 0, efsig = 0, elai = 0, elaf = 0;
    float d, z, r1, r2;
    float dem;

	n1celda = numcelda(n1);
	n0celda = numcelda(0);
	vecmax = maximo(matriz[ n0celda ].nvecinos,matriz[ n1celda ].nvecinos);
	for(i=1;i<=matriz[ n1celda ].nvecinos;i++){
        if(matriz[n1celda].vecinos[i]==n0celda){
            sonvecinos = 1;
        }
	}

    for(j=1;j<=vecmax;j++){
        if(matriz[ n0celda ].vecinos[j]!=0){
            if(matriz[ matriz[ n0celda ].vecinos[j] ].nparticulas!=0){
                for(i=1;i<=matriz[ matriz[ n0celda ].vecinos[j] ].nparticulas;i++){
                    ni = matriz[ matriz[ n0celda ].vecinos[j] ].particulas[i];
                    if(n1!=ni){
                        d = distancia(0,ni);
                        if(d<1.0){
							rechazo = 1;
                            return(0);
                        }
                        z = distanciaz(0,ni);
                        r1 = sqrt(0.5+pow((z)/(w*1.0),2));
                        r2 = sqrt(0.25+pow((z)/(w*1.0),2));

                        eicoul += (qe*qe*part[0].carga*part[ni].carga)/(4*pi*epce*epsi*d*esc);
                        elai += ((-qe*qe*part[0].carga*part[ni].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
                    }
                }
            }
        }
        if(matriz[ n1celda ].vecinos[j]!=0){
            if(matriz[ matriz[ n1celda ].vecinos[j] ].nparticulas!=0){
                for(i=1;i<=matriz[ matriz[ n1celda ].vecinos[j] ].nparticulas;i++){
                    ni = matriz[ matriz[ n1celda ].vecinos[j] ].particulas[i];
                    if(n1!=ni){
                        d = distancia(n1,ni);
                        if(d<1.0){
                            rechazo = 1;
                            return(0);
                        }
                        z = distanciaz(n1,ni);
                        r1 = sqrt(0.5+pow((z)/(w*1.0),2));
                        r2 = sqrt(0.25+pow((z)/(w*1.0),2));

                        efcoul += (qe*qe*part[n1].carga*part[ni].carga)/(4*pi*epce*epsi*d*esc);
                        elaf += ((-qe*qe*part[n1].carga*part[ni].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
                    }
                }
            }
        }
    }
    for(i=1;i<=nceldas;i++){
        if(matriz[ n0celda ].vecinos[ contadorvec_i ]!=i){
            if((matriz[i].carga!=0)&&(i!=n1celda)){
                d = distanciapartcel(0,i);
                if(i == n0celda){
                    /*printf("\nVete a la roña pues v: 0");
                    getchar();*/
                }
                z = distanciazpartcel(0,i);
                r1 = sqrt(0.5+pow((z)/(w*1.0),2));
                r2 = sqrt(0.25+pow((z)/(w*1.0),2));

                //eicoul += (qe*qe*matriz[ncvi].carga*matriz[i].carga)/(4*pi*epce*epsi*d*esc);
                //elai += ((-qe*qe*matriz[ncvi].carga*matriz[i].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
                eicoul += (qe*qe*part[0].carga*matriz[i].carga)/(4*pi*epce*epsi*d*esc);
                elai += ((-qe*qe*part[0].carga*matriz[i].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
            }
        }
        else{
            contadorvec_i++;
        }
        if(matriz[ n1celda ].vecinos[ contadorvec_f ]!=i){
            if((matriz[i].carga!=0)&&(i!=n0celda)){
                d = distanciapartcel(n1,i);
                if(i == n1celda){
                    /*printf("\nVete a la roña pues v: n1");
                    getchar();*/
                }
                z = distanciazpartcel(n1,i);
                r1 = sqrt(0.5+pow((z)/(w*1.0),2));
                r2 = sqrt(0.25+pow((z)/(w*1.0),2));

                //efcoul += (qe*qe*matriz[ncnf].carga*matriz[i].carga)/(4*pi*epce*epsi*d*esc);
                //elaf += ((-qe*qe*matriz[ncnf].carga*matriz[i].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
                efcoul += (qe*qe*part[n1].carga*matriz[i].carga)/(4*pi*epce*epsi*d*esc);
                elaf += ((-qe*qe*part[n1].carga*matriz[i].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
            }
        }
        else{
            contadorvec_f++;
        }
    }
    if(sonvecinos==0){
        d = distanciapartcel(0,ncni);
        z = distanciazpartcel(0,ncni);
        r1 = sqrt(0.5+pow((z)/(w*1.0),2));
        r2 = sqrt(0.25+pow((z)/(w*1.0),2));

        //eicoul += (qe*qe*matriz[ncvi].carga*matriz[ncni].carga)/(4*pi*epce*epsi*d*esc);
        //elai += ((-qe*qe*matriz[ncvi].carga*matriz[ncni].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
        eicoul += (qe*qe*part[0].carga*matriz[ncni].carga)/(4*pi*epce*epsi*d*esc);
        elai += ((-qe*qe*part[0].carga*matriz[ncni].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));

        d = distanciapartcel( n1, ncvf);
        z = distanciazpartcel( n1, ncvf);
        r1 = sqrt(0.5+pow((z)/(w*1.0),2));
        r2 = sqrt(0.25+pow((z)/(w*1.0),2));

        //efcoul += (qe*qe*matriz[ncnf].carga*matriz[ncvf].carga)/(4*pi*epce*epsi*d*esc);
        //elaf += ((-qe*qe*matriz[ncnf].carga*matriz[ncvf].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
        efcoul += (qe*qe*part[n1].carga*matriz[ncvf].carga)/(4*pi*epce*epsi*d*esc);
        elaf += ((-qe*qe*part[n1].carga*matriz[ncvf].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));
    }
    eisig = (part[0].carga*sigma*part[0].z*qe*esc)/(2*epce*epsi);
    efsig = (part[n1].carga*sigma*part[n1].z*qe*esc)/(2*epce*epsi);
    dem = efcoul + efsig + elaf - eicoul - eisig - elai;
    //printf("\n matriz eic: %e elai: %e eis: %e efc: %e elaf: %e efs: %e de: %e",eicoul,elai,eisig,efcoul,elaf,efsig,dem);
    //calcular_de_mov();
    return(dem);
}
////////////////////////////////////////////////////////////////////////////////////
float calcular_de_mov(){
    int i;
    float eicoul = 0, efcoul = 0, eisig = 0, efsig = 0, elai = 0, elaf = 0;
    float d, z, r1, r2;
    float dem;

    for(i=1;i<=npart;i++){
        if(n1!=i){
            d = distancia(0,i);
            if(d<1.0){
                rechazo = 1;
                return(0);
            }

            z = distanciaz(0,i);//fabs(part[0].z-part[i].z);
            r1 = sqrt(0.5+pow((z)/(w*1.0),2));
            r2 = sqrt(0.25+pow((z)/(w*1.0),2));

            eicoul += (qe*qe*part[0].carga*part[i].carga)/(4*pi*epce*epsi*d*esc);
            elai += ((-qe*qe*part[0].carga*part[i].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));

            d = distancia(n1,i);
            if(d<1.0){
                rechazo = 1;
                return(0);
            }
            z = distanciaz(n1,i);
            r1 = sqrt(0.5+pow((z)/(w*1.0),2));
            r2 = sqrt(0.25+pow((z)/(w*1.0),2));


            efcoul += (qe*qe*part[n1].carga*part[i].carga)/(4*pi*epce*epsi*d*esc);
            elaf += ((-qe*qe*part[n1].carga*part[i].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));

        }
    }
    eisig = (part[0].carga*sigma*part[0].z*qe*esc)/(2*epce*epsi);
    efsig = (part[n1].carga*sigma*part[n1].z*qe*esc)/(2*epce*epsi);
    dem = efcoul + efsig + elaf - eicoul - eisig - elai;

    //printf("\nclasica eic: %e elai: %e eis: %e efc: %e elaf: %e efs: %e de: %e",eicoul,elai,eisig,efcoul,elaf,efsig,dem);
    //getchar();

    return(dem);
}
////////////////////////////////////////////////////////////////////////////////////
void metropolis_mov(void){
    float zeta, argexp, emet;
    int i;
    zeta = alea();
    argexp = -de_mov()/(kb*tempe);

    if(argexp>=100)
    {
        emet = 2.0;
    }
    else
    {
        emet = exp(argexp);
    }

    /*if((emet>zeta)&&(rechazo==0)){

        matriz[ numcelda(n1) ] = matriz[ncnf];
        matriz[ numcelda(0) ] = matriz[ncvf];
        c_mova++;
        part[n1].inmovilidad=0;
    }
    else{
        part[n1] = part[0];
        rechazo_met_mov++;
    }*/

    if(rechazo==0){
        if(emet>zeta){
            posicion_promedio_celda(ncnf);
            matriz[ numcelda(n1) ] = matriz[ncnf];
            matriz[ numcelda(0) ] = matriz[ncvf];
            c_mova++;
            part[n1].inmovilidad=0;
            /*printf("\n\nSe acepto el movimiento.");
            imprimir_celda(numcelda(n1));
            imprimir_celda(ncnf);*/
        }
        else{
            if(particulainmovil==1){
                printf("\nNo hubo traslape");
                getchar();
            }
            part[n1] = part[0];
            rechazo_met_mov++;
            /*printf("\n\nSe rechazo el movimiento por el metropolis.");
            imprimir_celda(numcelda(n1));
            imprimir_celda(numcelda(0));*/
        }
    }
    else{
        if(particulainmovil==1){
            printf("\nHubo traslape");
            imprimir_particula(n1);
            esferas_duras_part(n1);
            //getchar();
        }
        rechazo_esf_mov++;
        /*printf("\n\nSe rechazo el movimiento por que hubo un traslape.");
        esferas_duras_part(n1);
        esferas_duras_part(0);
        imprimir_celda(numcelda(n1));
        imprimir_celda(numcelda(0));*/
        part[n1] = part[0];
    }
}

////////////////////////////////////////////////////////////////////////////////////
void gdr(void){
	int i,mp,mn;
	mp=mn=0;
	if(p>terma){
        for(i=1;i<=npart;i++){
            if(part[i].carga > 0 ){
                dgdrpx[(int)(part[i].x/(dw))+1]++;
                dgdrpy[(int)(part[i].y/(dw))+1]++;
                dgdrpz[(int)((part[i].z-0.5)/(dl))+1]++;
                mp++;
            }
            else{
                dgdrnx[(int)(part[i].x/(dw))+1]++;
                dgdrny[(int)(part[i].y/(dw))+1]++;
                dgdrnz[(int)((part[i].z-0.5)/(dl))+1]++;;
                mn++;
            }
        }
        if(mn!=mp-(np-nn)){
            printf("\nmp: %i, mn: %i", mp, mn);
            if(mis==1) printf(" mover");
            if(mis==2) printf(" sustraer");
            if(mis==3) printf(" insertar");
            getchar();
        }
	}
}
////////////////////////////////////////////////////////////////////////////////////
void sumas(void){
    int i;
    sumapx = 0; sumanx = 0; sumapy = 0; sumany = 0; sumapz = 0; sumanz = 0;
    for(i=1;i<=clases;i++){
        sumapx += dgdrpx[i];
        sumanx += dgdrnx[i];
        sumapy += dgdrpy[i];
        sumany += dgdrny[i];
        sumapz += dgdrpz[i];
        sumanz += dgdrnz[i];
    }

}
////////////////////////////////////////////////////////////////////////////////////
void salida(void){
	int i;
	if(((p%actu==0)&&(p>terma))){
		sprintf(salidac,"datos/gdr%i.dat",(p-terma)/actu);
		dat=fopen(salidac,"w");
		for(i=1;i<=clases;i++){
			fprintf(dat,"\n%f	%f	%f	%f	%f	%f	%f	%f	%f",((i-0.5)*dw),
            dgdrpx[i]/((p-terma)*distrmp),dgdrnx[i]/((p-terma)*distrmn),((i-0.5)*dw),
            dgdrpy[i]/((p-terma)*distrmp),dgdrny[i]/((p-terma)*distrmn),(0.5+(i-0.5)*dl),
            dgdrpz[i]/((p-terma)*distrmn),dgdrnz[i]/((p-terma)*distrmn));
            //fprintf(dat,"\n%i %f  %f  %f",i,(dl/2)+(i-1)*dl,dgdrpz[i]/(p*densidadmn),dgdrnz[i]/(p*densidadmn));
		}
		//fprintf(dat,"\nsuma:    %f  %f  %f  %f  %f  %f,%i",sumapx/((p-terma)*1.0),sumanx/((p-terma)*1.0),sumapy/((p-terma)*1.0),sumany/((p-terma)*1.0),sumapz/((p-terma)*1.0),sumanz/((p-terma)*1.0),p);
		//fprintf(dat,"\nsuma: %f, %f",sumapz/(p*1.0),sumanz/(p*1.0));
		fclose(dat);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void salidageo(void){
int m, mp, mn;
mp=mn=0;
sprintf(salidac,"temp/Salida%i.dat",p/actu);
        //printf("\n%s\n",salidac);
dat=fopen(salidac,"w");
        fprintf(dat,"Execute[{");
        for(m=1; m<=npart; m++){
                if(part[m].carga==1){
                        fprintf(dat,"\"A%i:(x-%lf)^2+(y-%lf)^2+(z-%lf)^2=%lf\",\"SetColor[A%i,Blue]\"", m, part[m].x, part[m].y, part[m].z, pow(0.5,2), m);
                        mp++;
                }else if (part[m].carga==-1){
                        fprintf(dat,"\"A%i:(x-%lf)^2+(y-%lf)^2+(z-%lf)^2=%lf\",\"SetColor[A%i,Red]\"", m, part[m].x, part[m].y, part[m].z, pow(0.5,2), m);
                        mn++;
                }
                if(m<npart)fprintf(dat,",");
        }
        fprintf(dat,"}]");
        fclose(dat);
//printf("\n mp: %i mn: %i \n", mp, mn);
}
////////////////////////////////////////////////////////////////////////////////////
float de_insus(int a, int b){
    int i;
    float de, ecoul = 0, esiga = 0, esigb = 0, ela = 0;
    float d, z, r1, r2;

    d = distancia(a,b);
    if(d<1.0){
        //printf("TRASLAPE ENTRE A Y B  a: %i b: %i d: %f",a,b,d);
        //getchar();
        rechazo = 1;
        return(0);
    }

    z = distanciaz(a,b);
    r1 = sqrt(0.5+pow((z)/(w*1.0),2));
    r2 = sqrt(0.25+pow((z)/(w*1.0),2));

    ecoul += (qe*qe*part[a].carga*part[b].carga)/(4*pi*epce*epsi*d*esc);
    ela += ((-qe*qe*part[a].carga*part[b].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));


    for(i=1;i<=npart;i++){
        if(i!=a&&i!=b){
            d = distancia(a,i);
            if(d<1.0){
                //printf("TRASLAPE i: %i a: %i d: %f",i,a,d);
				//getchar();
            	rechazo = 1;
                return(0);
            }

            z = distanciaz(a,i);
            r1 = sqrt(0.5+pow((z)/(w*1.0),2));
            r2 = sqrt(0.25+pow((z)/(w*1.0),2));

            ecoul += (qe*qe*part[a].carga*part[i].carga)/(4*pi*epce*epsi*d*esc);
            ela += ((-qe*qe*part[a].carga*part[i].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));

            d = distancia(b,i);
            if(d<1.0){
                //printf("TRASLAPE i: %i b: %i d: %f",i,b,d);
                //getchar();
                rechazo = 1;
                return(0);
            }
            z = distanciaz(b,i);
            r1 = sqrt(0.5+pow((z)/(w*1.0),2));
            r2 = sqrt(0.25+pow((z)/(w*1.0),2));


            ecoul += (qe*qe*part[b].carga*part[i].carga)/(4*pi*epce*epsi*d*esc);
            ela += ((-qe*qe*part[b].carga*part[i].carga)/(pi*epce*epsi*w*w*esc*esc))*(w*esc*log((0.5+r1)/(r2))+z*esc*atan((4*z*esc*r1)/(w*esc)));

        }
    }

    esiga = (part[a].carga*sigma*part[a].z*qe*esc)/(2*epce*epsi);
    esigb = (part[b].carga*sigma*part[b].z*qe*esc)/(2*epce*epsi);
    de = esiga+esigb+ela+ecoul;
	/*if(p>=000){
			if(c==1){
	        printf("\nSUS ");
	    }
	    else{
	        printf("\nINS ");
	    }
		printf("esiga: %e esigb: %e ecoul: %e de: %e",esiga, esigb, ecoul, de);
		getchar();
	}*/
    return(de);
}
////////////////////////////////////////////////////////////////////////////////////
void sustraer_particulas(){
    n1 = (int)(np*alea()+1);
    if(n1==np+1) n1 = np;
    n2 = (int)(nn*alea()+1);
    if(n2==nn+1) n2 = nn;
    n2+=np;
}
////////////////////////////////////////////////////////////////////////////////////
void metropolis_sustraer(){
    float argexp, emet, factor;

    argexp = -B + de_insus(n1,n2)/(kb*tempe);
    factor = np*nn;

    if(argexp>=100){
        emet = 2.0;
    }
    else{
        emet = factor*exp(argexp);
    }
    if(emet>alea()&&rechazo==0){
        c_susa++;
        intercambiar(n1,npart);
        intercambiar(n1,np);
        intercambiar(n2,npart-1);
        np-=1;
        nn-=1;
        npart=np+nn;
    }
}
////////////////////////////////////////////////////////////////////////////////////
void insertar_particulas(){
    //if(npart+2<=MAXPART){
        part[npart+1].x = alea()*w;
        if(part[npart+1].x==w) part[npart+1].x = 0;
        part[npart+1].y = alea()*w;
        if(part[npart+1].y==w) part[npart+1].y = 0;
        part[npart+1].z = alea()*(l-1)+0.5;
        if(part[npart+1].z==(l-0.5)) part[npart+1].z = 0.5;

        part[npart+2].x = alea()*w;
        if(part[npart+2].x==w) part[npart+2].x = 0;
        part[npart+2].y = alea()*w;
        if(part[npart+2].y==w) part[npart+2].y = 0;
        part[npart+2].z = alea()*(l-1)+0.5;
        if(part[npart+2].z==(l-0.5)) part[npart+2].z = 0.5;

        part[npart+1].carga = v1;
        part[npart+2].carga = v2;
    //}
}
////////////////////////////////////////////////////////////////////////////////////
void metropolis_insertar(){
    float argexp, emet, factor;

    argexp = B - de_insus(npart+1,npart+2)/(kb*tempe);
    factor = 1/((np+1)*(nn+1)*1.0);

    if(argexp>=100){
        emet = 2.0;
    }
    else{
        emet = factor*exp(argexp);
    }
    //printf("\nEXP %e de %e factor %e np %i nn %i",emet,de_insus(npart+1,npart+2,0),factor,np,nn);
    //getchar();
    if(rechazo==0){
	    if(emet>alea()){
	        c_insa++;
			intercambiar(npart+1,np+1);
	        np+=1;
	        nn+=1;
	        npart=np+nn;
	        //printf("\nINSERTAR ACEPTADO");
	        //getchar();
	    }else{
	    	rechazo_met++;
		}
	}
	else{
		rechazo_tras++;
	}
}
////////////////////////////////////////////////////////////////////////////////////
float encuentrab(){
    int i,cont_prom=0,contador_b1=1,contador_b2=2;
    int pr=0,salir=0,i_corte=10,clasei=0.5*clases+1,clasef=0.75*clases;
    int actu_b=100000;
    float db=0.5;
    float sumadgdrpx,sumadgdrnx,sumadgdrpy,sumadgdrny,sumadgdrpz,sumadgdrnz;
    float promediob,sumasg,promedio_b,sumasb;
    float promediog[11]={0},promediobe[11]={0};

    p=1;
    printf("actu_be: %i\n",actu_b);
    sprintf(salidac,"temp/be.dat");
	dat=fopen(salidac,"w");
	fprintf(dat,"B		Bprom		G		Gprom		lngammaprom");
    fclose(dat);

    //for(p=1;p<=50000000;p++){
    while(salir==0){
        rechazo = 0;
        if(alea()<0.8&&npart>0){
            c_mov++;
            elegir_particula();
            mover();
            condiciones_periodicas();
            if(rechazo == 0){
                metropolis_mov();
            }
            else{
                rechazo_per++;
                part[n1] = part[0];
            }
        }
        else{
            if(alea()<0.5&&nn>=1){
                c_sus++;
                sustraer_particulas();
                metropolis_sustraer();
            }
            else{
                c_ins++;
                insertar_particulas();
                metropolis_insertar();
            }
        }
        if(p%1000==0)printf("\r Paso %i",p);
        if(p>terma){
            pr++;
            for(i=1;i<=npart;i++){
                if(part[i].carga > 0 ){
                    dgdrpx[(int)(part[i].x/(dw))+1]++;
                    dgdrpy[(int)(part[i].y/(dw))+1]++;
                    dgdrpz[(int)((part[i].z-0.5)/(dl))+1]++;
                }
                else
                {
                    dgdrnx[(int)(part[i].x/(dw))+1]++;
                    dgdrny[(int)(part[i].y/(dw))+1]++;
                    dgdrnz[(int)((part[i].z-0.5)/(dl))+1]++;;
                }
            }
            if(p%actu_b==0){
                cont_prom++;
                if(cont_prom>i_corte) cont_prom=1;
                contador_b1 = contador_b2;
                sumadgdrpx=sumadgdrnx=sumadgdrpy=sumadgdrny=sumadgdrpz=sumadgdrnz=0;
                //for(i=1;i<=clases;i++){
                //for(i=0.5*clases+1;i<=clases*0.75;i++){
                for(i=clasei;i<=clasef;i++){
                    /*sumadgdrpx += dgdrpx[i];
                    sumadgdrnx += dgdrnx[i];
                    sumadgdrpy += dgdrpy[i];
                    sumadgdrny += dgdrny[i];*/
                    sumadgdrpz += dgdrpz[i];
                    sumadgdrnz += dgdrnz[i];
                }
                //promediog[cont_prom] = (sumadgdrpx+sumadgdrnx+sumadgdrpy+sumadgdrny+sumadgdrpz+sumadgdrnz)/(6.0*distrmn*pr*clases);
                promediog[cont_prom] = (sumadgdrpz+sumadgdrnz)/(2.0*distrmn*pr*(clasef-clasei));
				promediobe[cont_prom] = B;
                sumasg = 0;
                sumasb = 0;
                for(i=1;i<=i_corte;i++){
                    sumasg += promediog[i];
                    sumasb += promediobe[i];
                }
				if(p<=terma+i_corte*actu_b){
                    promediob = sumasg/(1.0*cont_prom);
                    promedio_b = sumasb/(1.0*cont_prom);
				}
				else{
                    promediob = sumasg/(1.0*i_corte);
                    promedio_b = sumasb/(1.0*i_corte);
				}

				sprintf(salidac,"temp/be.dat");
				dat=fopen(salidac,"a");
    			fprintf(dat,"\n%f	%f	%f	%f 	%f",B,promedio_b,promediog[cont_prom],promediob,(promedio_b-log(nnr*nnr))*0.5);

                if(db<0.0001){
                    printf("\nENCONTRO B! B = %f y entonces lngamma = %f",B,(B-log(nnr*nnr))*0.5);
                    getchar();
                    return(B);
                }
                else{
                    if(promediob>1.0){
                        contador_b2 = -1;
                    }
                    else{
                        contador_b2 = 1;
                    }
                    if(contador_b2*contador_b1<0){
                        if(contador_b1!=2)db = db/2;
                    }
                    //printf(" db: %f B: %f\n",contador_b2*db,B);

                    printf("            B: %f db: %f prom: %f prombar: %f \n", B, contador_b2*db, promediog[cont_prom],promediob);
                    B += contador_b2*db;
                    for(i=0;i<=clases;i++){
                        dgdrpx[i]=dgdrpy[i]=dgdrpz[i]=dgdrnx[i]=dgdrny[i]=dgdrnz[i]=0;
                    }
                    pr = 0;
                }
            }
        }
        p++;
        fclose(dat);
    }
}
