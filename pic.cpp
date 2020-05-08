#include <epot_bicgstabsolver.hpp>
#include <particlediagplotter.hpp>
#include <epot_gssolver.hpp>
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "epot_field.hpp"
#include "particledatabase.hpp"
#include "meshvectorfield.hpp"
#include <scalarfield.hpp>
#include "particledatabase.hpp"
#include "geomplotter.hpp"
//#include <gtkplotter.hpp>
#include "ibsimu.hpp"
#include "error.hpp"
#include <iostream>
#include <fstream>
#include "constants.hpp"
#include <time.h>
#include <random.hpp>
#define PI 3.14159265359
#define KB 1.38064852E-23
using namespace std;

//fprintf(fp11," %9.6f ," ,epot(Xce));
float radio = 0.03;
bool cilindro(double x, double y, double z){
	return( ( x*x + y*y ) >= radio*radio  );
}
bool esfera(double x, double y, double z){
	return( ( x*x + y*y + z*z ) >= radio*radio  );
}
bool solenoides(int a, float b){
	/*float z = a*b;
	cout << "a: " << a;
	getchar();*/
	return(a%100<=50);
	//return((z>=1.0&&z<1.2)||(z>=2.0&&z<2.2)||(z>=3.0&&z<3.2)||(z>=4.0&&z<4.2));
	//return(true);
}
double alea(void);
double alea_ab(double a, double b);
float signo(float a);
float distribucion_normal_5(float mu, float sigma);
void prueba (ParticleDataBase3D &pdb);
double calcular_vprom(ParticleDataBase3D *pdb, int a);
void calcular_temperatura(int a, ParticleDataBase3D pdb, double E, Vec3D& temp, ostream& fout);
void crear_haz_gaussiano(double IQ, double q, double m, double l, double E_p, double T_n, int32_t npart, ParticleDataBase3D* pdb);
void colision_pared(ParticleDataBase3D *pdb, size_t k, Particle3D& a);
void scattering_e_0(ParticleDataBase3D *pdb, Particle3D pp, ParticleP3D& a, ostream& fout);
void ionization_e_1(ParticleDataBase3D *pdb, Particle3D pp, ParticleP3D& a, ParticleP3D& b, ParticleP3D& c, ostream& fout, ostream& fout2, int contador);
void scattering_i_0(ParticleDataBase3D *pdb, Particle3D pp, ParticleP3D& a, ostream& fout);
void hacer_rms(int b, ParticleDataBase3D pdb, ofstream &fout);
void leer_archivo_salida(int a, ParticleDataBase3D *pdb, ifstream &fin, ifstream &fin_size);
void imprimir_geometria(Geometry *geo);
void imprimir_pdb(int a, int b, Geometry *geo, EpotField potencial, EpotEfield campoelectrico, ParticleDataBase3D pdb, MeshScalarField tdens);
void hacer_histograma_v_r_E(int a, int b, ParticleDataBase3D pdb);
void calcular_potencial(int a, EpotField potencial);
void calcular_prob_null(int a, ParticleDataBase3D pdb);
void crear_archivo_salida(int a, ParticleDataBase3D pdb, ostream& fout);

void calc_cross_sections(void);
double cross_sec(int a, int b, double E);


//GLOBALES
int32_t cont_vr, cont_scat=0, cont_ref=0, cont_tot=0;
int32_t actuplot = 10;
int32_t npart = 5000, size_arreglo_e = npart;
double E_p = 10000.0, I = 1.0e-1, carga = -1.0, masa = 1.0/1836.0, temp_p = 0.0, temp_n = 4000.0, radio_haz = 0.005, J = -I/(PI*radio_haz*radio_haz);
double Bmax = sqrt( 2*masa*MASS_U*(E_p+temp_p+temp_n)/fabs(carga*CHARGE_E) )/(radio-radio_haz);
//int32_t npart2 = 5000;
//double E_p2 = 50000.0, I2 = 1e-2, carga2 = 1.0, masa2 = 1.00727647, temp_p2 = 0.0, temp_n2 = 20.0, radio_haz2 = 0.005, J2 = -I2/(PI*radio_haz2*radio_haz2);
//double Bmax2 = sqrt( 2*masa2*MASS_U*(E_p2+temp_p2+temp_n2)/fabs(carga2*CHARGE_E) )/(radio-radio_haz2);
int tipo_part[1000000]={0};
double masa_argon = 39.948;
double dt_global = 9e-10;
int contador_ionizadas=0;
Vec3D size(0.08,0.08,0.01), origen(-size[0]/2.0,-size[1]/2.0,0.0);
/*origen[0]=-size[0]/2.0;
origen[1]=-size[1]/2.0;
origen[2]=0.0;*/
int int_size[3];
float cell_size = 0.0025;
float z_max = 10.0;
int32_t niteraciones = z_max/size[2];
double masas[10], cargas[10];
int i_comp, cont_para_ionizacion=0;

double cross_section[5][5], neutral_gas_density = 1e19, prob_null[5]={0}, cs_vmax[5]={0};
int collision_types[5], part_types, contador_scattering_0 = 0, contador_scattering_1 = 0, contador_ionizadas_0 = 0;

//MAIN
int main( int argc, char **argv ){

	srand((unsigned)time(NULL));/*sembrando la semilla*/
	system("mkdir dump");
	system("cp null_collision.cpp dump/");


	/*Geometry g(MODE_3D, Int3D (21, 21, 21), Vec3D (-0.05, -0.05, 0.0), 0.005 );
	ParticleDataBase3D particulas_prueba( g );
	particulas_prueba.add_cylindrical_beam_with_energy( 20000, -0.5, -1.0, 1.0/1836.0, 50000, 0.0, 10.0, Vec3D(0.0, 0.0, 0.0), Vec3D(1.0, 0.0, 0.0), Vec3D(0.0, 1.0, 0), 0.001 );
	scattering2_5(particulas_prueba, 1);*/

	double /*prob_ion, prob_scat = atof(argv[1]), prob_ion, cs_ion = atof(argv[2]), */campoB = atof(argv[1]);
	Vec3D temp;

	/*prob_ion = 1-exp(-size[2]*cs_ion*neutral_gas_density);
	cout << endl << prob_ion << endl;*/
	//E_p = 10; temp_n = 20;

    ifstream fin( "tipos_part.txt" );
    while(!fin.eof()){
        fin>>i_comp;
        fin>>masas[i_comp]>>cargas[i_comp];
        cout << i_comp << " " << masas[i_comp] << " " << cargas[i_comp] << endl;
    }
    part_types = i_comp+1;
    fin.close();
    collision_types[0]=2;
    collision_types[1]=0;
    collision_types[2]=1;

    calc_cross_sections();
	//return(0);

	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );

	int_size[0] = size[0]/cell_size + 1; int_size[1] = size[1]/cell_size + 1; int_size[2] = size[2]/cell_size + 1;
	/*cout << "size[" << 0 << "] " << size[0] << " size[" << 1 << "] " << size[1] << " size[" << 2 << "] " << size[2] << endl;
	cout << "int_size[" << 0 << "] " << int_size[0] << " int_size[" << 1 << "] " << int_size[1] << " int_size[" << 2 << "] " << int_size[2] << endl;
	printf("limite: %f",size[0]-cell_size);
	getchar();*/
	bool fout[3] = {true, true, true};

	Geometry cubo(MODE_3D, Int3D ( int_size[0], int_size[1], int_size[2] ), origen, cell_size );
	Mesh mesh1(MODE_3D, Int3D ( int_size[0], int_size[1], int_size[2] ), origen, cell_size );
	Solid *s2 = new FuncSolid( cilindro );
	//Solid *s2 = new FuncSolid( esfera );
	cubo.set_solid( 7, s2);										//Utilizar siempre el numero mas pequeno mayor a 7 sin usar.
	cubo.set_boundary(7, Bound(BOUND_DIRICHLET, 0.0));
	cubo.build_mesh();

	EpotBiCGSTABSolver solcubo(cubo, 1.0e-4, 10000, 1.0e-4, 10, true);
	//EpotGSSolver solcubo( cubo );

	EpotField potencial(cubo);
	MeshScalarField scharge_beam(cubo);
	MeshScalarField scharge_sec(cubo);
	MeshScalarField scharge_tot(cubo);
	EpotEfield campoelectrico(potencial);
	MeshVectorField campomagnetico(mesh1, fout);
	MeshScalarField tdensi(cubo);
	MeshScalarField tdens_beam(cubo);
	MeshScalarField tdens_sec(cubo);
	MeshScalarField tdens_tot(cubo);
	ParticleDataBase3D electrones(cubo);
	ParticleDataBase3D iones(cubo);
	ParticleDataBase3D particulas(cubo);
	ParticleDataBase3D particulas_sec(cubo);

	particulas.set_max_time( dt_global );
	particulas.set_max_steps( 10000 );
	particulas_sec.set_max_time( dt_global );
	particulas_sec.set_max_steps( 10000 );

	for(int32_t i=0; i<int_size[0];i++){
		for(int32_t j=0; j<int_size[1];j++){
			for(int32_t k=0; k<int_size[2];k++){
				//campomagnetico.set(i, j, k, Vec3D (0.0,0.0,0.0) );
				//campomagnetico.set(i, j, k, Vec3D (0.0, 0.0, Bmax*sqrt( (i-10.0)*(i-10.0) + (j-10.0)*(j-10.0) )/10.0) );
				campomagnetico.set(i, j, k, Vec3D (0.0, 0.0, campoB));
			}
		}
	}

	imprimir_geometria(&cubo);

	ofstream hist_vr( "fraccion_vr.txt" );
	hist_vr << "# fraccion_vr.txt" << endl;
	ofstream ofemitanciax( "emitanciax.txt" );
	ofemitanciax << "# emitanciax.txt alfa\tbeta\tgamma\tepsilon" << endl;
	ofstream ofemitanciay( "emitanciay.txt" );
	ofemitanciay << "# emitanciay.txt alfa\tbeta\tgamma\tepsilon" << endl;

	string filenamevprom = "vprom.txt";
	ofstream vpromedio( filenamevprom.c_str() );
	vpromedio<<"# vprom.txt"<<endl<<"#X\tY"<<endl;

	ofstream rms_elec( "rms_elec.dat" );
	rms_elec <<"#Z X_P Y_P X_M Y_M\n";

	ofstream rms_ion( "rms_ion.dat" );
	rms_ion <<"#Z X_P Y_P X_M Y_M\n";

	ofstream filetemp( "temperatura.txt" );
	filetemp << "# temperatura.txt" << endl;

	crear_haz_gaussiano(I,cargas[0],masas[0],radio_haz,E_p,temp_n,npart,&particulas);
    for(size_t i=0; i<particulas.size(); i++)tipo_part[i]=0;

    //imprimir_pdb(0,&cubo,potencial,campoelectrico,particulas);

    ofstream fout_size("pdb_size.txt");
    fout_size<<0<<" "<<particulas.size()<<endl;
    fout_size<<2<<" "<<particulas_sec.size()<<endl;
    fout_size.close(); fout_size.clear();
    {
        ofstream file_part_out;
        file_part_out.open( "electrones_out.txt" );
        crear_archivo_salida(0,particulas, file_part_out);
        file_part_out.close(); file_part_out.clear();
    }
    {
        ofstream file_part_out;
        file_part_out.open( "iones_out.txt" );
        crear_archivo_salida(2,particulas_sec, file_part_out);
        file_part_out.close(); file_part_out.clear();
    }

    calcular_temperatura(-1,particulas,E_p,temp,filetemp);

    particulas.clear();
	//INICIO
	for(int a=0;a<niteraciones;a++){
	//for(int a=0;a<10;a++){
		cont_vr = 0;
		contador_ionizadas=0;
		//if(a>0)getchar();
		/*for(size_t k=0; k<particulas.size(); k++){
            if(tipo_part[k]!=0){
                Particle3D pp = particulas.particle(k);
                cout << pp.velocity() << pp.location() << endl;
                getchar();
            }
		}*/

		ibsimu.message(1) << "Major cycle " << a << "\n";
		ibsimu.message(1) << "-----------------------\n";
        {
            ifstream fin, fin_size;
            fin.open("electrones_out.txt");
            fin_size.open("pdb_size.txt");
            leer_archivo_salida(0,&particulas,fin,fin_size);
            fin.close(); fin_size.close();
        }
        {
            ifstream fin, fin_size;
            fin.open("iones_out.txt");
            fin_size.open("pdb_size.txt");
            leer_archivo_salida(2,&particulas_sec,fin,fin_size);
            fin.close(); fin_size.close();
        }

        /*cout << "particulas_sec.size: " << particulas_sec.size() << endl;
        for(size_t i=0; i<particulas_sec.size(); i++){
            Particle3D psec = particulas_sec.particle(i);
            for(int j=0; j<7; j++)cout << "psec_" << j << ": " << psec(j) << " ";
            cout << endl;
        }
        getchar();*/

		//if(solenoides(a,size[2])){
		/*nif(solenoides(a,size[2])||1>0){
			for(int32_t i=0; i<int_size[0];i++){
				for(int32_t j=0; j<int_size[1];j++){
					for(int32_t k=0; k<int_size[2];k++){
						//campomagnetico.set(i, j, k, Vec3D (0.0,0.0,0.0) );
						//campomagnetico.set(i, j, k, Vec3D (0.0, 0.0, Bmax*sqrt( (i-10.0)*(i-10.0) + (j-10.0)*(j-10.0) )/10.0) );
						//campomagnetico.set(i, j, k, Vec3D (0.0, 0.0, 10*campoB*exp(-pow((a*size[2]-z_max/2),2))/2*(2*2)));
                        campomagnetico.set(i, j, k, Vec3D (0.0, 0.0, campoB) );
					}
				}
			}
		}else{
			for(int32_t i=0; i<int_size[0];i++){
				for(int32_t j=0; j<int_size[1];j++){
					for(int32_t k=0; k<int_size[2];k++){
						//campomagnetico.set(i, j, k, Vec3D (0.0,0.0,0.0) );
						//campomagnetico.set(i, j, k, Vec3D (0.0, 0.0, Bmax*sqrt( (i-10.0)*(i-10.0) + (j-10.0)*(j-10.0) )/10.0) );
						campomagnetico.set(i, j, k, Vec3D (0.0, 0.0, 0.0));
					}
				}
			}
		}*/

		//densidad.clear();
		solcubo.solve(potencial, scharge_tot);
		/*if( solcubo.get_iter() == 0 && a>0){
			ibsimu.message(1) << "No iterations, breaking major cycle\n";
			break;
		}*/
		campoelectrico.recalculate();

		particulas.iterate_trajectories( scharge_beam, campoelectrico, campomagnetico );
		particulas.build_trajectory_density_field(tdensi);

		particulas_sec.iterate_trajectories( scharge_sec, campoelectrico, campomagnetico );
		particulas_sec.build_trajectory_density_field(tdens_sec);

		hacer_rms(a,particulas,rms_elec);
		hacer_rms(a,particulas_sec,rms_ion);

		scharge_tot = scharge_beam;
		scharge_tot += scharge_sec;

        cout << "Antes del for, size: " << particulas.size() << endl;
        cout << "scat_0: " << contador_scattering_0 << " scat_1: " << contador_scattering_1 << " ion_0: " << contador_ionizadas_0 << endl;
		for( size_t k=0; k<particulas.size(); k++ ){
            Particle3D PPPP_ = particulas.particle( k );
			particle_status_e pstat = PPPP_.get_status();
			if(pstat == PARTICLE_BADDEF ){
			//if(pp(5)>0.0995){//||pp(5)<0.0005){
                cout << endl << "BADDEF" << endl << k << " posicion: " << PPPP_.location() << " rho: " << sqrt( PPPP_[1]*PPPP_[1]+PPPP_[3]*PPPP_[3] ) << " velocidad: " << PPPP_.velocity() << " masa: " << PPPP_.m() << " tipo: " << tipo_part[k] <<endl;

                //Particle3D PPPP_2 = particulas.particle( k + 1 );
                //cout << k + 1 << " posicion: " << PPPP_2.location() << " rho: " << sqrt( PPPP_2[1]*PPPP_2[1]+PPPP_2[3]*PPPP_2[3] ) << " velocidad: " << PPPP_2.velocity() << " masa: " << PPPP_2.m() << PPPP_.m() << " tipo: " << tipo_part[k+1] << endl;

                //PPPP_2 = particulas.particle( k + 2 );
                //cout << k + 2 << " posicion: " << PPPP_2.location() << " rho: " << sqrt( PPPP_2[1]*PPPP_2[1]+PPPP_2[3]*PPPP_2[3] ) << " velocidad: " << PPPP_2.velocity() << " masa: " << PPPP_2.m() << " tipo: " << tipo_part[k+2] << endl;

                //PPPP_2 = particulas.particle( k + 3 );
                //cout << k + 3 << " posicion: " << PPPP_2.location() << " rho: " << sqrt( PPPP_2[1]*PPPP_2[1]+PPPP_2[3]*PPPP_2[3] ) << " velocidad: " << PPPP_2.velocity() << " masa: " << PPPP_2.m() << " tipo: " << tipo_part[k+3] << endl;

                //PPPP_2 = particulas.particle( k + 4 );
                //cout << k + 4 << " posicion: " << PPPP_2.location() << " rho: " << sqrt( PPPP_2[1]*PPPP_2[1]+PPPP_2[3]*PPPP_2[3] ) << " velocidad: " << PPPP_2.velocity() << " masa: " << PPPP_2.m() << " tipo: " << tipo_part[k+4] << endl;
                getchar();
            }
        }
		for( size_t k=0; k<particulas_sec.size(); k++ ){
            Particle3D PPPP_ = particulas_sec.particle( k );
			particle_status_e pstat = PPPP_.get_status();
			if(pstat == PARTICLE_BADDEF ){
			//if(pp(5)>0.0995){//||pp(5)<0.0005){
                cout << endl << "BADDEF_SEC" << endl << k << " posicion: " << PPPP_.location() << " rho: " << sqrt( PPPP_[1]*PPPP_[1]+PPPP_[3]*PPPP_[3] ) << " velocidad: " << PPPP_.velocity() << " masa: " << PPPP_.m() << " tipo: " << tipo_part[size_arreglo_e+k] <<endl;

                //Particle3D PPPP_2 = particulas.particle( k + 1 );
                //cout << k + 1 << " posicion: " << PPPP_2.location() << " rho: " << sqrt( PPPP_2[1]*PPPP_2[1]+PPPP_2[3]*PPPP_2[3] ) << " velocidad: " << PPPP_2.velocity() << " masa: " << PPPP_2.m() << PPPP_.m() << " tipo: " << tipo_part[k+1] << endl;

                //PPPP_2 = particulas.particle( k + 2 );
                //cout << k + 2 << " posicion: " << PPPP_2.location() << " rho: " << sqrt( PPPP_2[1]*PPPP_2[1]+PPPP_2[3]*PPPP_2[3] ) << " velocidad: " << PPPP_2.velocity() << " masa: " << PPPP_2.m() << " tipo: " << tipo_part[k+2] << endl;

                //PPPP_2 = particulas.particle( k + 3 );
                //cout << k + 3 << " posicion: " << PPPP_2.location() << " rho: " << sqrt( PPPP_2[1]*PPPP_2[1]+PPPP_2[3]*PPPP_2[3] ) << " velocidad: " << PPPP_2.velocity() << " masa: " << PPPP_2.m() << " tipo: " << tipo_part[k+3] << endl;

                //PPPP_2 = particulas.particle( k + 4 );
                //cout << k + 4 << " posicion: " << PPPP_2.location() << " rho: " << sqrt( PPPP_2[1]*PPPP_2[1]+PPPP_2[3]*PPPP_2[3] ) << " velocidad: " << PPPP_2.velocity() << " masa: " << PPPP_2.m() << " tipo: " << tipo_part[k+4] << endl;
                getchar();
            }
        }
        cout << "Despues for";

		//particulas.export_path_manager_data("particulas.dat", 50000, 1.0, 1.0/1800.0, Vec3D (0.05, 0.05, 0.0999), Vec3D (1.0, 0.0, 0.0), Vec3D (0.0, 1.0, 0.0) );


		ParticleDiagPlotter pDiagPlotterX(cubo, particulas, AXIS_Z, origen[2]+size[2]-0.1*cell_size, PARTICLE_DIAG_PLOT_SCATTER, DIAG_X, DIAG_XP);
		ParticleDiagPlotter pDiagPlotterY(cubo, particulas, AXIS_Z, origen[2]+size[2]-0.1*cell_size, PARTICLE_DIAG_PLOT_SCATTER, DIAG_Y, DIAG_YP);

		Emittance emitanciax = pDiagPlotterX.calculate_emittance();
		ofemitanciax << a << "\t" << emitanciax.alpha() << "\t" << emitanciax.beta() << "\t" << emitanciax.gamma() << "\t" << emitanciax.epsilon() << "\t" << endl;

		Emittance emitanciay = pDiagPlotterY.calculate_emittance();
		ofemitanciay << a << "\t" << emitanciay.alpha() << "\t" << emitanciay.beta() << "\t" << emitanciay.gamma() << "\t" << emitanciay.epsilon() << "\t" << endl;

		if((a+1)%actuplot==0||a==0){

            imprimir_pdb(a,0,&cubo,potencial,campoelectrico,particulas,tdensi);
            imprimir_pdb(a,1,&cubo,potencial,campoelectrico,particulas_sec,tdens_sec);

			//ParticleDiagPlotter pDiagPlotterY(cubo, particulas, AXIS_Z, 0.000009, PARTICLE_DIAG_PLOT_HISTO2D, DIAG_Y, DIAG_YP);
			ParticleDiagPlotter pDiagPlotterX(cubo, particulas, AXIS_Z, origen[2]+size[2]-0.1*cell_size, PARTICLE_DIAG_PLOT_SCATTER, DIAG_X, DIAG_XP);
			ParticleDiagPlotter pDiagPlotterY(cubo, particulas, AXIS_Z, origen[2]+size[2]-0.1*cell_size, PARTICLE_DIAG_PLOT_SCATTER, DIAG_Y, DIAG_YP);
			ParticleDiagPlotter pDiagPlotterXY(cubo, particulas, AXIS_Z, origen[2]+size[2]-0.1*cell_size, PARTICLE_DIAG_PLOT_SCATTER, DIAG_X, DIAG_Y);

			//pDiagPlotterX.set_ranges(-0.04, -0.05, 0.04, 0.05);
			//pDiagPlotterY.set_ranges(-0.04, -0.05, 0.04, 0.05);
			pDiagPlotterXY.set_ranges(origen[0], origen[1], origen[0]+size[0], origen[1]+size[1]);
			pDiagPlotterX.plot_png( "dump/emittance_beamx"+to_string(1.0*a)+".png");
			pDiagPlotterY.plot_png( "dump/emittance_beamy"+to_string(1.0*a)+".png");
			pDiagPlotterXY.plot_png( "dump/emittance_beamxy"+to_string(1.0*a)+".png");

			hacer_histograma_v_r_E(a, 0, particulas);
			hacer_histograma_v_r_E(a, 2, particulas_sec);

			calcular_potencial(a, potencial);
		}

		cout << endl << "CALCULANDO TEMPERATURA" << endl;
		calcular_temperatura(a, particulas, E_p, temp, filetemp);
		//filetemp << a << "\t";
		//for(int i = 0; i < 3; i++)filetemp << setw(12) << temp(i) << "\t";
		//filetemp << endl;

		/*std::vector<Vec3D> loc_rem;
		std::vector<Vec3D> vel_rem;
   		std::vector<double> charge_rem;
		std::vector<double> mass_rem;

		particle_status_e pstat = pp.get_status();
		if ( (pstat==PARTICLE_OK) || (pstat==PARTICLE_TIME) ) {
			charge_rem.push_back( charge_aux );
			loc_rem.push_back( loc );
			vel_rem.push_back( vel );
			mass_rem.push_back( massnn );
			remainder_particles=remainder_particles + 1;
		}


		for( size_t ae = 0; ae < remainder_particles; ae++ ) {
			double auxcurrent;
			auxcurrent=current;
			if(mass_rem[ae]<0.5){
				auxcurrent=-current;
			}
			pdb.add_particle( auxcurrent, charge_rem[ae],mass_rem[ae],ParticleP3D( 0.0,loc_rem[ae][0], vel_rem[ae][0] ,loc_rem[ae][1],vel_rem[ae][1] ,loc_rem[ae][2],vel_rem[ae][2] ) );
		}*/
		/*for(size_t k = 0; k < particulas_sec.size(); k++){
			Particle3D PP_sec = particulas_sec.particle(k);
			particulas.add_particle(PP_sec);
		}*/
        //CALCULAR PROB_NULL
        calcular_prob_null(0,particulas);
        calcular_prob_null(2,particulas_sec);
		//getchar();
		cout << endl << "CREANDO ARCHIVO DE SALIDA" << endl;
		ParticleP3D PP, PP_elec, PP_ion;

		double r, vr, v_norma, E_pp, aleat;
		ofstream fileout, fileout2;
		fileout.open( "electrones_out.txt" );
		fileout2.open("iones_out.txt");

		for(size_t k = 0; k < particulas.size(); k++){
			Particle3D pp = particulas.particle(k);
			r = sqrt( pp(1)*pp(1) + pp(3)*pp(3) );
			vr = sqrt( pp(2)*pp(2) + pp(4)*pp(4) );
			v_norma = sqrt( pp.velocity()*pp.velocity() );
			E_pp = 0.5*pp.m()*(pp.velocity()*pp.velocity())/CHARGE_E;
			if(vr/pp(6)>0.1)cont_vr++;


			if(r>radio){//-0.01*cell_size){
                //cout << "\npp antes colision\nposicion: " << pp(1) << " " << pp(3) << " " << pp(5) << " r: " << sqrt(pp(1)*pp(1)+pp(3)*pp(3)) ;
                //cout << "\nvelocidad: " << pp(2) << " " << pp(4) << " " << pp(6);
				colision_pared(&particulas, k, pp);
                //cout << "\npp despues colision\nposicion: " << pp(1) << " " << pp(3) << " " << pp(5) << " r: " << sqrt(pp(1)*pp(1)+pp(3)*pp(3)) ;
                //cout << "\nvelocidad: " << pp(2) << " " << pp(4) << " " << pp(6);
                //getchar();
				//cout << "\nparticula: " << k << " r: " << r << " r_nuevo: " << sqrt( pp(1)*pp(1) + pp(3)*pp(3) ) << endl;
				//getchar();
			}
			if(alea()<=prob_null[ 0 ]){
                aleat = alea();
                if(aleat<v_norma*cross_sec(0,0,E_pp)/cs_vmax[0]){
                    contador_scattering_0++;
                    scattering_e_0(&particulas, pp, PP, fileout);
                }else if( ( aleat<v_norma*( cross_sec(0,0,E_pp) + cross_sec(0,1,E_pp) )/cs_vmax[0] ) && E_pp > 11.55 ){
                    contador_ionizadas_0++;
                    cont_para_ionizacion++;
                    ionization_e_1(&particulas, pp, PP, PP_elec, PP_ion, fileout, fileout2, cont_para_ionizacion);
                    if(cont_para_ionizacion%100==0){
                        contador_ionizadas++;
                        cont_para_ionizacion=0;
                        size_arreglo_e++;
                    }
                }else{
                    fileout << "\n";
                    fileout << setw(12) << pp.IQ() << " " << setw(12) << cargas[0] << " " << setw(12) << masas[0] << " ";
                    for(int i = 0; i < 7; i++)fileout << " " << setw(12) << pp(i);
                    fileout << " " << 0;
                }
                //if(aleat<prob_scat){
                /*if(aleat<-1){
                    scattering_e_0(&particulas, k, PP, fileout);
                //}else if(aleat<prob_scat+prob_ion && E_pp > 11.55){
                }else if(aleat<0.001 && E_pp > 11.55){
                    contador_ionizadas++;
                    ionization_e_1(&particulas, k, PP, PP_elec, PP_ion, fileout);
                }else{
                    fileout << "\n";
                    fileout << setw(12) << pp.IQ() << " " << setw(12) << cargas[tipo_part[k]] << " " << setw(12) << masas[tipo_part[k]] << " ";
                    for(int i = 0; i < 7; i++)fileout << " " << setw(12) << pp(i);
                    fileout << " " << tipo_part[k];
                }*/
			}else{
                fileout << "\n";
                fileout << setw(12) << pp.IQ() << " " << setw(12) << cargas[0] << " " << setw(12) << masas[0] << " ";
                for(int i = 0; i < 7; i++)fileout << " " << setw(12) << pp(i);
                fileout << " " << 0;
			}
		}
		for(size_t k = 0; k < particulas_sec.size(); k++){
			Particle3D pp = particulas_sec.particle(k);
			r = sqrt( pp(1)*pp(1) + pp(3)*pp(3) );
			vr = sqrt( pp(2)*pp(2) + pp(4)*pp(4) );
			v_norma = sqrt( pp.velocity()*pp.velocity() );
			E_pp = 0.5*pp.m()*(pp.velocity()*pp.velocity())/CHARGE_E;
			if(vr/pp(6)>0.1)cont_vr++;


			if(r>radio){//-0.01*cell_size){
                //cout << "\npp antes colision SEC\nposicion: " << pp(1) << " " << pp(3) << " " << pp(5) << " r: " << sqrt(pp(1)*pp(1)+pp(3)*pp(3)) ;
                //cout << "\nvelocidad: " << pp(2) << " " << pp(4) << " " << pp(6);
				colision_pared(&particulas_sec, k, pp);
                //cout << "\npp despues colision SEC\nposicion: " << pp(1) << " " << pp(3) << " " << pp(5) << " r: " << sqrt(pp(1)*pp(1)+pp(3)*pp(3)) ;
                //cout << "\nvelocidad: " << pp(2) << " " << pp(4) << " " << pp(6);
                //getchar();
				//cout << "\nparticula: " << k << " r: " << r << " r_nuevo: " << sqrt( pp(1)*pp(1) + pp(3)*pp(3) ) << endl;
				//getchar();
			}
			if(alea()<=prob_null[ 2 ]){
                aleat = alea();
                if(aleat<v_norma*cross_sec(1,0,E_pp)/cs_vmax[1]){
                    contador_scattering_1++;
                    scattering_i_0(&particulas_sec, pp, PP, fileout2);
                }else{
                    fileout2 << "\n";
                    fileout2 << setw(12) << pp.IQ() << " " << setw(12) << cargas[2] << " " << setw(12) << masas[2] << " ";
                    for(int i = 0; i < 7; i++)fileout2 << " " << setw(12) << pp(i);
                    fileout2 << " " << 2;
                }
                //if(aleat<prob_scat){
                /*if(aleat<-1){
                    scattering_e_0(&particulas, k, PP, fileout);
                //}else if(aleat<prob_scat+prob_ion && E_pp > 11.55){
                }else if(aleat<0.001 && E_pp > 11.55){
                    contador_ionizadas++;
                    ionization_e_1(&particulas, k, PP, PP_elec, PP_ion, fileout);
                }else{
                    fileout << "\n";
                    fileout << setw(12) << pp.IQ() << " " << setw(12) << cargas[tipo_part[k]] << " " << setw(12) << masas[tipo_part[k]] << " ";
                    for(int i = 0; i < 7; i++)fileout << " " << setw(12) << pp(i);
                    fileout << " " << tipo_part[k];
                }*/
			}else{
                fileout2 << "\n";
                fileout2 << setw(12) << pp.IQ() << " " << setw(12) << cargas[2] << " " << setw(12) << masas[2] << " ";
                for(int i = 0; i < 7; i++)fileout2 << " " << setw(12) << pp(i);
                fileout2 << " " << 2;
			}
		}
		fileout.close(); fileout.clear();
		fileout2.close(); fileout2.clear();
		{
		    ofstream fout_size("pdb_size.txt");
            fout_size << 0 << " " << particulas.size()+contador_ionizadas<< endl;
            fout_size << 2 << " " << particulas_sec.size()+contador_ionizadas << endl;
		}

		particulas.clear();
		particulas_sec.clear();
		if(a%actuplot==0){
			vpromedio << a << "\t" << calcular_vprom(&particulas, a) << endl;
		}
		hist_vr << a << "\t" <<1.0*cont_vr/(particulas.size())<< "\t" << particulas.size() << endl;
		cout << particulas.size() << endl;
	}
	//FINAL
	filetemp.close();
	hist_vr.close();
	ofemitanciax.close();
	ofemitanciay.close();
	vpromedio.close();
	rms_elec.close();
	rms_ion.close();
	cout << "reflejadas: " << cont_ref << " dispersadas: " << cont_scat << " totales: " << cont_tot << endl;
	cout << endl << "ionizadas: " << contador_ionizadas << endl;
	//system("free -m");
	return( 0 );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
double alea(void){//genera numeros aleatorios entre 0 y 1
return((double)rand()/RAND_MAX); //32767 dependiente de la libreria
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
double alea_ab(double a, double b){//genera numeros aleatorios entre a y b
return( a + (double)rand()/RAND_MAX*(b-a) ); //32767 dependiente de la libreria
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
float distribucion_normal_5(float mu, float sigma){
	float b, fb;
	do{
		b = alea_ab(-5*sigma + mu, 5*sigma + mu);
		fb = exp( -(b-mu)*(b-mu)/(2*sigma*sigma) )/(sqrt(2*PI)*sigma);
	}while(alea()/(sqrt(2*PI)*sigma)>fb);
	return(b);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void colision_pared(ParticleDataBase3D *pdb, size_t k, Particle3D& a){
	Particle3D pp = pdb->particle(k);
	Vec3D vel = pp.velocity(), loc = pp.location();
	double theta = atan2( pp(3), pp(1) );
	a(0) = pp(0);
	a(1) = pp(1)-1.0*signo(pp(1))*cell_size;
	a(2) = -pp(4)*sin(2*theta) - pp(2)*cos(2*theta);
	a(3) = pp(3)-1.0*signo(pp(3))*cell_size;
	a(4) = -pp(2)*sin(2*theta) + pp(4)*cos(2*theta);
	a(5) = pp(5);
	a(6) = pp(6);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void scattering_e_0(ParticleDataBase3D *pdb, Particle3D pp, ParticleP3D& a, ostream& fout){
	//Particle3D pp = pdb->particle(k);
	double ang_chi, ang_phi, ang_theta,  v_norma, E_inc;
	Vec3D vinc_i, vinc_j, vinc_k;
	v_norma = sqrt(pp.velocity()*pp.velocity());

	E_inc = 0.5*pp.m()*(pp.velocity()*pp.velocity())/CHARGE_E;									//ENERGIA EN ELECTRON-VOLTS
	ang_chi = acos( (2.0 + E_inc - 2.0*pow((1.0+E_inc), alea()))/E_inc );
	//ang_chi = acos(1-2.0*alea());																//Angulo para E muy pequeno
	ang_phi = 2.0*PI*alea();

	vinc_i = pp.velocity();
	vinc_i.normalize();
	vinc_j = cross( vinc_i, Vec3D (1.0, 0.0, 0.0) );
	vinc_k = cross( vinc_i, -vinc_j );
	ang_theta = acos( vinc_i*Vec3D (1.0, 0.0, 0.0) );
	v_norma *= sqrt(1+(2*masa/masa_argon)*(1-cos(ang_chi)));

	a(0) = pp(0);
	a(1) = pp(1);
	a(2) = v_norma*(vinc_i(0)*cos(ang_chi) + vinc_j(0)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(0)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));
	a(3) = pp(3);
	a(4) = v_norma*(vinc_i(1)*cos(ang_chi) + vinc_j(1)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(1)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));
	a(5) = pp(5);
	a(6) = v_norma*(vinc_i(2)*cos(ang_chi) + vinc_j(2)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(2)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));

    fout << "\n";
    fout << setw(12) << pp.IQ() << " " << setw(12) << cargas[0] << " " << setw(12) << masas[0] << " ";
    for(int i = 0; i < 7; i++)fout << " " << setw(12) << a(i);
    fout << " " << 0;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void scattering_i_0(ParticleDataBase3D *pdb, Particle3D pp, ParticleP3D& a, ostream& fout){
	//Particle3D pp = pdb->particle(k);
	double ang_chi, ang_phi, ang_theta,  v_norma, E_inc, E_scat;
	Vec3D vinc_i, vinc_j, vinc_k;
	v_norma = sqrt(pp.velocity()*pp.velocity());

	E_inc = 0.5*pp.m()*(pp.velocity()*pp.velocity())/CHARGE_E;									//ENERGIA EN ELECTRON-VOLTS
	ang_chi = sqrt(1-alea());
	E_scat = E_inc*cos(ang_chi)*cos(ang_chi);
	//ang_chi = acos(1-2.0*alea());																//Angulo para E muy pequeno
	ang_phi = 2.0*PI*alea();

	vinc_i = pp.velocity();
	vinc_i.normalize();
	vinc_j = cross( vinc_i, Vec3D (1.0, 0.0, 0.0) );
	vinc_k = cross( vinc_i, -vinc_j );
	ang_theta = acos( vinc_i*Vec3D (1.0, 0.0, 0.0) );
	v_norma = sqrt(2*E_scat*CHARGE_E/pp.m());
	//v_norma *= sqrt(1+(2*masa/masa_argon)*(1-cos(ang_chi)));

	a(0) = pp(0);
	a(1) = pp(1);
	a(2) = v_norma*(vinc_i(0)*cos(ang_chi) + vinc_j(0)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(0)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));
	a(3) = pp(3);
	a(4) = v_norma*(vinc_i(1)*cos(ang_chi) + vinc_j(1)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(1)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));
	a(5) = pp(5);
	a(6) = v_norma*(vinc_i(2)*cos(ang_chi) + vinc_j(2)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(2)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));

    fout << "\n";
    fout << setw(12) << pp.IQ() << " " << setw(12) << cargas[2] << " " << setw(12) << masas[2] << " ";
    for(int i = 0; i < 7; i++)fout << " " << setw(12) << a(i);
    fout << " " << 2;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void ionization_e_1(ParticleDataBase3D *pdb, Particle3D pp, ParticleP3D& a, ParticleP3D& b, ParticleP3D& c, ostream& fout, ostream& fout2, int contador){

	//Particle3D pp = pdb->particle(k);
	double ang_chi, ang_phi, ang_theta, v_norma, v_enorma, E_inc;
	Vec3D vinc_i, vinc_j, vinc_k;
	v_norma = sqrt(pp.velocity()*pp.velocity());
	double z_alea = alea_ab(0.01*cell_size,size[2]-0.01*cell_size);

	E_inc = 0.5*pp.m()*(pp.velocity()*pp.velocity())/CHARGE_E;									//ENERGIA EN ELECTRON-VOLTS
	ang_chi = acos( (2.0 + E_inc - 2.0*pow((1.0+E_inc), alea()))/E_inc );
	//ang_chi = acos(1-2.0*alea());																//Angulo para E muy pequeno
	ang_phi = 2.0*PI*alea();
	vinc_i = pp.velocity();
	vinc_i.normalize();
	vinc_j = cross( vinc_i, Vec3D (1.0, 0.0, 0.0) );
	vinc_k = cross( vinc_i, -vinc_j );
	ang_theta = acos( vinc_i*Vec3D (1.0, 0.0, 0.0) );
	v_norma *= sqrt(1+(2*masa/masa_argon)*(1-cos(ang_chi)));

	a(0) = pp(0);
	a(1) = pp(1);
	a(2) = v_norma*(vinc_i(0)*cos(ang_chi) + vinc_j(0)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(0)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));
	a(3) = pp(3);
	a(4) = v_norma*(vinc_i(1)*cos(ang_chi) + vinc_j(1)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(1)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));
	a(5) = pp(5);
	a(6) = v_norma*(vinc_i(2)*cos(ang_chi) + vinc_j(2)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(2)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));


	/*if(sqrt( a(2)*a(2) + a(4)*a(4) + a(6)*a(6))>1e10){

		cout << "Particula  pp super luminica k: " << k << endl;
		cout << "\n";
		for(int i = 0; i < 7; i++)cout << " " << setw(12) << a(i);
		cout << endl << "vr: " << sqrt( a(2)*a(2) + a(4)*a(4) + a(6)*a(6));
		getchar();
	}*/


	double B_einc=10, E_ion = 11.55, E_elec;													//En electronvolts
	double aleatt = alea();

	E_inc = 0.5*pp.m()*(pp.velocity()*pp.velocity())/CHARGE_E;									//ENERGIA EN ELECTRON-VOLTS
	E_elec = B_einc*tan( aleatt*atan( (E_inc-E_ion)/(2*B_einc) ) );								//ENERGIA EN ELECTRON-VOLTS
	v_enorma = sqrt( (2*E_elec*CHARGE_E)/(pp.m()) );
	/*if( (2*E_elec*CHARGE_E)/(pp.m()) < 0){

		cout << endl << "VELOCIDAD IMAGINARIA X)" << endl;
		cout << "aleatt: " << aleatt << " B_einc: " << B_einc << " E_inc: " << E_inc << " E_ion: " << E_ion << endl;

		getchar();

	}*/

	ang_chi = acos( (2.0 + E_inc - 2.0*pow((1.0+E_inc), alea()))/E_inc );
	//ang_chi = acos(1-2.0*alea());																//Angulo para E muy pequeno
	ang_phi = 2.0*PI*alea();

	vinc_i = pp.velocity();
	vinc_i.normalize();
	vinc_j = cross( vinc_i, Vec3D (1.0, 0.0, 0.0) );
	vinc_k = cross( vinc_i, -vinc_j );
	ang_theta = acos( vinc_i*Vec3D (1.0, 0.0, 0.0) );

	b(0) = pp(0);
	b(1) = pp(1);
	b(2) = v_enorma*(vinc_i(0)*cos(ang_chi) + vinc_j(0)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(0)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));
	b(3) = pp(3);
	b(4) = v_enorma*(vinc_i(1)*cos(ang_chi) + vinc_j(1)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(1)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));
	b(5) = z_alea;
	b(6) = v_enorma*(vinc_i(2)*cos(ang_chi) + vinc_j(2)*(sin(ang_chi)*sin(ang_phi)/sin(ang_theta)) + vinc_k(2)*(sin(ang_chi)*cos(ang_phi)/sin(ang_theta)));

	/*if(sqrt( b(2)*b(2) + b(4)*b(4) + b(6)*b(6))>1e10){

		cout << "Particula  pp_elec super luminica k: " << k << endl;
		cout << "\n";
		for(int i = 0; i < 7; i++)cout << " " << setw(12) << b(i);
		cout << endl << "vr: " << sqrt( b(2)*b(2) + b(4)*b(4) + b(6)*b(6));
		getchar();
	}*/

	double /*vi, fvi, temp_neutro = 300, */sigma=0, energia_iones=0.0;
	//sigma = sqrt((KB*temp_neutro)/(masas[2]*MASS_U));
	sigma = sqrt( 2.0*1000.0*CHARGE_E/(3.0*masas[2]*MASS_U));
	c(0) = pp(0);
	c(1) = pp(1);
	/*do{
		vi = (alea()-0.5)*10.0*sigma;
		fvi = exp(-(vi*vi)/(2.0*sigma*sigma));
	}while(fvi<alea());
	c(2) = vi;*/
	c(2) = distribucion_normal_5(0,sigma);
	c(3) = pp(3);
	/*do{
		vi = (alea()-0.5)*10.0*sigma;
		fvi = exp(-(vi*vi)/(2.0*sigma*sigma));
	}while(fvi<alea());
	c(4) = vi;*/
	c(4) = distribucion_normal_5(0,sigma);
	c(5) = z_alea;
	/*do{
		vi = (alea()-0.5)*10.0*sigma;
		fvi = exp(-(vi*vi)/(2.0*sigma*sigma));
	}while(fvi<alea());
	c(6) = vi;*/
	c(6) = distribucion_normal_5(0,sigma);
	energia_iones = 0.5*masas[2]*MASS_U*(c(2)*c(2) + c(4)*c(4) + c(6)*c(6))/CHARGE_E;
	//energia_iones = 0.0;
	if(energia_iones<-100.0){
		cout << "energia ion: " << energia_iones << endl;
		getchar();
	}
	/*if(sqrt( c(2)*c(2) + c(4)*c(4) + c(6)*c(6))>1e10){

		cout << "Particula  pp_ion super luminica k: " << k << endl;
		cout << "\n";
		for(int i = 0; i < 7; i++)cout << " " << setw(12) << c(i);
		cout << endl << "vr: " << sqrt( c(2)*c(2) + c(4)*c(4) + c(6)*c(6));
		getchar();
	}*/
	if(contador%100==0){
        fout << "\n";
        fout << setw(12) << -fabs(pp.IQ())*100 << " " << setw(12) << cargas[0] << " " << setw(12) << masas[0] << " ";
        for(int i = 0; i < 7; i++)fout << " " << setw(12) << b(i);
        fout << " " << 0;
        fout2 << "\n";
        fout2 << setw(12) << fabs(pp.IQ())*100 << " " << setw(12) << cargas[2] << " " << setw(12) << masas[2] << " ";
        for(int i = 0; i < 7; i++)fout2 << " " << setw(12) << c(i);
        fout2 << " " << 2;
	}
    fout << "\n";
    fout << setw(12) << pp.IQ() << " " << setw(12) << cargas[0] << " " << setw(12) << masas[0] << " ";
    for(int i = 0; i < 7; i++)fout << " " << setw(12) << a(i);
    fout << " " << 0;

    /*cout << "\nPARTICULA A";
    cout << "\n";
    cout << setw(12) << pp.IQ() << " " << setw(12) << cargas[tipo_part[k]] << " " << setw(12) << masas[tipo_part[k]] << " ";
    for(int i = 0; i < 7; i++)cout << " " << setw(12) << a(i);
    cout << " " << tipo_part[k];
    cout << " " << 2;
    cout << "\nPARTICULA B";
    cout << "\n";
    cout << setw(12) << -fabs(pp.IQ()) << " " << setw(12) << cargas[0] << " " << setw(12) << masas[0] << " ";
    for(int i = 0; i < 7; i++)cout << " " << setw(12) << b(i);
    cout << " " << 0;
    cout << "\nPARTICULA C";
    cout << "\n";
    cout << setw(12) << fabs(pp.IQ()) << " " << setw(12) << cargas[2] << " " << setw(12) << masas[2] << " ";
    for(int i = 0; i < 7; i++)cout << " " << setw(12) << c(i);
    getchar();*/
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void crear_haz_gaussiano(double IQ, double q, double m, double r, double E_p, double T_n, int32_t npart, ParticleDataBase3D* pdb){
double fxi, xi, sigma = r/5.0;
double x, vx, y, vy, z, vz;
int clases=200;
double ancho = 2.0*r/clases;
int32_t gdx[5001]={0}, gdy[5001]={0};

for(int i=1;i<=npart;i++){
	do{
		xi = 2.0*(alea()-0.5)*r;
		fxi = (1.0/(sigma*sqrt(2.0*PI)))*(exp(-(xi*xi)/(2.0*sigma*sigma)));
	}while(fxi<alea()*(1.0/(sigma*sqrt(2.0*PI))));
	x = xi;
	do{
		xi = 2.0*(alea()-0.5)*r;
		fxi = (1.0/(sigma*sqrt(2.0*PI)))*(exp(-(xi*xi)/(2.0*sigma*sigma)));
	}while(fxi<alea()*(1.0/(sigma*sqrt(2.0*PI))));
	y = xi;
	z = 0.0;
	vx = distribucion_normal_5(0,sqrt((T_n*CHARGE_E)/(m*MASS_U) ) );
	vy = distribucion_normal_5(0,sqrt((T_n*CHARGE_E)/(m*MASS_U) ) );
	vz = sqrt( 2.0*(E_p*CHARGE_E)/(m*MASS_U) );
	gdx[(int)((x+r)/(ancho))+1]++;
	gdy[(int)((y+r)/(ancho))+1]++;
	//histo<<setw(12)<<x<<"\t"<<y<<endl;

	pdb->add_particle( IQ/npart, q, m, ParticleP3D(0,x,vx,y,vy,z,vz) );
}
ofstream histox( "histogramax.txt" );
histox<<"# histogramax.txt"<<endl<<"#X\tY"<<endl;
for(int i=1;i<=clases;i++){
histox<<-r+i*ancho<<"\t"<<gdx[i]<<endl;
}
histox.close();

/*ofstream histoy( "histogramay.txt" );
histoy<<"# histogramay.txt"<<endl<<"#X\tY"<<endl;
for(int i=1;i<=clases;i++){
	histoy<<-0.04+i*ancho<<"\t"<<gdy[i]<<endl;
}
histoy.close();
getchar();*/
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void prueba (ParticleDataBase3D &pdb){
Particle3D pp = pdb.particle( 0 );
Particle3D pp2 = pdb.particle( 1 );
/*Vec3D vel1=pp.velocity(), vel2=pp2.velocity();
cout << "Dentro de prueba bclear size: " << pdb.size() << endl;
cout << "vel1: " << vel1 << " vel2: " << vel2 << endl;
	cout << " Dentro de prueba bclear 00 " << pp(0) << " " << pp(1) << " " << pp(2) << " " << pp(3) << " " << pp(4) << " " << pp(5) << " " << pp(6) << endl;
	cout << " Dentro de prueba bclear 11 " << pp2(0) << " " << pp2(1) << " " << pp2(2) << " " << pp2(3) << " " << pp2(4) << " " << pp2(5) << " " << pp2(6) << endl;*/

pdb.clear();

/*Particle3D pp1 = pdb.particle( 0 );
Particle3D pp3 = pdb.particle( 1 );
vel1=pp1.velocity(), vel2=pp3.velocity();
	cout << "Dentro de prueba aclear size: " << pdb.size() << endl;
	cout << "vel1: " << vel1 << " vel2: " << vel2 << endl;
	cout << " Dentro de prueba aclear 00 " << pp1(0) << " " << pp1(1) << " " << pp1(2) << " " << pp1(3) << " " << pp1(4) << " " << pp1(5) << " " << pp1(6) << endl;
	cout << " Dentro de prueba aclear 11 " << pp3(0) << " " << pp3(1) << " " << pp3(2) << " " << pp3(3) << " " << pp3(4) << " " << pp3(5) << " " << pp3(6) << endl;
	getchar();*/

pdb.add_particle( pp.IQ(), pp.q(), pp.m(), ParticleP3D(pp(0),pp(1),pp(2),pp(3),pp(4),pp(5),100) );
pdb.add_particle( pp2.IQ(), pp2.q(), pp2.m(), ParticleP3D(pp2(0),pp2(1),pp2(2),pp2(3),pp2(4),pp2(5),200) );

}
/////////////////////////////////////////////////////////////////////////////////////////////////////FALTA LA DENSIDAD
double calcular_vprom(ParticleDataBase3D *pdb, int a){
	double vprom = 0;
	Vec3D vel;
	for(size_t k = 0; k < pdb->size(); k++){
		Particle3D pp = pdb->particle(k);
		vel = pp.velocity();
		vprom += vel*vel;
	}
	vprom = vprom/(pdb->size()+1.0);
	return(vprom);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_temperatura(int a, ParticleDataBase3D pdb, double E, Vec3D& temp, ostream& fout){
double vref, masa;
double sum_x=0, sum_y=0, sum_z=0, I_total=0;
int contadorparticulas=0;
for(size_t k = 0; k < pdb.size(); k++){
	Particle3D pp = pdb.particle(k);
	if(pp.m()<2.0*MASS_U){
		Vec3D vel = pp.velocity();
		vref = sqrt( (2*E*fabs(pp.q()))/(pp.m()) );

		sum_x += fabs(pp.IQ())*vel(0)*vel(0);
		sum_y += fabs(pp.IQ())*vel(1)*vel(1);
		sum_z += fabs(pp.IQ())*(vel(2)-vref)*(vel(2)-vref);

		masa = pp.m();
		I_total += fabs(pp.IQ());
		contadorparticulas++;
	}
}
temp(0) = (masa*sum_x)/(CHARGE_E*I_total);
temp(1) = (masa*sum_y)/(CHARGE_E*I_total);
temp(2) = (masa*sum_z)/(CHARGE_E*I_total);

if(contadorparticulas==0){
	cout << "Que shows v:" << endl;
	getchar();
}


fout << a << "\t";
for(int i = 0; i < 3; i++)fout << setw(12) << temp(i) << "\t";
fout << endl;

/*cout << endl << "Temp_x: " << temp(0) << endl;
cout << "Temp_y: " << temp(1) << endl;
cout << "Temp_z: " << temp(2) << endl;*/
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void hacer_rms(int b, ParticleDataBase3D pdb, ofstream& fout){
	double rms_x=0, rms_y=0, ms_x=0, ms_y=0;
	float contador = 0;
	for(size_t k = 0; k<pdb.size(); k++){
		Particle3D pp = pdb.particle(k);
        ms_x += pp.IQ()*pp(1)*pp(1);
        ms_y += pp.IQ()*pp(3)*pp(3);
        contador += pp.IQ();
	}
	/*for(size_t k = 0; k<10; k++){
		ms_x += k*k;
		ms_y += 4*k*k;
	}*/
	rms_x = sqrt(ms_x/contador);
	rms_y = sqrt(ms_y/contador);
	fout << b*size[2] << " " << rms_x << " " << rms_y << endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void calc_cross_sections(void){
    for(int i=0; i<part_types; i++){
        for(int j=0; j<collision_types[i]; j++){
            cross_section[i][j]=cross_sec(i,j,1);
        }
    }
    for(int i=0; i<part_types; i++){
        for(int j=0; j<collision_types[i]; j++){
            cout << "cs_" << i << "_" << j << ": " << cross_section[i][j] << "\t";
        }
        cout << endl;
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
double cross_sec(int a, int b, double E){
    if(a==0){
        if(b==0){
            return(1e-19);
        }else if(b==1){
            return(2e-20);
        }else if(b==2){
            return(1e-20);
        }else{
            cout << "ALGO SALIO MAL, CS() A";
            return(0);
        }
    }else if(a==1){
        if(b==0){
            return(6.5e-17);
        }else if(b==1){
            return(1e-20);
        }else if(b==2){
            return(1e-20);
        }else{
            cout << "ALGO SALIO MAL, CS() B";
            return(0);
        }
    }else if(a==2){
        if(b==0){
            return(1e-20);
        }else if(b==1){
            return(1e-20);
        }else if(b==2){
            return(1e-20);
        }else{
            cout << "ALGO SALIO MAL, CS() C";
            return(0);
        }
    }
    cout << "ALGO SALIO MAL, CS()";
    return(0);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
float signo(float a){
    return(a>0?1:-1);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void leer_archivo_salida(int a, ParticleDataBase3D *pdb, ifstream& fin, ifstream& fin_size ){
    cout << endl << "LEYENDO ARCHIVO DE SALIDA" << endl;
    double IQ, m, q, t, x, vx, y, vy, z, vz;
    int tipo=100, size_pdb;
    do{
        fin_size >> tipo >> size_pdb;
    }while(tipo!=a);

    if(size_pdb != 0){
        //Leer las particulas
        //ifstream fin( "electrones_out.txt" );
        //while(!fin.eof()){
        for(int i=0; i<size_pdb; i++){
            fin>>IQ>>q>>m>>t>>x>>vx>>y>>vy>>z>>vz>>tipo;
            if(vx>1e9||vy>1e9||vz>1e9){
                cout << "\nvx: " << vx << " vy: " << vy << " vz: " << vz;
                getchar();
            }
            if(IQ!=0||q!=0||m!=0||t!=0||x!=0||vx!=0||y!=0||vy!=0||z!=0||vz!=0){
                if(z>size[2]*0.995){//-cell_size){
                    if(vz>0){
                        pdb->add_particle( IQ, q, m, ParticleP3D(0.0,x,vx,y,vy,z-(size[2]*0.995),vz) );
                    }else{
                        pdb->add_particle( IQ, q, m, ParticleP3D(0.0,x,vx,y,vy,z-size[2]*0.005,vz) );
                    }
                }else if(z<size[2]*0.005){//cell_size){
                    if(vz<0){
                        pdb->add_particle( IQ, q, m, ParticleP3D(0.0,x,vx,y,vy,z+(size[2]*0.995),vz) );
                    }else{
                        pdb->add_particle( IQ, q, m, ParticleP3D(0.0,x,vx,y,vy,z+size[2]*0.005,vz) );
                    }
                }else{
                    pdb->add_particle( IQ, q, m, ParticleP3D(0.0,x,vx,y,vy,z,vz) );
                }
            }
            //}
            if(pdb->size()>0){
                Particle3D pp = pdb->particle( pdb->size()-1 );
                if(pp(5)>size[2]*0.995){//>size[2]-cell_size){
                    if(pp(6)>0){
                        cout << "MAYORMAYOR\n" << pp.location() << pp.velocity() << " " << pdb->size()-1 << "\t" << vz << endl;
                        getchar();
                    }
                }
                if(pp(5)<size[2]*0.005){//cell_size){
                    if(pp(6)<0){
                        cout << "MENORMENOR\t" << pp.location() << pp.velocity() << " " << pdb->size()-1 << "\t" << vz << endl;
                        getchar();
                    }
                }
            }
        }
    }
    cout << "a: " << a << " pdb.size: " << pdb->size();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void imprimir_geometria(Geometry *geo){

	GeomPlotter geomplotter2( *geo );
	geomplotter2.set_size( 750, 750 );
	string nombre2 = "dump/geoxy.png";
	geomplotter2.set_view(VIEW_XY);
	geomplotter2.plot_png( nombre2 );
	nombre2 = "dump/geozx.png";
	geomplotter2.set_view(VIEW_ZX);
	geomplotter2.plot_png( nombre2 );

}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void imprimir_pdb(int a, int b, Geometry *geo, EpotField potencial, EpotEfield campoelectrico, ParticleDataBase3D pdb, MeshScalarField tdens){

    GeomPlotter geomplotter( *geo );
    geomplotter.set_size( 750, 750 );
    geomplotter.set_epot( &potencial );
    geomplotter.set_efield( &campoelectrico );
    geomplotter.set_particle_database( &pdb );
    geomplotter.set_particle_div(pdb.size()/100);

    string nombre;
    if(b==0) nombre = "dump/plot_eleczx"+to_string(1.0*a)+".png";
    else nombre = "dump/plot_ionzx"+to_string(1.0*a)+".png";
    geomplotter.set_view(VIEW_ZX);
    geomplotter.plot_png( nombre );

    if(b==0) nombre = "dump/plot_eleczy"+to_string(1.0*a)+".png";
    else nombre = "dump/plot_ionzy"+to_string(1.0*a)+".png";
    geomplotter.set_view(VIEW_ZY);
    geomplotter.plot_png( nombre );

    if(b==0) nombre = "dump/plot_elecxy"+to_string(1.0*a)+".png";
    else nombre = "dump/plot_ionxy"+to_string(1.0*a)+".png";
    geomplotter.set_view(VIEW_XY);
    geomplotter.plot_png( nombre );


    GeomPlotter geomplotter2(*geo);
    geomplotter2.set_size( 750, 750 );
    geomplotter2.set_efield( &campoelectrico );
    geomplotter2.set_trajdens( &tdens );
    geomplotter2.set_fieldgraph_plot( FIELD_TRAJDENS );
    if(b==0) nombre = "dump/tdens_beam"+to_string(1.0*a)+".png";
    else nombre = "dump/tdens_sec"+to_string(1.0*a)+".png";
    geomplotter2.set_view(VIEW_ZX);
    geomplotter2.plot_png( nombre );

}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void hacer_histograma_v_r_E(int a, int b, ParticleDataBase3D pdb){

    if(pdb.size()>1){
        cout << endl << "HACIENDO HISTOGRAMAS" << endl;
        //if(a>2)getchar();
        int clases=500, clasesr = 200;
        float vxmin, vxmax, vymin, vymax, rmin=0, rmax=0.05, E_min=0, E_max;
        double anchovx, anchovy, anchor, anchoE;
        int32_t gdvx[5001]={0}, gdvy[5001]={0};
        double I_total=0, I_total_E = 0, gdr[5001]={0}, gdE[5001]={0};

        string filename;
        Particle3D pp = pdb.particle( 0 );
        Vec3D vel = pp.velocity();
        double r_i, E_i;
        vxmax = vel[0]; vxmin = vel[0];
        vymax = vel[1]; vymin = vel[1];
        for( size_t i = 0; i < pdb.size(); i++ ) {
            Particle3D pp = pdb.particle( i );
            vel = pp.velocity();
            if(b==0){
                E_i = 0.5*pp.m()*(vel*vel)/CHARGE_E;
                if(E_i>E_max)E_max=E_i;
            }
            if(vel[0]<vxmin)vxmin=vel[0];
            if(vel[0]>vxmax)vxmax=vel[0];
            if(vel[1]<vymin)vymin=vel[1];
            if(vel[1]>vymax)vymax=vel[1];
        }
        anchovx = (vxmax - vxmin)/clases;
        anchovy = (vymax - vymin)/clases;
        anchor = (rmax - rmin)/clasesr;
        anchoE = (E_max - E_min)/clasesr;

        for( size_t i = 0; i < pdb.size(); i++ ) {
            Particle3D pp = pdb.particle( i );
            vel = pp.velocity();
            r_i = sqrt( pp(1)*pp(1)+pp(3)*pp(3) );
            //cout << "ANTES " << r_i << "\t" << (int)((r_i-rmin)/(anchor))+1;//fabs( pp.IQ() ) << "\t" << I_total[ tipo_part[i] ] << "\t" << gdr[(int)((r_i-rmin)/(anchovy))+1][ tipo_part[i] ];
            if(( (int)((vel[0]-vxmin)/(anchovx))+1 )<0 ){
                cout << "\nint_0: " << ( (int)((vel[0]-vxmin)/(anchovx))+1 ) <<  " vel_0: " << vel[0] << " vxmin: " << vxmin << " ancho: " << anchovx << endl;
                getchar();
            }
            if( ( (int)((vel[1]-vymin)/(anchovy))+1 )<0 ){
                cout << "\nint_1: " << ( (int)((vel[1]-vymin)/(anchovy))+1 ) <<  " vel_1: " << vel[1] << " vymin: " << vymin << " ancho: " << anchovy << endl;
                getchar();
            }
            if( ( (int)((r_i-rmin)/(anchor))+1)<0){
                cout << "\nint_2: " << ( (int)((r_i-rmin)/(anchor))+1) <<  " r_i: " << r_i << " rmin: " << rmin << " anchor: " << anchor << endl;
                getchar();
            }
            gdvx[(int)((vel[0]-vxmin)/(anchovx))+1]++;
            gdvy[(int)((vel[1]-vymin)/(anchovy))+1]++;
            gdr[(int)((r_i-rmin)/(anchor))+1]+=fabs( pp.IQ() );
            I_total += fabs( pp.IQ() );
            if(b==0){
                E_i = 0.5*pp.m()*(vel*vel)/CHARGE_E;
                gdE[(int)((E_i-E_min)/(anchoE))+1] += fabs( pp.IQ() );
                I_total_E += fabs( pp.IQ() );
            }
            //cout << "\nDESPUES " << fabs( pp.IQ() ) << "\t" << I_total[ tipo_part[i] ] << "\t" << gdr[(int)((r_i-rmin)/(anchovy))+1][ tipo_part[i] ];
            //getchar();
        }

        if(b==0) filename = "dump/histv_elec"+to_string(1.0*a)+".txt";
        else filename = "dump/histv_ion"+to_string(1.0*a)+".txt";
        ofstream histo( filename.c_str() );

        if(b==0) filename = "dump/dens_elec"+to_string(1.0*a)+".dat";
        else filename = "dump/dens_ion"+to_string(1.0*a)+".dat";
        ofstream densidad( filename.c_str() );

        if(b==0){
            filename = "dump/eedf"+to_string(1.0*a)+".dat";
            ofstream eedf( filename.c_str() );
            for(int i=1;i<=clasesr;i++){
                eedf<<E_min+i*anchoE<<"\t"<<gdE[i]/I_total<<endl;
            }
            eedf.close();
        }
        histo<<"# histograma.txt"<<endl<<"#X\tY"<<endl;
        for(int i=1;i<=clases;i++){
            histo<<vxmin+i*anchovx<<"\t"<<gdvx[i]<<"\t"<<vymin+i*anchovy<<"\t"<<gdvy[i]<<endl;
        }
        for(int i=1;i<=clasesr;i++){
            densidad<<rmin+i*anchor<<"\t"<<gdr[i]/I_total_E<<endl;
        }
        histo.close();densidad.close();
    }else{
        cout << "\nPdb.size(): " << pdb.size();
    }

}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_potencial(int a, EpotField potencial){

    cout << endl << "CREANDO ARCHIVO POTENCIAL" << endl;
    Vec3D r_epot;
    //double ancho_z = size[2]/100;
    double ancho_x = size[0]/100;
    string filename = "dump/epot"+to_string(1.0*a)+".txt";
    ofstream epotfile( filename.c_str() );
    epotfile << "\n";
    for(int k = 0; k < 100; k++){
        r_epot(0) = k*ancho_x;
        r_epot(1) = 0;
        r_epot(2) = 0;//k*ancho_z;
        epotfile << r_epot(0) << "\t" << potencial(r_epot) << endl;
    }
    epotfile.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void calcular_prob_null(int a, ParticleDataBase3D pdb){
    double cs_v[5]={0}, E_eV, v_norma;// = 1-exp(-neutral_gas_density*size[2]*(cross_section[0][0]+cross_section[0][1]+cross_section[1][0]));

    for(size_t k=0; k<pdb.size(); k++){
        Particle3D part_null = pdb.particle(k);
        v_norma = sqrt( part_null.velocity()*part_null.velocity() );
        E_eV = 0.5*part_null.m()*v_norma*v_norma/CHARGE_E;
        cs_v[ a ]=0;

        for(int i=0; i<collision_types[ a ]; i++ )cs_v[ a ] += v_norma*cross_sec(a,i,E_eV);
        if(cs_v[ a ]>cs_vmax[ a ] ) cs_vmax[ a ] = cs_v[ a ];
        //for(int i=0; i<collision_types[ tipo_part[k] ]; i++ )cout << " cs_" << tipo_part[k] << "_" << i << ": " << cross_sec(tipo_part[k],i,E_eV);
        //getchar();
    }
    cout << endl << "dt*n*cs_vmax_" << a << ": " << dt_global*neutral_gas_density*cs_vmax[a];
    prob_null[a] = 1 - exp( - dt_global*neutral_gas_density*cs_vmax[a] );
    cout << endl << "prob_null_" << a << ": " << prob_null[a];

}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void crear_archivo_salida(int a, ParticleDataBase3D pdb, ostream& fout){
    for(size_t i=0; i<pdb.size(); i++){
        Particle3D pout=pdb.particle(i);
        fout << "\n";
        fout << setw(12) << pout.IQ() << " " << setw(12) << cargas[a] << " " << setw(12) << masas[a] << " ";
        for(int j = 0; j < 7; j++)fout << " " << setw(12) << pout(j);
        fout << " " << a;
    }
    cout << "\nparticulas.size: " << pdb.size();
}
