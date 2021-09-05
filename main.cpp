//#include"stdafx.h"

#include"search_for.h"
#include <vector>
#include<stdio.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <fstream>

using namespace std;
int str_to_int(string str)
{
	return stoi(str);
}
long long int str_to_long_long(string str)
{
	return stoll(str);
}
double str_to_d(string str)
{
	return stod(str);
}
double ourrandom()
{
	return rand() / (RAND_MAX + 1.);
}
void read_file(string infilename, int &L, double &J, double &alpha, double &beta, int &num, int &outNum, double &Jthird, double &T, long long int &nsamples,
	int &seed, string &outfilename)
{
	string str_ret;
	bool found;

	search_for(string("L"), infilename, str_ret, found);
	if (found) { L = str_to_int(str_ret); }
	else { L = 4; }

	search_for(string("num"), infilename, str_ret, found);
	if (found) { num = str_to_int(str_ret); }
	else { num = 100; }

	search_for(string("outNum"), infilename, str_ret, found);
	if (found) { outNum = str_to_int(str_ret); }
	else { outNum = 0; }

	//cout<<str_ret<<L<<endl;
	search_for(string("J"), infilename, str_ret, found);
	if (found) { J = str_to_d(str_ret); }
	else { J = 0.0; }

	search_for(string("alpha"), infilename, str_ret, found);
	if (found) { alpha = str_to_d(str_ret); }
	else { alpha = 0.0; }
	
	search_for(string("beta"), infilename, str_ret, found);
	if (found) { beta = str_to_d(str_ret); }
	else { beta = 0.0; }

	//cout<<"J"<<str_ret<<J<<endl;
	search_for(string("Jthird"), infilename, str_ret, found);
	if (found) { Jthird = str_to_d(str_ret); }
	else { Jthird = 0.0; }
	//cout<<"Jthird"<<str_ret<<Jthird<<endl;
	search_for(string("T"), infilename, str_ret, found);
	if (found) { T = str_to_d(str_ret); }
	else { cout << "T not found" << endl; T = 0.0; }

	search_for(string("nsamples"), infilename, str_ret, found);
	if (found) { nsamples = str_to_long_long(str_ret); }
	else { nsamples = 1000; }


	search_for(string("seed"), infilename, str_ret, found);
	if (found) { seed = str_to_int(str_ret); }
	else { seed = 0; }


	search_for(string("outfilename"), infilename, outfilename, found);
	if (!found) { outfilename = "C.txt"; }
}

int randint(int nsites)
{
	return rand() % nsites;
}

vector<vector<double>> create_config(int L)
{
	vector<vector<double>> config(L*L);

	for (int i = 0; i < L*L; i++) {
		double x1, x2;
		double new_x, new_y, new_z;
		while (true) {
			x1 = (ourrandom() - 0.5) *2.0;
			x2 = (ourrandom() - 0.5) *2.0;
			if (x1*x1 + x2*x2 < 1.0) break;
		}
		new_x = 2.0 * x1 * sqrt(1.0 - x1*x1 - x2*x2);
		new_y = 2.0 * x2 * sqrt(1.0 - x1*x1 - x2*x2);
		new_z = 1.0 - 2.0 * (1.0 - x1*x1 - x2*x2);
		config[i] = { new_x,new_y,new_z };
	} 
	return config;
}
int Right(int i,int L){return (i + 1) % L + i / L*L;}
int Left(int i, int L){return (i - 1 + L) % L + i / L*L;}
int TopLeft(int i, int L) { return (i + L) % (L*L); }
int TopRight(int i,int L){return (TopLeft(i,L) + 1) % L + TopLeft(i,L) / L*L;}
int TopRightRight(int i, int L){	return (TopLeft(i,L) + 2) % L + TopLeft(i,L) / L*L;}
int BotRight(int i, int L) { return (i - L + L*L) % (L*L); }
int BotRightRight(int i, int L){ return (BotRight(i,L) + 1) % L + BotRight(i,L) / L*L; }
int TopTopLeft(int i, int L){ return(i + 2 * L) % (L*L); }
int TopTop(int i, int L){ return(TopTopLeft(i,L) + 1) % L + TopTopLeft(i,L) / L*L; }

void make_nn(int &L, vector< vector<int> > &neighbors)
{
	vector<int>    nnentry;

	//////////////////////////////////////////////
	//             Triangular lattice
	//////////////////////////////////////////////
	// Make sites on a quadrilateral with length L on both sides.

	for (int i = 0; i < L*L; i++) {

		///////////////////////////////////////////////////////////////////////////////////////////////
		int right = (i + 1) % L+i/L*L;
		int left = (i - 1 + L) % L+i/L*L;
		int topLeft = (i + L) % (L*L);
		int topRight = (topLeft + 1) % L+topLeft/L*L;
		int botRight = (i - L + L*L) % (L*L);
		int botLeft = (botRight - 1 + L) % L+botRight/L*L;
		nnentry.push_back(right);
		nnentry.push_back(left);
		nnentry.push_back(topLeft);
		nnentry.push_back(topRight);
		nnentry.push_back(botRight);
		nnentry.push_back(botLeft);
		neighbors.push_back(nnentry);
		nnentry.clear();
	}


}

void make_3rdn(int &L, vector< vector<int> > &third_neighbors)
{
	vector<int>    nnentry;

	for (int i = 0; i < L*L; i++) {

		///////////////////////////////////////////////////////////////////////////////////////////////
		int right = (i + 2) % L + i / L*L;
		int left = (i - 2 + L) % L + i / L*L;
		int topLeft = (i + 2*L) % (L*L);
		int topRight = (topLeft + 2) % L + topLeft / L*L;
		int botRight = (i - 2*L + L*L) % (L*L);
		int botLeft = (botRight - 2 + L) % L + botRight / L*L;
		nnentry.push_back(right);
		nnentry.push_back(left);
		nnentry.push_back(topLeft);
		nnentry.push_back(topRight);
		nnentry.push_back(botLeft);
		nnentry.push_back(botRight);
		third_neighbors.push_back(nnentry);
		nnentry.clear();
	}


}


vector<vector<double>> disorder(int&L, int&num, int &nsites, double &beta, vector<vector<int>> &neighbors, vector<vector<vector<double>>> &e)
{
	vector<vector<double>> nn = { { 1,0 },{ -1,0 },{ -0.5,sqrt(3)*0.5 },{ 0.5,sqrt(3)*0.5 },{ 0.5,-sqrt(3)*0.5 },{ -0.5,-sqrt(3)*0.5 } };

	vector<vector<double>> result(nsites, { 5,5,5,5,5,5 });  //Initalize as 5 
	for (int i = 0; i < nsites; i++) {  //There is a problem. dj is not symmetric. 
		for (int j = 0; j < 6; j++) {
			double dj = 0;
			int i_i = i / L;
			int i_j = i % L;
			double exx = e[i_i][i_j][0];
			double eyy = e[i_i][i_j][1];
			double exy = e[i_i][i_j][2];
			double nx = nn[j][0];
			double ny = nn[j][1];
			dj = dj + exx * nx*nx + eyy * ny*ny + 2 * exy*nx*ny;//this is dl;
			dj = sqrt(beta) * dj;
            if (dj > 1){dj = 1;}
            else if (dj < -1){dj = -1;}
			result[i][j] = dj;
			int n = neighbors[i][j];
			for (int k = 0; k < 6; k++){
				if (neighbors[n][k] == i){
					if (result[n][k] > 4){continue;}
					else {
						result[i][j] = (result[n][k] + dj) * 0.5 ;
						result[n][k] = result[i][j] ; 
						}
				}
			}
		}
	}
	return result;
}
double nnenergy(int&nsites, int&L, double &J, vector<vector<double>> &config, vector< vector<int> > &neighbors, vector< vector<double> > &djm)
{


	double energycalc = 0.0;

	for (int i = 0; i<nsites; i++)
	{

		for (int k = 0; k < neighbors[i].size(); k++)
		{
			//cout <<i<<"	"<< neighbors[i][k] << endl;
			energycalc = energycalc + (djm[i][k] + J)*(double(config[i][0]) * double(config[neighbors[i][k]][0]) +
				double(config[i][1]) * double(config[neighbors[i][k]][1]) +
				double(config[i][2]) * double(config[neighbors[i][k]][2]));
			//cout<< double(config[neighbors[i][k]])<<endl;		
		}
	}
	return energycalc / (2.0);
}
double third_energy(int&nsites, int&L, double &Jthird, vector<vector<double>> &config, vector< vector<int> > &third_neighbors)
{


	double energycalc = 0.0;

	for (int i = 0; i<nsites; i++)
	{

		for (int k = 0; k < third_neighbors[i].size(); k++)
		{
			energycalc = energycalc + (Jthird)*(double(config[i][0]) * double(config[third_neighbors[i][k]][0]) +
				double(config[i][1]) * double(config[third_neighbors[i][k]][1]) +
				double(config[i][2]) * double(config[third_neighbors[i][k]][2]));

			//cout<< double(config[neighbors[i][k]])<<endl;		
		}
	}
	return energycalc / (2.0);
}

double diffenergy_nn(double &J, vector<vector<double>> &config,
	vector< vector<int> > &neighbors, int &chosensite, vector<double> new_old, vector< vector<double> > &djm)
{


	double energycalc = 0.0;
	for (int k = 0; k < neighbors[chosensite].size(); k++)
	{
		energycalc = energycalc + (djm[chosensite][k] + J)*(double(new_old[0]) * double(config[neighbors[chosensite][k]][0]) +
			double(new_old[1]) * double(config[neighbors[chosensite][k]][1]) +
			double(new_old[2]) * double(config[neighbors[chosensite][k]][2]));
	}
	return energycalc;
}
double diffenergy_third(double &Jthird, vector<vector<double>> &config,
	vector< vector<int> > &third_neighbors, int &chosensite, vector<double> new_old)
{


	double energycalc = 0.0;
	for (int k = 0; k < third_neighbors[chosensite].size(); k++)
	{
		energycalc = energycalc + Jthird*(double(new_old[0]) * double(config[third_neighbors[chosensite][k]][0]) +
			double(new_old[1]) * double(config[third_neighbors[chosensite][k]][1]) +
			double(new_old[2]) * double(config[third_neighbors[chosensite][k]][2]));
	}
	return energycalc;
}

vector<double> rotate(const vector<double> &v, const vector<double> &n, double &theta) {


    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    vector<double> cross = {n[1] * v[2] - n[2] * v[1], n[2] * v[0] - n[0] * v[2], n[0] * v[1] - n[1] * v[0]};
    double dot = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];

    vector<double> term1 = {v[0] * cos_theta, v[1] * cos_theta, v[2] * cos_theta};
    vector<double> term2 = {cross[0] * sin_theta, cross[1] * sin_theta, cross[2] * sin_theta};
    vector<double> term3 = {n[0] * dot * (1 - cos_theta), n[1] * dot * (1 - cos_theta), n[2] * dot * (1 - cos_theta)};

    vector<double> rotated = {term1[0] + term2[0] + term3[0], term1[1] + term2[1] + term3[1],
                              term1[2] + term2[2] + term3[2]};


    return rotated;
}
vector<double> cone_update(double cone_size, vector<double> &old)
{

    double A = cone_size;
    double cos_theta = (1.0 - cos(A)) * ourrandom() + cos(A);
    double theta_rot = acos(cos_theta);

    double x = old[0];
    double y = old[1];
    double z = old[2];
    double theta_vector = acos(z);
    double phi_vector = atan2(y, x);

    double phi_rot = 2.0 * M_PI * ourrandom();
    double theta_new = theta_vector - theta_rot;
    double xn = sin(theta_new) * cos(phi_vector);
    double yn = sin(theta_new) * sin(phi_vector);
    double zn = cos(theta_new);
    vector<double> new_vector = {xn, yn, zn};
    //vector<double> v_rot = new_vector*cos(phi_rot) + cross(old, new_vector)*sin(phi_rot) + old*dot(old, new_vector)*(1-std::cos(phi_rot));
    // if (isnan(new_vector[0] )) {cout<<z<<" "<<acos(z)<<endl;}

    return rotate(new_vector, old, phi_rot);
}


int main(int argc, char *argv[])
{

	double J = 0.0;
	double Jthird = 0.0;
	double alpha = 0.0;
	double beta = 0.0;
	long long int nsamples = 10;

	double T;
	int seed;
	int L;
	int num;
	int outNum;


	string outfilename;
	string infilename;
	string sample;
	infilename = argv[1];
	sample = argv[2];
	read_file(infilename, L, J,alpha, beta ,num, outNum, Jthird, T, nsamples, seed, outfilename);

	vector<vector<vector<double> > > e(L, vector<vector<double>>(L, vector<double>(3)));  //matrix of {exx, eyy, exy}
	ifstream file(sample);   
	
	if (! file.is_open()){
		cout<<"sample can't open"<<endl;
		return 1; 
		}

	for (int k = 0; k < 3; k++) {
		for (int i = 0; i < L; i++) {
			for (int j = 0; j < L; j++) {
				file >> e[i][j][k];
			}
		}
	}
	file.close(); 
	std::srand(seed);
//	cout<<L<<J<<alpha<<Jthird<<T<<endl;
	nsamples = (long long int)120000 * (long long int)(L*L);
	//nsamples = 10;
	long long int nwates =nsamples/12*2; //(long long int)200000 * (long long int)(L*L);
	int nsites = L*L;
	int pointNumber = 0;
	
	T = T*Jthird;
	double final_t=T;
	double t=0.5*Jthird;

	std::vector< std::vector<int> > neighbors;

	std::vector< std::vector<int> > third_neighbors;

	make_nn(L, neighbors);

	make_3rdn(L, third_neighbors);

	double zero = 0.0;
	vector<vector<double>> djm = disorder(L,num, nsites, zero, neighbors, e);
	ofstream out;
    string name = outfilename+"bias_field.txt";
	out.open(name);
	for (int i = 0; i < nsites; i++) { out << i << "=" << djm[i][0] << "=" << djm[i][1] << "=" << djm[i][2] << "=" << djm[i][3] << "=" << djm[i][4] << "=" << djm[i][5] << endl; }
	out.close();
	std::srand(seed);

	ofstream Eoutfile;
	name = outfilename + ".txt";
	Eoutfile.open(name);


	ofstream outfile[outNum];

	double energySum = 0.0;
	double energySquareSum = 0.0;
	cout<<outNum<<endl;

	//// Create configuration and measure energy

	vector<vector<double>> config = create_config(L);
	//for (int i=0;i<nsites;i++) cout<<config[i];

	//cout<<config.size()<<endl;
	//	cout << "start Metropolis" << endl;
	//	cout<<nsamples<<endl;
	// Run Metropolis
	double energycurrent = nnenergy(nsites, L, J, config, neighbors,djm) + third_energy(nsites, L, Jthird, config, third_neighbors);
	//cout << energycurrent/nsites << endl;
	
	for (long long int n = 0; n < nsamples; n++)

	{	/*if (n == 0) { T = t; }
		if (n == nsamples / 20-1) { T = (t-final_t)*0.8+final_t; }
		if (n == nsamples / 20*2-1) { T = (t-final_t)*0.7+final_t; }
		if (n == nsamples / 20*3-1) { T = (t-final_t)*0.6+final_t; }
		if (n == nsamples / 20*4-1) { T = (t-final_t)*0.5+final_t; }
		if (n == nsamples / 20*5-1) { T = (t-final_t)*0.4+final_t; }
		if (n == nsamples / 20*6-1) { T = (t-final_t)*0.3+final_t; }
		if (n == nsamples / 20*7-1) { T = (t-final_t)*0.2+final_t; }
		if (n == nsamples / 20*8-1) { T = (t-final_t)*0.1+final_t; }
		if (n == nsamples / 20*9-1) { T = final_t; }
	  //	if (n == nsamples ) { T = t*0.02; }*/
		if (n == nwates) {
			cout<<"add disorder."<<endl;
			djm = disorder(L,num, nsites,beta, neighbors, e);}
			cout<<energycurrent<<endl;
			energycurrent = nnenergy(nsites, L, J, config, neighbors,djm) + third_energy(nsites, L, Jthird, config, third_neighbors);
			cout<<energycurrent<<endl;


		int chosensite = randint(nsites);

		double x1, x2;
		double new_x, new_y, new_z;
		while (true) {
			x1 = (ourrandom() - 0.5) *2.0;
			x2 = (ourrandom() - 0.5) *2.0;
			if (x1*x1 + x2*x2 < 1.0) break;
		}
		new_x = 2.0 * x1 * sqrt(1.0 - x1*x1 - x2*x2);
		new_y = 2.0 * x2 * sqrt(1.0 - x1*x1 - x2*x2);
		new_z = 1.0 - 2.0 * (1.0 - x1*x1 - x2*x2);

		/*
        vector<double> new_spin = cone_update(M_PI/5.0, config[chosensite]);
        double norm = sqrt(new_spin[0]*new_spin[0]+new_spin[1]*new_spin[1]+new_spin[2]*new_spin[2]);
        double new_x = new_spin[0]/norm;
        double new_y = new_spin[1]/norm;
        double new_z = new_spin[2]/norm;
		*/

		vector<double >new_old =
		{ new_x - config[chosensite][0], new_y - config[chosensite][1], new_z - config[chosensite][2] };

		double energydiff = diffenergy_nn(J, config, neighbors, chosensite, new_old,djm) +
			diffenergy_third(Jthird, config, third_neighbors, chosensite, new_old);
	    //	cout << energydiff << endl;
		double boltzmann = exp(-energydiff / T);
		double r = ourrandom();

		if (r < boltzmann)
		{
			//cout << "flip" << endl;
			energycurrent = energycurrent + energydiff;
			config[chosensite][0] = new_x;
			config[chosensite][1] = new_y;
			config[chosensite][2] = new_z;
			//accept += 1;
		}
		else
		{
		//	reject += 1;
			//cout<< "rejected"<< endl;
			double local_1_x = 0;
			double local_1_y = 0;
			double local_1_z = 0;
			for (int k = 0; k < neighbors[chosensite].size(); k++)
			{	
				local_1_x +=  -(J + djm[chosensite][k]) * config[neighbors[chosensite][k]][0];
				local_1_y +=  -(J + djm[chosensite][k]) * config[neighbors[chosensite][k]][1];
				local_1_z +=  -(J + djm[chosensite][k]) * config[neighbors[chosensite][k]][2];			
			}		
			

			double local_3_x = 0;
			double local_3_y = 0;
			double local_3_z = 0;
			for (int k = 0; k < third_neighbors[chosensite].size(); k++)
			{	
				local_3_x +=  config[third_neighbors[chosensite][k]][0];
				local_3_y +=  config[third_neighbors[chosensite][k]][1];
				local_3_z +=  config[third_neighbors[chosensite][k]][2];			
			}		
			local_3_x *= (-Jthird);
			local_3_y *= (-Jthird);
			local_3_z *= (-Jthird);

			vector<double> axis = {local_1_x + local_3_x, local_1_y + local_3_y, local_1_z + local_3_z};
			/*vector<double> local_field_3 = {0, 0, 0};
			for (int k = 0; k < third_neighbors[chosensite].size(); k++)
			{
				for (int c = 0; c < 3; c++) {local_field_3[c] = local_field_3[c] + config[third_neighbors[chosensite][k]][c];}	
				//std::transform (local_field_3.begin(), local_field_3.end(), config[third_neighbors[chosensite][k]].begin(), local_field_3.begin(), std::plus<double>());
			}		
			for (int i = 0; i < 3; i++) local_field_3[i] = local_field_3[i] * (-Jthird);
			*/			
			double beta = 2.0 * M_PI * ourrandom();
			double sx = config[chosensite][0];
			double sy = config[chosensite][1];
			double sz = config[chosensite][2];
			double Sx = axis[0];
			double Sy = axis[1];
			double Sz = axis[2];
			double S =sqrt(Sx*Sx+Sy*Sy+Sz*Sz);
			double cos_theta = Sz / S;
			double sin_theta = sqrt(1-cos_theta * cos_theta);
			double cos_phi = Sx / (S * sin_theta);
			double sin_phi = Sy / (S * sin_theta); 
			double cos_alpha = (sx*Sx + sy*Sy + sz*Sz)/S;
			double sin_alpha = sqrt(1-cos_alpha *cos_alpha);			
			double cos_beta = cos(beta);
			double sin_beta = sin(beta);

			double new_x = sin_alpha * cos_beta * cos_theta * cos_phi - sin_alpha * sin_beta *sin_phi + cos_alpha * sin_theta * cos_phi;
			double new_y = sin_alpha * cos_beta * cos_theta * sin_phi + sin_alpha * sin_beta *cos_phi + cos_alpha * sin_theta * sin_phi;
			double new_z = -sin_alpha * cos_beta* sin_theta + cos_alpha  * cos_theta;
				config[chosensite][0] = new_x;
				config[chosensite][1] = new_y;
				config[chosensite][2] = new_z;
			//cout<<new_x<<" "<<new_y <<" "<<new_z<<endl;
			//cout<<new_x*new_x + new_y*new_y + new_z*new_z<<endl;
			/*
			double phi_rot = -beta;
			vector<double> unit_vector = {Sx/S, Sy/S, Sz/S};
			vector<double> new_spin = rotate(config[chosensite], unit_vector, phi_rot);
			cout<<new_spin[0]<<" "<<new_spin[1]<<" "<<new_spin[2]<<endl;*/
		}
		//}

		//}
        if (n>=nwates && n % (nsites) == 0)
		{       
			energySum += energycurrent;
			energySquareSum += energycurrent*energycurrent;
			pointNumber++;
			Eoutfile<<energycurrent<<endl;
		//	cout<<energySum<<" "<<pointNumber<<endl;

		}
		if ( n>nwates && (n-nwates) % ((nsamples-nwates)/outNum) == 0)         //outfile after relaxation to get equilibrium spin map
        {	
	    cout<<djm[0][0]<<endl;
            int num = (n-nwates) / ((nsamples-nwates)/outNum);
            string name = outfilename+to_string(num-1) + ".txt";
            outfile[num].open(name);
            for (int i = 0; i<nsites; i++) { outfile[num] << i << "=" << config[i][0] << "=" << config[i][1] << "=" << config[i][2]  << endl; }
            outfile[num].close();
        }


		

	}
	if (outNum > 0) {
		int num = outNum - 1;
		name = outfilename+to_string(num) + ".txt";
		outfile[num].open(name);
		for (int i = 0; i<nsites; i++) { outfile[num] << i << "=" << config[i][0] << "=" << config[i][1] << "=" << config[i][2] << "="  << endl; }
		outfile[num].close();
	}
	

	double aveEsquare = energySquareSum / (double)pointNumber;
	double aveEnergy = energySum / (double)pointNumber;
	double specificHeat = (aveEsquare - aveEnergy*aveEnergy) / (double)(T*T*L*L);
	//Eoutfile<<aveEnergy<<" "<<specificHeat<<endl;
	return 0;
	//      cout <<"program finished"<<endl;
}
