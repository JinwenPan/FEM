#include "set"
#include "iostream"
#include "vector"
#include "fstream"
#include "sstream"
#include "cmath"
#include "ColsammSource/Colsamm.h"
#include "iterator"

struct Vertex{

	void add_neighbor(int new_id){
		if (neighbors.empty()){
			neighbors.push_back(id);
			neighbors.push_back(new_id);
		}
		else {
			for (size_t i=0; i<neighbors.size(); ++i){
				if (neighbors[i] == new_id){
					break;
				}
				else if (i==neighbors.size()-1){
					neighbors.push_back(new_id);
				}
			}
		}
	}
	int id;
	double x, y;
	std::vector<int> neighbors;
};

struct Triangle{
	Triangle(): local_M(3,std::vector<double>(3)), local_A(3,std::vector<double>(3)){} 
	Vertex* a;
	Vertex* b;
	Vertex* c;

	std::vector< std::vector< double > > local_M;
	std::vector< std::vector< double > > local_A;
};

double delta = 0.0;

double k_sq(double x, double y){
	return (100+delta)*exp(-50*(x*x + y*y)) -100;
}

int main(int argc, char* argv[]){
	delta = std::stod(argv[1]);
	double epsilon = std::stod(argv[2]);
	std::string meshfile_name = "unit_circle.txt";

//========================================= build vertices and faces vectors =========================
std::ifstream meshfile(meshfile_name);
std::string line;
std::string word;
std::getline(meshfile, line);	// Get first line
std::istringstream vertexline(line); // vertex num is here
vertexline>>word;		// First word in meshfile is number of vertices
int vertices_num=std::stoi(word);
std::getline(meshfile, line); // skip the next line

// ------------- vertex vector -------------------
std::vector<Vertex> vertices(vertices_num);

for (int i=0; i<vertices_num; ++i){
	std::getline(meshfile, line);
	std::istringstream inputline (line);

	inputline >> word;
	vertices[i].id = std::stoi(word);
	inputline >> word;
	vertices[i].x = std::stod(word);
	inputline >> word;
	vertices[i].y = std::stod(word);
}
//-------------------------------------------------------
std::getline(meshfile, line);
std::istringstream faceline(line);	// num of faces is here
faceline>>word;
int faces_num=std::stoi(word);
// std::cout<<faces_num<< std::endl;
std::string word0;
std::string word1;
std::string word2;
int id0, id1, id2;
std::getline(meshfile, line); // skip the next line
// -------------------- faces vector ---------------
std::vector<Triangle> faces(faces_num);
for (int i=0; i<faces_num; ++i){
	std::getline(meshfile, line);
	std::istringstream inputline (line);
	inputline>>word0;	inputline>>word1;	inputline>>word2;
	id0 = std::stoi(word0);	id1 = std::stoi(word1); id2 = std::stoi(word2); 

	faces[i].a=&vertices[id0];
	faces[i].b=&vertices[id1];
	faces[i].c=&vertices[id2];

}
// meshfile.close();
// ---------- add neighbors for each vertex -----------------
for (int i=0; i<faces_num; ++i){
		faces[i].a->add_neighbor(faces[i].b->id);
		faces[i].a->add_neighbor(faces[i].c->id);
		faces[i].b->add_neighbor(faces[i].a->id);
		faces[i].b->add_neighbor(faces[i].c->id);
		faces[i].c->add_neighbor(faces[i].a->id);
		faces[i].c->add_neighbor(faces[i].b->id);
}

// -----------------------------------------------some examples
std::cout<< "x value of vertex with id 7 = " << vertices[7].x<< std::endl;
std::cout<< "id of third vertex in face 0 = " <<faces[0].a->id<< std::endl;
std::cout<< "number of neighbors for the third vertex in face 1975 = " <<faces[1975].c->neighbors.size()<< std::endl;

std::cout<< "neighbors of second vertex in face 2: " << std::endl;
for (size_t i=0; i<faces[2].b->neighbors.size(); ++i){
	std::cout<<faces[2].b->neighbors[i] << std::endl;
}
std::cout<<"==============================" <<std::endl;

//--------------- find k2 ---------------------------------
std::vector<double> k2(vertices_num);
for (int i=0; i<vertices_num; ++i){
	k2[i]=k_sq(vertices[i].x, vertices[i].y);
}

// std::ofstream k2_output;
// k2_output.open("k2_file.txt");
// for (int i=0; i<vertices_num; ++i){
// 	k2_output << vertices[i].x << "\t" << vertices[i].y << "\t" << k2[i] << std::endl;
// }
// k2_output.close();
//================================================== /////// M matrix ////////// ============================

//-------------------------------- build local M --------------------------
for (int i=0; i<faces_num; ++i){

	std::vector<double> local_vertices(6,0.0);

	local_vertices[0] = faces[i].a->x;
	local_vertices[1] = faces[i].a->y;
	local_vertices[2] = faces[i].b->x;
	local_vertices[3] = faces[i].b->y;
	local_vertices[4] = faces[i].c->x;
	local_vertices[5] = faces[i].c->y;
	
	using namespace ::_COLSAMM_;
	ELEMENTS::Triangle ref_element;
	ref_element(local_vertices);
	faces[i].local_M = ref_element.integrate(v_() * w_());
}
//----------------------- build global M ---------------------
std::vector<std::vector<double>> global_M(vertices_num,std::vector<double>(vertices_num));
for (int i=0; i<faces_num; ++i){

	int a_id=faces[i].a->id;
	int b_id=faces[i].b->id;
	int c_id=faces[i].c->id;
	std::vector< std::vector< double > > local_M = faces[i].local_M;

	global_M[a_id][a_id] = global_M[a_id][a_id] + local_M[0][0];
	global_M[a_id][b_id] = global_M[a_id][b_id] + local_M[0][1];
	global_M[a_id][c_id] = global_M[a_id][c_id] + local_M[0][2];
	
	global_M[b_id][a_id] = global_M[b_id][a_id] + local_M[1][0];
	global_M[b_id][b_id] = global_M[b_id][b_id] + local_M[1][1];
	global_M[b_id][c_id] = global_M[b_id][c_id] + local_M[1][2];

	global_M[c_id][a_id] = global_M[c_id][a_id] + local_M[2][0];
	global_M[c_id][b_id] = global_M[c_id][b_id] + local_M[2][1];
	global_M[c_id][c_id] = global_M[c_id][c_id] + local_M[2][2];

}

//======================================================//////// A matrix ///////////=========================================
// ------------------------ build local A  -------------------------------

for (int i=0; i<faces_num; ++i){
	std::vector<double> local_vertices(6,0.0);

	local_vertices[0] = faces[i].a->x;
	local_vertices[1] = faces[i].a->y;
	local_vertices[2] = faces[i].b->x;
	local_vertices[3] = faces[i].b->y;
	local_vertices[4] = faces[i].c->x;
	local_vertices[5] = faces[i].c->y;

	using namespace ::_COLSAMM_;
	ELEMENTS::Triangle ref_element;
	ref_element(local_vertices);
	faces[i].local_A = ref_element.integrate(grad(v_()) * grad(w_()) - func<double>(k_sq) * v_() * w_());

}

// ----------------------- build global A ---------------------
std::vector<std::vector<double>> global_A(vertices_num,std::vector<double>(vertices_num));
for (int i=0; i<faces_num; ++i){

	int a_id=faces[i].a->id;
	int b_id=faces[i].b->id;
	int c_id=faces[i].c->id;
	std::vector< std::vector< double > > local_A = faces[i].local_A;

	global_A[a_id][a_id] = global_A[a_id][a_id] + local_A[0][0];
	global_A[a_id][b_id] = global_A[a_id][b_id] + local_A[0][1];
	global_A[a_id][c_id] = global_A[a_id][c_id] + local_A[0][2];
	
	global_A[b_id][a_id] = global_A[b_id][a_id] + local_A[1][0];
	global_A[b_id][b_id] = global_A[b_id][b_id] + local_A[1][1];
	global_A[b_id][c_id] = global_A[b_id][c_id] + local_A[1][2];

	global_A[c_id][a_id] = global_A[c_id][a_id] + local_A[2][0];
	global_A[c_id][b_id] = global_A[c_id][b_id] + local_A[2][1];
	global_A[c_id][c_id] = global_A[c_id][c_id] + local_A[2][2];

}

// for (int i=0; i<5 ; ++i){
// 	for (int j=0 ; j<vertices_num; ++j){
// 		std::cout<< global_A[i][j]<< " ";
// 	}
// 	std::cout<< " ////////////////////////////////////////// "<<std::endl;
// }
// 	std::cout<<faces[0].local_A[0][0]<<" "<<faces[0].local_A[0][1]<<" "<<faces[0].local_A[0][2]<<" "<<std::endl;
// 	std::cout<<faces[0].local_A[1][0]<<" "<<faces[0].local_A[1][1]<<" "<<faces[0].local_A[1][2]<<" "<<std::endl;
// 	std::cout<<faces[0].local_A[2][0]<<" "<<faces[0].local_A[2][1]<<" "<<faces[0].local_A[2][2]<<" "<<std::endl;

// std::cout<<global_A[118][905]<< "hhhhhhhhhhhhhhhhhhhhhh";
//=================================== Inverse power iteration ========================================

double lambda_new, lambda_old, lambda_res=1.0;
std::vector<double> u_h(vertices_num,0.008);
std::vector<double> u_h_old(vertices_num);
std::vector<double> f_h(vertices_num);

int ipi_count=0;
int ipi_iters = 5;
std::cout << "will do " << ipi_iters << " number of ipi iters " << std::endl;
while (lambda_res>1e-10 || lambda_res<-1e-10){
	lambda_old=lambda_new;

	if (ipi_count > ipi_iters){break;}
	 std::cout << "ipi = "<< ipi_count++ << " lambda_res = " << lambda_res << std::endl;

	for (int i=0; i<vertices_num; ++i){		// f = M  * u
		for (int j=0; j<vertices_num; ++j){
			f_h[i] = f_h[i] + global_M[i][j]*u_h[j];
		}
	}

	double GS_res=1.0;
	// int GS_count=0;

	// ---------------------------- GS solver ----------------- 
	while (GS_res>epsilon){	
		for (int i=0; i<vertices_num; ++i){

			u_h_old[i] = u_h[i];
		}
		
		for (int i=0; i<vertices_num; ++i){
			double sum=0.0;
			for (int j=0; j<vertices_num; ++j){
				if (i==j){continue;}
				sum = sum + global_A[i][j]*u_h[j]; 
			}
			u_h[i] = (1.0/global_A[i][i])*(f_h[i]-sum);
		}

		for (int i=0; i<vertices_num; ++i){
			u_h_old[i] = u_h[i] - u_h_old[i]; // for residual
		}	

		GS_res=0.0;
		for (int i=0; i<vertices_num; ++i){

			GS_res = GS_res + u_h_old[i]*u_h_old[i];
		}
		GS_res= sqrt(GS_res);
		// std::cout << "GS iter = "<< GS_count++ << "\t |  GS res = " << GS_res << std::endl;

	} // end GS solver

	//-------------- normalise u_h --------------
	double u_h_norm=0.0;
	for (int i=0; i<vertices_num; ++i){
		u_h_norm = u_h_norm + u_h[i]*u_h[i];
	}
	u_h_norm= sqrt(u_h_norm);
	for (int i=0; i<vertices_num; ++i){
		u_h[i] = u_h[i]/u_h_norm;
	}

	//---------------- find lambda -------------------------------
	double lambda_numerator=0.0;
	double lambda_denominator=0.0;
	for (int i=0; i<vertices_num; ++i){

		double sum_num=0.0;
		double sum_denom=0.0;
		for (int j=0; j<vertices_num; ++j){
			sum_num = sum_num + global_A[i][j] * u_h[j];
			sum_denom = sum_denom + global_M[i][j] * u_h[j];

		}

		lambda_numerator = lambda_numerator + sum_num * u_h[i];
		lambda_denominator = lambda_denominator + sum_denom * u_h[i];
		// std::cout << i <<" " << lambda_numerator << std::endl;
	}
	lambda_new = lambda_numerator / lambda_denominator;
	lambda_res=(lambda_new - lambda_old)/lambda_old;

}// end ipi

// for (int i=0; i<vertices_num; ++i){
// 	std::cout<< i << " " << u_h[i] << std::endl;
// }

std::ofstream u_h_output;
u_h_output.open("eigenmode.txt");
for (int i=0; i<vertices_num; ++i){
	u_h_output << vertices[i].x << "\t" << vertices[i].y << "\t" << u_h[i] << std::endl;
}
u_h_output.close();


}