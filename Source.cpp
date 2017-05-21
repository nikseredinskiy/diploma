#include <iostream>
#include <fstream>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

double alpha_j(int j){
	if (j == 0){
		return (double)3 / 4.;
	}
	else {
		return (double)j * (2 * j + 1) / ((j + 1) * (j + 1));
	}
}

double beta_j(int j){
	if (j == 0){
		return (double)1 / 4.;
	}
	else {
		return (double)j / ((2 * j - 1) * (j + 1) * (j + 1));
	}
}

double gamma_j(int j){
	if (j == 0){
		return (double)0;
	}
	else {
		return (double)((2 * j + 1) * (j - 1) * (j - 1)) / ((2 * j - 1) * (j + 1) * (j + 1));
	}
}

//DONE: Created a function that will return a proper boundary value if needed
double u_ij(int i, int j, int N, MatrixXd u){
	if (i >= 0 && j >= 0 && i < N && j < N){
		return u(i, j);
	}
	if (j == -1){
		return 1 - (i + 1) / (double)N;
	}
	else {
		if (i == -1){
			return 1 - (j + 1)/(double)N;
		}
	}
	if (j == N){
		return 0;
	}
	else {
		if (i == N){
			return 0;
		}
	}
}

double f_ij(int N, int i, int j, MatrixXd u, double sum_p_q){
	double result = N * N * (u_ij(i - 1, j, N, u) + u_ij(i + 1, j, N, u) - 4 * u_ij(i, j, N, u) + u_ij(i, j - 1, N, u) + u_ij(i, j + 1, N, u));
	result = result - 10 * (sum_p_q * sum_p_q) / (N*N*N*N);
	return (double)result;
}

void f(MatrixXd x, int N, MatrixXd &result_f){
	double sum_p_q = 0;
	for (int p = 1; p < N; p++){
		for (int q = 1; q < N; q++){
			sum_p_q = sum_p_q + cosh(x(p, q));
		}
	}
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			//TODO: remove filling with zeros on the boundary of the matrix
			//TEST: filling with zeros has been removed
			result_f(i, j) = f_ij(N, i, j, x, sum_p_q);
		}
	}
}

void fi(MatrixXd x, double w, int N, MatrixXd &result_f, MatrixXd &result_fi){
	f(x, N, result_f);
	result_fi = x + w * result_f;
}

void F_s_with_prev(MatrixXd x, int s, MatrixXd f_s_1, MatrixXd f_s_2, double w, int N, MatrixXd &f_s, MatrixXd &result_f, MatrixXd &result_fi){
	if (s == 0){
		f_s = x;
	}
	if (s == 1){
		//DONE: changed alpha_j/beta_j/gamma_j parameter from s-1 to s(due to a page 587,Fadeev)
		fi(f_s_1, w, N, result_f, result_fi);
		f_s = alpha_j(s - 1)* result_fi + beta_j(s - 1) * f_s_1;
	}
	if (s > 1){
		fi(f_s_1, w, N, result_f, result_fi);
		f_s = alpha_j(s - 1)*result_fi + beta_j(s - 1) * f_s_1 - gamma_j(s - 1) * f_s_2;
	}
}

MatrixXd r(MatrixXd x, double w, int N, MatrixXd &result_f,MatrixXd &result_fi){
	fi(x, w, N, result_f, result_fi);
	return result_fi - x;
}


//Seems like it works now
void createMatrixV(Array<VectorXd, Dynamic, 1> &r_i, int N, int M, MatrixXd &V){
	for (int i = 0; i < M; i++){
		V.col(i) = r_i(i) - r_i(M);
	}
}

void countC(VectorXd &c, MatrixXd &V, Map<VectorXd> &r_temp_vector){
	c = V.colPivHouseholderQr().solve(r_temp_vector);
}

void lsdamp(Array<MatrixXd, Dynamic, 1> &x_i, int M, double w, int N, MatrixXd &u, MatrixXd &result_f, MatrixXd &result_fi, MatrixXd &V, VectorXd &c, Array<VectorXd, Dynamic, 1> &r_i){
	MatrixXd tempRI;
	for (int i = 0; i < M + 1; i++){
		tempRI = r(x_i(i), w, N, result_f,result_fi);
		Map<VectorXd> r_i_vector(tempRI.data(), tempRI.size());
		r_i(i) = r_i_vector;
	}
	createMatrixV(r_i, N, M, V);
	
	MatrixXd r_temp = -r(x_i(M), w, N, result_f,result_fi);
	Map<VectorXd> r_temp_vector(r_temp.data(), r_temp.size());
	
	countC(c, V, r_temp_vector);
	
	u.setZero();
	for (int i = 0; i < M; i++){
		u = u + c(i)*x_i(i);
	}
	u = u + (1 - c.sum())*x_i(M);

	VectorXd(M + 1).swap(c);
	MatrixXd(N*N, M).swap(V);

}

void main(){
	int N = 25;
	int M = 4;
	int s = 30;
	double E = 0.000000000001;
	double w = 1/((double)8*N*N);
	
	//cout << "Enter N = ";
	//cin >> N;
	MatrixXd u(N, N);

	//Initialization for u
	//TODO: remove certain conditions from u
	//DONE: removed certain conditions from u
	//TODO: replace by setZero() method
	//DONE
	u.setZero();
	
	//Initialization for Fs-1, Fs-2 and Fs
	//TODO:replace by setZero() method
	//DONE
	MatrixXd f_s_1(N, N), f_s_2(N, N), f_s(N, N), result_f(N, N), result_fi(N, N), V(N*N, M);
	VectorXd c(M + 1);
	Array<VectorXd, Dynamic, 1> r_i(M + 1);
	f_s_1.setZero();
	f_s_2.setZero();
	f_s.setZero();
	result_f.setZero();
	result_fi.setZero();
	
	Array<MatrixXd, Dynamic, 1> x_i(M + 1);
	int countOfTheProcessRun = 0;
	int iterationsCount = 0;
	ofstream resultsFile;
	resultsFile.open("results_25_4_30_with.csv");

	while (r(u, w, N, result_f,result_fi).norm() > E){
		x_i(0) = u;
		countOfTheProcessRun++;
		//is there need to put M + 1 instead of M
		//due to x_i will not have x_i(M) item
		for (int k = 1; k < M + 1; k++){
			for (int _s = 0; _s < s; _s++){
				iterationsCount++;
				if (_s == 0){
					F_s_with_prev(u, _s, f_s_1, f_s_2, w, N, f_s, result_f, result_fi);
				}
				if (_s == 1){
					f_s_1 = f_s;
					F_s_with_prev(u, _s, f_s_1, f_s_2, w, N, f_s, result_f, result_fi);
				}
				if (_s > 2){
					f_s_2 = f_s_1;
					f_s_1 = f_s;
					F_s_with_prev(u, _s, f_s_1, f_s_2, w, N, f_s, result_f, result_fi);
				}
				resultsFile << r(f_s, w, N, result_f, result_fi).norm() << ";" << iterationsCount << endl;
				cout << iterationsCount << endl;
			}
			u = f_s; //x_i that will be passed to lsdamp function
			//cout << "M = " << k << " | iterNum = " << s << endl;
			//cout << "Matrix u:" << endl;
			//cout << u << endl;

			x_i(k) = f_s;
			MatrixXd(N,N).swap(f_s);
			MatrixXd(N, N).swap(f_s_1);
			MatrixXd(N, N).swap(f_s_2);
		}
		//now we have x_0...x_m and can send it to lsdamp
		lsdamp(x_i, M, w, N, u, result_f, result_fi,V, c, r_i);
	}
	iterationsCount++;
	resultsFile << r(u, w, N, result_f, result_fi).norm() << ";" << iterationsCount << endl;
	cout << "Process has been run " << countOfTheProcessRun << " times" << endl;
	resultsFile.close();


	cout << endl;
	system("pause");
}