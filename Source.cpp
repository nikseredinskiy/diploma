#include <iostream>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

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

double f_ij(int N, int i, int j, MatrixXd u){
	double result = N * N * (u_ij(i - 1, j, N, u) + u_ij(i + 1, j, N, u) - 4 * u_ij(i, j, N, u) + u_ij(i, j - 1, N, u) + u_ij(i, j + 1, N, u));
	double sum_p_q = 0;
	for (int p = 1; p < N; p++){
		for (int q = 1; q < N; q++){
			sum_p_q = sum_p_q + cosh(u(p, q));
		}
	}
	result = result - 10 * (sum_p_q * sum_p_q) / (N*N*N*N);
	return (double)result;
}

MatrixXd f(MatrixXd x, int N){
	MatrixXd temp(N, N);
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			//TODO: remove filling with zeros on the boundary of the matrix
			//TEST: filling with zeros has been removed
			temp(i, j) = f_ij(N, i, j, x);
		}
	}
	return temp;
}

MatrixXd fi(MatrixXd x, double w, int N){
	return x + w * f(x, N);
}

MatrixXd F_s_with_prev(MatrixXd x, int s, MatrixXd f_s_1, MatrixXd f_s_2 ,double w, int N){
	if (s == 0){
		return x;
	}
	if (s == 1){
		//DONE: changed alpha_j/beta_j/gamma_j parameter from s-1 to s(due to a page 587,Fadeev)
		return alpha_j(s-1)*fi(f_s_1, w, N) + beta_j(s-1) * f_s_1;
	}
	if (s > 1){
		return alpha_j(s-1)*fi(f_s_1, w, N) + beta_j(s-1) * f_s_1 - gamma_j(s-1) * f_s_2;
	}
}

void main(){
	int N = 5;
	int M = 14;
	int s = 101;
	double E = 0.00000001;
	double w = 1/((double)8*N*N);
	
	//cout << "Enter N = ";
	//cin >> N;
	MatrixXd u(N, N);

	//Initialization for u
	//TODO: remove certain conditions from u
	//DONE: removed certain conditions from u
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			u(i, j) = 0;
		}
	}

	cout << u << endl;
	
	//Initialization for Fs-1, Fs-2 and Fs
	MatrixXd f_s_1(N, N), f_s_2(N, N), f_s(N, N);
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			f_s_1(i, j) = 0;
			f_s_2(i, j) = 0;
			f_s(i, j) = 0;
		}
	}

	for (int k = 1; k < M; k++){
		for (int _s = 0; _s <= s; _s++){
			if (_s == 0){
				f_s = F_s_with_prev(u, s, f_s_1, f_s_2, w, N);
			}
			if (_s == 1){
				f_s_1 = f_s;
				f_s = F_s_with_prev(u, s, f_s_1, f_s_2, w, N);
			}
			if (_s > 2){
				f_s_2 = f_s_1;
				f_s_1 = f_s;
				f_s = F_s_with_prev(u, s, f_s_1, f_s_2, w, N);
			}
		}
		u = f_s;

		cout << "Martix u:" << endl;
		cout << u << endl;
	} 

	cout << endl;
	system("pause");
}