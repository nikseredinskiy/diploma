#include <iostream>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;
typedef Eigen::Triplet<double> Trip;

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

MatrixXd r(MatrixXd x, double w, int N){
	return fi(x, w, N) - x;
}


//Seems like it works now
MatrixXd createMatrixV(Array<VectorXd, Dynamic, 1> r_i, int N, int M){
	MatrixXd temp(N*N, M);
	for (int i = 0; i < M; i++){
		temp.col(i) = r_i(i) - r_i(M);
	}
	return temp;
}
/*
SparseMatrix <double> createSparseMatrixV(Array<MatrixXd, Dynamic, 1> r_i, int N, int M) {
	SparseMatrix<double> temp(N, M);
	vector<Trip> tripletList;
	//N * M = estimation_of_entries
	tripletList.reserve(N * M);

	for (int i = 0; i < M; i++){
		for (int j = 0; j < M; j++){
			tripletList.push_back(Trip(i, j, r_i(i)(j)));
		}
	}

	temp.setFromTriplets(tripletList.begin(), tripletList.end());

	return temp;
}*/

MatrixXd lsdamp(Array<MatrixXd, Dynamic, 1> x_i, int M, double w, int N){
	Array<VectorXd, Dynamic, 1> r_i(M + 1);
	for (int i = 0; i < M + 1; i++){
		MatrixXd tempRI = r(x_i(i), w, N);
		Map<VectorXd> r_i_vector(tempRI.data(), tempRI.size());
		r_i(i) = r_i_vector;
	}
	MatrixXd V(N*N, M);
	V = createMatrixV(r_i, N, M);

	MatrixXd r_temp = -r(x_i(M), w, N);
	Map<VectorXd> r_temp_vector(r_temp.data(), r_temp.size());

	VectorXd c = V.colPivHouseholderQr().solve(r_temp_vector);

	MatrixXd newU(N, N);
	newU.setZero();

	for (int i = 0; i < M; i++){
		newU = newU + c(i)*x_i(i);
	}

	newU = newU + (1 - c.sum())*x_i(M);

	cout << "Matrix after lsdamp:" << endl;
	cout << newU << endl;

	return newU;
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
	//TODO: replace by setZero() method
	//DONE
	u.setZero();
	
	//Initialization for Fs-1, Fs-2 and Fs
	//TODO:replace by setZero() method
	//DONE
	MatrixXd f_s_1(N, N), f_s_2(N, N), f_s(N, N);
	f_s_1.setZero();
	f_s_2.setZero();
	f_s.setZero();

	
	Array<MatrixXd, Dynamic, 1> x_i(M + 1);
	x_i(0) = u;
	int countOfTheProcessRun = 0;
	while (r(u, w, N).norm() > E){
		countOfTheProcessRun++;
		//is there need to put M + 1 instead of M
		//due to x_i will not have x_i(M) item
		for (int k = 1; k < M + 1; k++){
			int numberOfIterations = 0;
			for (int _s = 0; _s < s; _s++){
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
			u = f_s; //x_i that will be passed to lsdamp function
			cout << "M = " << k << " | iterNum = " << s << endl;
			cout << "Matrix u:" << endl;
			cout << u << endl;

			x_i(k) = f_s;
			//
		}

		//now we have x_0...x_m and can send it to lsdamp
		u = lsdamp(x_i, M, w, N);
	}
	
	cout << "Process has been run " << countOfTheProcessRun << " times" << endl;

	cout << endl;
	system("pause");
}