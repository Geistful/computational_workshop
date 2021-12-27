#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<double> vd;
typedef vector<long long> vl;
typedef pair<int,int> pi;

vector<vd> E_matr(int deg) {
	vector<vd> e(deg, vd (deg, 0));
	for (int i=0; i<deg; i++) {
		e[i][i]=1;
	}
	return e;
}

int sign (double temp) {
	if (temp>0) return 1;
	if (temp<0) return -1;
	return 0;
}

double set_p_q_max(vector<vd> m, int &p, int &q) {
	double maxvalue=0;
	for (int i=0; i<m.size(); i++) {
		for (int j=i+1; j<m.size(); j++) {
			if (abs(m[i][j])>maxvalue) {
				p=i, q=j;
				maxvalue = abs(m[i][j]);
			}
		}
	}
	return maxvalue;
}

vector<vd> matr_x_matr (vector<vd> m1, vector<vd> m2) {
  vector<vd> mult(m1.size(), vd (m2[0].size(), 0));
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m2[0].size(); i++) {
		for (int j = 0; j < m1[0].size(); j++) {
			mult[k][i] += m1[k][j]*m2[j][i];
		}
	}
  }
  return mult;
} 

vd matr_x_vect (vector<vd> m1, vd m2) {
  vd mult(m1.size(), 0);
  for (int k=0; k<m1.size(); k++) {
	for (int i=0; i<m1[0].size(); i++) {
		mult[k] += m1[k][i]*m2[i];
	}
  }
  return mult;
}

double get_d (vector<vd> m, int p, int q) {
	return sqrt((m[p][p]-m[q][q])*(m[p][p]-m[q][q])+4*m[p][q]*m[p][q]);
}

double get_c (vector<vd> m, int p, int q) {
	return sqrt(0.5*(1+(abs(m[p][p]-m[q][q])/get_d(m, p, q))));
}

double get_s (vector<vd> m, int p, int q) {
	return sign(m[p][q]*(m[p][p]-m[q][q]))*
		   sqrt(0.5*(1-(abs(m[p][p]-m[q][q])/get_d(m, p, q))));
}

vector<vd> a_get_next (vector<vd> a, int p, int q, double c, double s) {
	vector<vd> newa(a.size(), vd (a[0].size(), 0));
  	copy(begin(a), end(a), begin(newa));
	for (int i=0; i<a.size(); i++) {
		for (int j=0; j<a.size(); j++) {
			if (i!=p && i!=q && j!=p && j!=q)
				newa[i][j]=a[i][j];
			else if (i!=p && i!=q && j==p) {
				newa[i][j]=c*a[i][p]+s*a[i][q];   //change sign before s
				newa[j][i]=newa[i][j];
			}
			else if (i!=p && i!=q && j==q) {
				newa[i][j]=-s*a[i][p]+c*a[i][q];  //change sign before s
				newa[j][i]=newa[i][j];
			}
			else if (i==p && j==p)
				newa[i][j]=c*c*a[p][p]+2*c*s*a[p][q]+s*s*a[q][q];
			else if (i==q && j==q)
				newa[i][j]=s*s*a[p][p]-2*c*s*a[p][q]+c*c*a[q][q];
			else if (i==p && j==q) {
				newa[i][j]=0/*(c*c-s*s)*a[p][q]+c*s*(a[q][q]-a[p][p])*/;
				newa[j][i]=newa[i][j];
			}
		}
	}
	return newa;
}

void print_matr(vector<vd> m) {
	for (int i=0; i<m.size(); i++) {
		for (int j = 0; j < m.size(); j++) {
			if (j==m.size()-1)
				cout << m[i][j] << "\n";
			else cout << m[i][j] << " ";
		}
	}
}

vector<vd> x_get_next (vector<vd> x, int p, int q, double c, double s) {
	vector<vd> newx(x.size(), vd (x[0].size(), 0));
  	copy(begin(x), end(x), begin(newx));
	for (int i=0; i<newx.size(); i++) { //col
		for (int j=0; j<newx.size(); j++) { //row
			if (i!=p && i!=q)
				newx[j][i]=x[j][i];
			else if (i==p)
				newx[j][i]=c*x[j][p]+s*x[j][q];
			else if (i==q)
				newx[j][i]=-s*x[j][p]+c*x[j][q];
		}
	}
	return newx;
}

vector<vd> transpose_matr (vector<vd> m) {
	for (int i = 0; i < m.size(); ++i) {
		for (int j=i+1; j<m.size(); j++) {
			double temp = m[i][j];
			m[i][j] = m[j][i];
			m[j][i] = temp;
		}	
	}
	return m;
}

void print_vect(vd v) {
	for (int i=0; i<v.size(); i++) {
		if (i==0) 
			cout << "(" << v[i] << ", ";
		else if (i==v.size()-1)
			cout << v[i] << ")\n";
		else
			cout << v[i] << ", ";
	}
}

double vect_twonorm (vd v) {
	double sum=0;
	for (int i=0; i<v.size(); i++) {
		sum += v[i]*v[i];
	}
	sum = sqrt(sum);
	return sum;
}

vd vect_plus_vect (vd v1, vd v2) {
	vd v3(v1.size(), 0);
	for (int i=0; i<v1.size(); i++) {
		v3[i]=v1[i]+v2[i];
	}
	return v3;
}

vd double_x_vect (double num, vd v) {
	vd res(v.size(), 0);
	for (int i=0; i<v.size(); i++) {
		res[i]=num*v[i];
	}
	return res;
}

double vect_metric (vd v1, vd v2) {
	return vect_twonorm(vect_plus_vect(v1, double_x_vect(-1.0, v2)));
}

double apost_eval (vector<vd> a, double lambda, vd y) {
	return vect_twonorm(vect_plus_vect(matr_x_vect(a, y), double_x_vect(-1*lambda, y)))/vect_twonorm(y);
}

vd normalize_vect (vd v) {
	if (vect_twonorm(v)==1) return v;
	return double_x_vect(1.0/vect_twonorm(v), v);
}

double vect_absmax (vd v) {
	double maxval=0;
	for (int i=0; i<v.size(); i++) {
		maxval = max(maxval, abs(v[i]));
	}
	return maxval;
}

vd max_val_normalize (vd v) {
	double del = vect_absmax(v);
	for (int i = 0; i < v.size(); ++i) {
		v[i] /= del;	
	}
	return v;
}

double scalar_vect (vd v1, vd v2) {
	double res=0;
	for (int i=0; i<v1.size(); i++) {
		res += v1[i]*v2[i];
	}
	return res;
}

int main()
{
	vector<vd> matr {
  	{-1.48213, -0.03916, 1.08254},
	{-0.03916, 1.13958, 0.01617},
	{1.08254, 0.01617, -1.48271}
  	};
  	double eps = 0.001, c, s, bigeival=0;
  	int p, q;
  	vector<vd> x = E_matr(matr.size());
  	vector<vd> newmatr(matr.size(), vd (matr[0].size(), 0));
  	copy(begin(matr), end(matr), begin(newmatr));
  	while (set_p_q_max(newmatr, p, q)>=eps) {
  		c=get_c(newmatr, p, q);
		s=get_s(newmatr, p, q);
  		newmatr = a_get_next(newmatr, p, q, c, s);
  		x = x_get_next(x, p, q, c, s);
  	}
  	for (int i=0; i<newmatr.size(); i++) {
  		if (abs(newmatr[i][i])>bigeival) bigeival = newmatr[i][i];
  		cout << "lambda " << i+1 << " = " << newmatr[i][i] << "\n";
  	}
  	cout << "Eigenvectors:\n";
  	vector<vd> row_x = transpose_matr(x);
  	for (int i=0; i<x.size(); i++) {
  		cout << "E.v. " << i+1 << ": "; 
  		print_vect(row_x[i]);
  		cout << "|vect("<< i+1 <<")|_2= " << vect_twonorm(row_x[i]) <<"\n";
  		cout << " --> ";
  		if (vect_twonorm(row_x[i])==1) cout << "Already normalized\n";
  		else {
	  		print_vect(normalize_vect(row_x[i]));
  		}
  	}
  	vd y (matr.size(), -1.0);
  	vd prev = y;
  	double lamb1 = 0;
  	int iter = 0;
  	int fact_iter = 0;
  	while (apost_eval(matr, lamb1, y)>=eps) {
  		iter++;
  		prev = y;
  		y = matr_x_vect(matr, y);
  		lamb1 = y[0]/prev[0];
  		if (abs(lamb1-bigeival)>eps) fact_iter++;
  	}
  	cout << "Approximation eval: " <<apost_eval(matr, lamb1, y) << "\n";
  	cout << "Iterations needed = " << iter << "\n";
  	cout << "Factual Iterations = " << fact_iter << "\n";
  	cout << "Approx first lambda = " << lamb1 << "\n";
  	cout << "Approx first vector = ";
  	print_vect(y);
  	cout << "After normalization: ";
  	vd y_norm = normalize_vect(y);
  	print_vect(y_norm);
  	cout << "Checker:\n";
  	cout << "A*y = ";
  	print_vect(matr_x_vect(matr, y_norm));
  	cout << "lamb1*y = ";
  	print_vect(double_x_vect(lamb1, y_norm));
  	//SCALAR METHOD//
  	cout << "\nScalar Method:\n";
  	vd scal_appr(matr.size(), -1.0);
  	lamb1=0;
  	iter=0;
  	fact_iter=0;
  	while (apost_eval(matr, lamb1, scal_appr)>=eps*eps) {
  		iter++;
  		prev = scal_appr;
  		scal_appr = matr_x_vect(matr, scal_appr);
  		lamb1 = scalar_vect(scal_appr, prev)/scalar_vect(prev, prev);
  		if (abs(lamb1-bigeival)>eps*eps) fact_iter++;
  	}
  	cout << "Approximation eval: " << apost_eval(matr, lamb1, scal_appr) << "\n";
  	cout << "Iterations needed = " << iter << "\n";
  	cout << "Factual Iterations = " << fact_iter << "\n";
  	cout << "Approx first lambda = " << lamb1 << "\n";
  	cout << "Approx first vector = ";
  	print_vect(scal_appr);
  	cout << "After normalization: ";
  	vd scalar_appr_norm = normalize_vect(scal_appr);
  	print_vect(normalize_vect(scal_appr));
  	cout << "Checker:\n";
  	cout << "A*y = ";
  	print_vect(matr_x_vect(matr, y_norm));
  	cout << "lamb1*y = ";
  	print_vect(double_x_vect(lamb1, y_norm));
	return 0;
}
