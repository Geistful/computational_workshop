#pragma GCC optimize("Ofast")
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<double> vd;
typedef vector<long long> vl;
typedef pair<int,int> pi;

void set_p_q_max(vector<vd> m, int k, int &p, int &q) { //p и q это позиция наибольшего элемента в матрице
	double maxvalue=0;
	for (int i=k; i<m.size(); i++) {
		for (int j=k; j<m[0].size()-1; j++) {
			if (abs(m[i][j])>abs(maxvalue)) {
				p=i, q=j;
				maxvalue = m[i][j];
			}
		}
	}
	return;
}

vector<vd> swap_row(vector<vd> &m, int p, int k) { //Меняем строки
	if (p==k) return m;
	double tmp;
	for (int i=0; i<m[0].size(); i++) {
		tmp = m[p][i];
		m[p][i] = m[k][i];
		m[k][i] = tmp;
	}
	return m;
}

vector<vd> swap_col(vector<vd> &m, int q, int k, vi &xorder) { //Меняем столбцы
	if (q==k) return m;
	double tmp;
	for (int i=0; i<m.size(); i++) {
		tmp = m[i][q];
		m[i][q] = m[i][k];
		m[i][k] = tmp;
	}
	tmp = xorder[q];
	xorder[q]=xorder[k];
	xorder[k]=tmp;
	return m;
}

void print_long_matr(vector<vd> m) { //Печать матрицы вида A|B (надо напечатать только A)
	for (int i=0; i<m.size(); i++) {
		for (int j = 0; j < m[0].size(); j++) {
			if (j==m[0].size()-1)
				cout << "| " << m[i][j] << "\n";
			else cout << m[i][j] << " ";
		}
	}
}

void print_matr(vector<vd> m) {		//Просто печать матрицы
	for (int i=0; i<m.size(); i++) {
		for (int j = 0; j < m.size(); j++) {
			if (j==m.size()-1)
				cout << m[i][j] << "\n";
			else cout << m[i][j] << " ";
		}
	}
}

void print_vect(vd v) {				//Печать вектора
	for (int i=0; i<v.size(); i++) {
		if (i==0) 
			cout << "(" << v[i] << ", ";
		else if (i==v.size()-1)
			cout << v[i] << ")\n";
		else
			cout << v[i] << ", ";
	}
}

vi gauss_forw_elim (vector<vd> &newmatr) {	//"Прямой ход" метода Гаусса
	vi xorder;
 	for (int i=0; i<newmatr.size(); i++) {
  		xorder.push_back(i);
  	}
  	int p, q;
  	set_p_q_max(newmatr, 0, p, q);
  	double temp, maxval;
  	for (int k=0; k<newmatr.size(); k++) {
	  	maxval=0;
	  	swap_col(newmatr, q, k, xorder);
	  	swap_row(newmatr, p, k);
	  	temp=newmatr[k][k];
	  	for (int j=k; j<newmatr[0].size(); j++) {
	  		newmatr[k][j]/=temp;
	  	}
	  	for (int i=k+1; i<newmatr.size(); i++) {
	  		temp = newmatr[i][k];
	  		for (int j=k; j<newmatr[0].size(); j++) {
	  			newmatr[i][j]-=newmatr[k][j]*temp;
	  			if (abs(newmatr[i][j])>abs(maxval) && j!=newmatr[0].size()-1) {
	  				p=i, q=j;
	  				maxval = newmatr[i][j];
	  			}
	  		}
	  	}
	}
	return xorder;
}

vd gauss_back_subs (vector<vd> newmatr, vi xorder) { //Обратный ход+сортировка корней
	vd x(newmatr.size(), 0);
  	x[newmatr.size()-1]=newmatr[newmatr.size()-1][newmatr[0].size()-1];
  	for (int i=newmatr.size()-2; i>=0; i--) {
  		x[i]=newmatr[i][newmatr[0].size()-1];
  		for (int j=i+1; j<newmatr[0].size()-1; j++) {
  			x[i]-=newmatr[i][j]*x[j];
  		}	
  	}
  	vd x_ordered(newmatr.size(), 0);
  	for (int i=0; i<newmatr.size(); i++) {
  		x_ordered[i]=x[xorder[i]];
  	}
  	return x_ordered;
}

vd gauss_method (vector<vd> &newmatr) {	 	//Собственно метод Гаусса - прямой + обратный ход	
	vi xorder = gauss_forw_elim(newmatr);
	vd x = gauss_back_subs(newmatr, xorder);
	return x;
}

vd jordan_method (vector<vd> &newmatr) {	//Метод Жордана - обратный ход оставляет только диагональ
	vi xorder = gauss_forw_elim(newmatr);
	vd coef;
	for (int k=1; k<newmatr.size(); k++) {
		for (int i=k; i<newmatr[0].size(); i++) {
			for (int j=0; j<k; j++) {
				if (i==k) {
					if (j==k-1) coef.push_back(newmatr[j][i]);
					else coef[j]=newmatr[j][i];
				}
				newmatr[j][i]-=newmatr[k][i]*coef[j];
			}
		}
	}
	vd x(newmatr.size(), 0);		//Здесь тоже сортируем корни
  	for (int i=0; i<newmatr.size(); i++) {
  		x[i]=newmatr[xorder[i]][newmatr[0].size()-1];
  	}
	return x;
}

void sol_abs_diff (vector<vd> m, vd x) { //Функция для модуля невязки Ax и B
  for (int i=0; i<m.size(); i++) {
  	double sum=0;
  	for (int j=0; j<m[0].size()-1; j++) {
  		sum+=m[i][j]*x[j];
  	}
  	cout << "Sum of " << i+1 << " row = " << sum 
  		 << " | Abs diff from expected = " << abs(sum-m[i][m[0].size()-1]) << "\n";
  }
}

vector<vd> mult_matr (vector<vd> m1, vector<vd> m2) { //Произведение двух матриц
  if (m1[0].size()!=m2.size())
  	throw std::invalid_argument("Dimensions doesn't match up");
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

vector<vd> rev_matrix (vector<vd> matr) {		//Получение обратной матрицы (подаем на вход расширенную!)
	vector<vd> rev (matr.size(), vd (matr.size(), 0));
	for (int i=0; i<matr.size(); i++) {
		vector<vd> temp = matr;
		for (int j=0; j<temp.size(); j++) {
			if (j==i) temp[j][temp[0].size()-1]=1;
			else temp[j][temp[0].size()-1]=0;
		}
		vd x = jordan_method(temp);		//Используем метод Жодана для i=0..n
		for (int k=0; k<matr.size(); k++) {
			rev[k][i] = x[k];
		}
	}
	return rev;
}

class lu_matrices {		//Исходя из смысловой связи, создал класс, объединяющий матрицы L и U
public:
	vector<vd> l;
	vector<vd> u;
};

lu_matrices lu_decomp (vector<vd> &matr) { //Разложение матрицы на L и U
	int cnt = matr.size();
	lu_matrices lu;
	lu.l.resize(cnt, vd (cnt, 0));
	lu.u.resize(cnt, vd (cnt, 0));
	for (int i=0; i<cnt; i++) {
		lu.u[i][i]=1;
		for (int j=i; j<cnt; j++) {
			lu.l[j][i]=matr[j][i];
			for (int k=0; k<i; k++) {
				lu.l[j][i]-=lu.l[j][k]*lu.u[k][i];
			}
			lu.u[i][j]=matr[i][j];
			for (int k=0; k<i; k++) {
				lu.u[i][j]-=lu.l[i][k]*lu.u[k][j];
			}
			lu.u[i][j]/=lu.l[i][i];
		}
	}
	return lu;
}

double one_norm (vector<vd> m) { //Подсчет нормы ||.||_1 для матрицы
	double maxvalue=0;
	for (int i=0; i<m.size(); i++) {
		double colsum=0;
		for (int j=0; j<m.size(); j++) {
			colsum+=abs(m[j][i]);
		}
		maxvalue=max(maxvalue, colsum);
	}
	return maxvalue;
}

double vect_norm (vd v) {	//Норма для вектора
	double res = 0;
	for (int i=0; i<v.size(); i++) {
		res += v[i]*v[i];
	}
	return sqrt(res);
}

vd num_vec_prod (double k, vd v) { //Число умножается на вектор
	for (int i=0; i<v.size(); i++) {
		v[i]*=k;
	}
	return v;
} 

vd vect_sum (vd a, vd b) {		//Сумма векторов
	for (int i=0; i<a.size(); i++) {
		a[i]+=b[i];
	}
	return a;
}

double get_cond (vector<vd> m) {	//Число обусловленности
	vector<vd> rev = rev_matrix(m);
	return one_norm(rev)*one_norm(m);
}

vector<vd> unite_a_b (vector<vd> a, vd b) {		//Соединяем матрицы (A+B=A|B)
	vector<vd> newa(a.size(), vd(a.size()+1, 0));
	for (int i=0; i<newa.size(); i++) {
		for (int j=0; j<newa[0].size(); j++) {
			if (j==newa[0].size()-1)
				newa[i][j]=b[i];
			else newa[i][j]=a[i][j];
		}
	}
	return newa;
}

int main()
{
  //ios::sync_with_stdio(0);
  //cin.tie(0);
  vector<vd> matr {
  	{8.29381, 0.995516, -0.560617, 0.766522},
	{0.995516, 6.298198, 0.595772, 3.844422},
	{-0.560617, 0.595772, 4.997407, 5.239231}
  };

  cout << "Initial matrix:\n-------------\n";
  print_long_matr(matr);

  vector<vd> newmatr(matr.size(), vd (matr[0].size(), 0));
  copy(begin(matr), end(matr), begin(newmatr));
  vector<vd> jord_matr(matr.size(), vd (matr[0].size(), 0));
  copy(begin(matr), end(matr), begin(jord_matr));

  cout << "--------------\nJordan method\n";
  cout << "Resulting matrix:\n-------------\n";
  vd jord_roots = jordan_method(jord_matr);
  print_long_matr(jord_matr);

  cout << "--------------\nGauss method\n";
  cout << "Forward elimination results:\n-------------\n";
  vd gau_roots = gauss_method(newmatr);
  print_long_matr(newmatr);

  cout << "--------------\nRoots in order after backward substitution:\n-------------\n";
  for (int i=0; i<newmatr.size(); i++) {
  	cout << "x[" << i+1 << "] = " << gau_roots[i] <<"\n";
  }

  cout << "\nChecking roots:\n";
  sol_abs_diff(matr, gau_roots);

  cout << "\n-----------\nRev Matrix\n------------\n";
  vector<vd> rev = rev_matrix(matr);
  print_matr(rev);

  vector<vd> small_matr (matr.size(), vd(matr.size(), 0));
  for (int i = 0; i < matr.size(); i++) {
  	for (int j = 0; j < matr.size(); j++) {
  		small_matr[i][j]=matr[i][j];
  	}
  }
  cout << "-----------\nResult of multiplication of Rev. matr and initial matrix:\n";
  print_matr(mult_matr(small_matr, rev));

  lu_matrices dec = lu_decomp(matr);
  cout << "\nL matrix:\n";
  print_matr(dec.l);
  cout << "\nU matrix:\n";
  print_matr(dec.u);
  cout << "\nL*U:\n";
  vector<vd> luprod = mult_matr(dec.l, dec.u);
  print_matr(luprod);
  cout << "\nSum of abs. diff. between corresponding elements of A and LU:\n";
  double sum=0;
  for (int i=0; i<matr.size(); i++) {
  	for (int j = 0; j< matr.size(); j++) {
  		sum+=abs(luprod[i][j]-matr[i][j]);
  	}
  }
  cout << sum << "\n";

  cout << "\nSince A=LU => det(A)=det(L)=";
  double diag_l_prod =1;
  for (int i = 0; i < dec.l.size(); ++i) {
  	diag_l_prod*=dec.l[i][i];
  }
  cout << diag_l_prod << "\n";

  //------------------------------------------
  vector<vd> a {
  	{-401.43, 200.19},
	{1201.14, -601.62}
  };
  vd b {200.1, -600.1};
  vd delta_b {199.1, -601.1};
  
  vector<vd> ab = unite_a_b(a, b);
  vector<vd> rev_a = rev_matrix(a);
  double cond_a = get_cond(ab);
  vd x_gau = gauss_method(ab);
  vector<vd> a_delta_b = unite_a_b(a, delta_b);
  vd delta_x_gau = gauss_method(a_delta_b);
  cout << "\n---------------------------\n";
  cout << "Condition number task:\n";
  cout << "Matrix A:\n";
  print_matr(a);
  cout << "\nb = ";
  print_vect(b);
  cout << "delta b = ";
  print_vect(delta_b);
  cout << "\nSol. for b:\nx = ";
  print_vect(x_gau);
  cout << "Sol. for delta_b:\ndelta x = ";
  print_vect(delta_x_gau);
  cout << "\ncond(A) = " << cond_a << "\n";
  cout << "Factual difference for solutions = " 
  	   << vect_norm(vect_sum(delta_x_gau, 
  	   				num_vec_prod(-1.0, x_gau)))/
  	   	  vect_norm(x_gau) << "\n";
  cout << "Theoretical difference = " 
       << cond_a*(vect_norm(vect_sum(delta_b, num_vec_prod(-1.0, b)))/vect_norm(b)) 
       << "\n";
  return 0;
}