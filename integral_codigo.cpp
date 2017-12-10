/** Thiago Henrique Neves Coelho
  * Integracao Dupla Numerica
  */

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

/* Soma todos os numeros da lista e coloca o resultado em v[0] */
/* Opera a soma de forma alternada para manter os termos em mesma ordem de grandeza e reduzir o erro numerico */
void soma(vector<double> &v) {
	int A = v.size();
	while(A != 1) {
		for(int i = 0; i < A/2; i++) v[i] = v[2*i] + v[2*i + 1];
		if(A % 2 != 0) v[0] += v[A-1];
		A /= 2;
	}
}

/* Metodo de Simpson */
/* Calcula a integral dupla utilizando interpolacoes com parabolas */
double simpson(double (*f)(double,double), double ax, double bx, double (*infX)(double), double (*supX)(double), int N, int M) {
	/* Garante que o numero de particoes seja par */
	if(N % 2 == 1) N++;
	if(M % 2 == 1) M++;
	
	double hx = (bx - ax)/N; // Tamanho do passo no eixo x
	
	/* Calcula a integral em y para valores fixos de x */
	vector<double> integ, integFinal;
	for(int i = 0; i < N + 1; i++) {
		double hy = (supX(ax + i*hx) - infX(ax + i*hx))/M; // Tamanho do passo no eixo y
		double ay = (infX(ax + i*hx)); // Posicao inicial do y
		
		/* Fixa um valor de x e integra em y */
		/* Valor de x fixado eh ax + i*hx */
		/* hy*(f0 + 4*f1 + 2*f2 + 4*f3 + ... + 4*fn-1 + fn)/3 */
		for(int j = 0; j < M + 1; j++) {
			if(j == 0 || j == M) integ.push_back(f(ax + i*hx, ay + j*hy)*hy/3.0);
			else {
				if(j % 2 == 1) integ.push_back(4.0*f(ax + i*hx, ay + j*hy)*hy/3.0);
				else integ.push_back(2.0*f(ax + i*hx, ay + j*hy)*hy/3.0);
			}
		}
		
		/* Soma todos as partes das integrais das interpolacoes para encontrar o valor da integral em y com um x fixo e guarda esse valor */
		soma(integ);
		integFinal.push_back(integ[0]);
		integ.clear();
	}

	/* Calcula a integral em x */
	/* hx*(f0 + 4*f1 + 2*f2 + 4*f3 + ... + 4*fn-1 + fn)/3 */
	for(int i = 0; i < N + 1; i++) {
		if(i == 0 || i == N) integ.push_back(integFinal[i]*hx/3.0);
		else {
			if(i % 2 == 1) integ.push_back(4.0*integFinal[i]*hx/3.0);
			else integ.push_back(2.0*integFinal[i]*hx/3.0);
		}
	}
	soma(integ); // Soma integrais das interpolacoes
	
	double ans = integ[0];
	integ.clear();
	integFinal.clear();
	
	return ans;
}

/* Metodo do retangulo */
/* Calcula a integral somando os valores da funcao multiplicados ao tamanho do intervalo */
double retangulo(double (*f)(double,double), double ax, double bx, double (*infX)(double), double (*supX)(double), int N, int M) {
	double hx = (bx - ax)/N; // Tamanho do passo no eixo x

	/* Calcula a integral em y para valores fixos de x */
	vector<double> integ, integFinal;
	for(int i = 0; i < N; i++) {
		double hy = (supX(ax + i*hx) - infX(ax + i*hx))/M; // Tamanho do passo no eixo y
		double ay = (infX(ax + i*hx)); // Valor inicial do y
		
		/* Fixa o x em ax + i*hx e calcula a integral em y */
		/* hy*(f0 + f1 + f2 + ... + fn) */
		for(int j = 0; j < M; j++) {
			integ.push_back(f(ax + i*hx, ay + j*hy)*hy);
		}
		
		/* Soma as partes da integral e guarda o valor total */
		soma(integ);
		integFinal.push_back(integ[0]);
		integ.clear();
	}
	
	/* Calcula a integral em x */
	/* hx*(f0 + f1 + f2 + ... + fn) */
	for(int i = 0; i < N; i++) {
		integ.push_back(integFinal[i]*hx);
	}
	soma(integ); // Soma integrais das interpolacoes
	
	double ans = integ[0];
	integ.clear();
	integFinal.clear();
	return ans;
}

/* Metodo do trapezio */
/* Calcula a integral dupla utilizando interpolacoes com retas */
double trapezio(double (*f)(double,double), double ax, double bx, double (*infX)(double), double (*supX)(double), int N, int M) {
	double hx = (bx - ax)/N; // Tamanho do passo no eixo x
	
	/* Calcula a integral em y para valores fixos de x */
	vector<double> integ, integFinal;
	for(int i = 0; i < N + 1; i++) {
		double hy = (supX(ax + i*hx) - infX(ax + i*hx))/M; // Tamanho do passo no eixo y
		double ay = (infX(ax + i*hx)); // Valor inicial de y
		
		/* Calcula a integral em y para x = ax + i*hx */
		/* hy*(f0 + 2f1 + 2f2 + 2f3 + ... + 2fn-1 + fn)/2 */
		for(int j = 0; j < M + 1; j++) {
			if(j == 0 || j == M) integ.push_back(f(ax+i*hx, ay+j*hy)*hy/2.0);
			else integ.push_back(f(ax+i*hx, ay+j*hy)*hy);
		}
		
		/* Soma integrais das interpolacoes e guarda o resultado */
		soma(integ);
		integFinal.push_back(integ[0]);
		integ.clear();
	}
	
	/* Integra em x */
	/* hx*(f0 + 2f1 + 2f2 + 2f3 + ... + 2fn-1 + fn)/2 */
	for(int i = 0; i < N + 1; i++) {
		if(i == 0 || i == N) integ.push_back(integFinal[i]*hx/2.0);
		else integ.push_back(integFinal[i]*hx);
	}
	soma(integ); // Soma integrais das interpolacoes
	
	double ans = integ[0];
	integ.clear();
	integFinal.clear();
	
	return ans;
}

/* Limite inferior de y */
double inferiorFx(double x) {
	return 0.0;
}

/* Limite superior de y */
double superiorFx(double x) {
	return x*x;
}

/* Exemplo de funcao de duas variaveis para ser integrada */
double func_val(double x, double y) {
	return x*cos(y);
}

int main(void) {
	double ax = 0.0; // Limite inferior de integracao de x
	double bx = 7.0; // Limite superior de integracao de x
	
	int N = 200; // Numero de segmentos no eixo x
	int M = 200; // Numero de segmentos no eixo y
	
	printf("Metodo do retangulo: %.4lf\n", retangulo(func_val, ax, bx, inferiorFx, superiorFx, N, M));
	printf("Metodo do trapezio: %.4lf\n", trapezio(func_val, ax, bx, inferiorFx, superiorFx, N, M));
	printf("Metodo de Simpson: %.4lf\n", simpson(func_val, ax, bx, inferiorFx, superiorFx, N, M));
	
	return 0;
}
