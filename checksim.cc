#include<iostream>
#include<iomanip>
#include<fstream>

using namespace std;

void checkvel(){

	int i1, i2, i3, j1, j2, j3, Nd, id;
	double a1, a2, a3, b1, b2, b3;
	

	ifstream in1, in2;
	
	in1.open("serial_data.txt");
	in2.open("parallel_data.txt");
	
	while( !(in1.eof()) ){
		in1>>i1>>i2>>i3>>a1>>a2>>a3;
		in2>>j1>>j2>>j3>>b1>>b2>>b3;
		
		if(i1 == j1 && i2 == j2 && i3 == j3) {
	//		if( abs( a1-b1 ) > 1e-10 || abs( a2-b2 ) > 1e-10 || abs( a3-b3 ) > 1e-10  ) {
				cout<<setw(10)<<i1<<setw(10)<<i2<<setw(10)<<i3<<setw(20)<<abs( a1-b1 )<<setw(20)<<abs( a2-b2 )<<setw(20)<<abs( a3-b3 )<<endl;
	//		}
		}
		else{
			cout<<"DIMENTIONS UNMATCHED";
		}		
	}
}

void checkspectra(){

	int ishell_1, ishell_2;	
	double E1, W1, E2, W2;
	

	ifstream in1, in2;
	
	in1.open("serial_data.txt");
	in2.open("parallel_data.txt");
	
	while( !(in2.eof()) ){
		in1>>ishell_1>>E1>>W1;
		in2>>ishell_2>>E2>>W2;
		
		if(ishell_1 == ishell_2) {
//			if( abs( (E1-E2)/E1 ) > 1e-10 || abs( (W2-W2)/W1 ) > 1e-10  ) {
				cout<<setw(10)<<ishell_2<<setw(20)<<abs( (E1-E2))<<setw(20)<<abs( (W1-W2) )<<endl;
//			}
		}
//		else{
//			cout<<"DIMENTIONS UNMATCHED";
//		}		
	}
	
	in1.close();
	in2.close();
}

void check_EW(){

	double t1, t2, E1, W1, E2, W2;
	

	ifstream in1, in2;
	
	in1.open("serial_data.txt");
	in2.open("parallel_data.txt");
	
	while( !(in2.eof()) ){
		in1>>t1>>E1>>W1;
		in2>>t2>>E2>>W2;
		if( abs(t1-t2) > 1e-10 || abs(E1 - E2) > 1e-10 || abs(W1 - W2) > 1e-10   ){
			cout<<setw(10)<<abs(t1-t2)<<setw(20)<<abs( (E1-E2))<<setw(20)<<abs( (W1-W2) )<<endl;
		}
	}
	
	in1.close();
	in2.close();
}

int main(){

	checkvel();
//	checkspectra();
//	check_EW();
	return 0;
}

















