/////////////////////////////////////////////////////////////
//   Main.cpp - Gaussian Elimination Program
//
//   Written by Jonathan Dunder on October 6, 2005
//
/////////////////////////////////////////////////////////////

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "matrix.h"

//Gaussian elimination program
//Written by Jonathan Dunder on October 6, 2005
matrix<double> GaussianElimination(matrix<double>A, std::vector<double>b, std::vector<double>x);

//Main function, calls gaussian elimination function
//Written by Jonathan Dunder on October 6, 2005
int main(){
	//Get a filename from the user
	std::string filename;
	std::cout << "Enter filename of matrix data\n";
	std::cin >> filename;
	std::cout << "\n\n";

	//Put the values from the file in a vector
	std::vector<double>values;
	unsigned int increment=1;
	unsigned int rows,columns;
	double temp;
	std::ifstream fin; //create input file stream
	fin.open(filename.c_str());
	if(fin) {
		std::cout << "Extracting values from configuration file...\n\n";
	} else {
		fin.close(); //Close file stream if file does not exist or is empty
		fin.clear();
		std::cout << "Fatal Error: Filename specified does not exist or is not accessible\n";
	}
	while(fin>>temp){
		if(increment==1){
			rows=temp;
		}else if(increment==2){
			columns=temp-1;
		}else {
			values.push_back(temp);
		}
		increment++;
	}
	
	//Put the values from the file into a matrix
	unsigned int currentrow=1;
	unsigned int currentcolumn=1;
	matrix<double>A(rows,columns); //A contains the coefficients for the equations
	std::vector<double>b; //b contains the answers to the equations
	std::vector<double>x; //X will contain the solution
	for(unsigned int i=0; i<values.size(); ++i){
		if(currentcolumn==columns+1){
			b.push_back(values[i]);
			currentcolumn=1;
			currentrow++;
		}else{
			A(currentrow,currentcolumn)=values[i];
			currentcolumn++;
		}
	}

	//Fill the solution vector with temporary values
	for(unsigned int i=1; i<=rows; ++i){
		x.push_back(1);
	}
	
	//Call the gaussian elimination program
	matrix<double>solution=GaussianElimination(A,b,x);

	return 0;
}

//Gaussian elimination program
//Written by Jonathan Dunder on October 6, 2005
matrix<double> GaussianElimination(matrix<double>A, std::vector<double>b, std::vector<double>x){
	//Start Gaussian Elimination
	int row=0;
	while(row<A.size()){
		//normalize
		double a=A(row+1,row+1);
		unsigned int col=0;
		while(col<A.size()){
			A(row+1,col+1)=A(row+1,col+1)/a;
			col++;
		}
		b[row]=b[row]/a;
		//end normalize
		unsigned int lowerRow=row+1;
		while(lowerRow<A.size()){
			//doLowerRow
			a=A(lowerRow+1,row+1);
			col=row;
			while(col<A.size()){
				A(lowerRow+1,col+1)=A(lowerRow+1,col+1)-a*A(row+1,col+1);
				col++;
			}
			b[lowerRow]=b[lowerRow]-a*b[row];
			//end doLowerRow
			lowerRow++;
		}
		row++;
	}
	//solve
	row=b.size()-1;
	while(row>=0){
		double sum=0.0;
		unsigned int col=row+1;
		while(col<b.size()){
			sum=sum+A(row+1,col+1)*x[col];
			col++;
		}
		x[row]=(b[row]-sum)/A(row+1,row+1);
		row--;
	}
	
	//Make a matrix to return
	matrix<double>returnmatrix(x.size(),1,0);

	//Fill the matrix with the solution to the system of equations
	for(unsigned int i=0; i<x.size(); ++i){
		std::cout << "x" << i+1 << "\=" << x[i] << std::endl;
		returnmatrix(i+1,1)=x[i];
	}
	
	return returnmatrix; //return the solution matrix
}