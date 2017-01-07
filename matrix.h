/////////////////////////////////////////////////////////////
//   Matrix.h - Matrix Class
//
//   Written by Jonathan Dunder on October 8, 2005
//
/////////////////////////////////////////////////////////////

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <assert.h>
#include <cmath>

//Matrix class
template <class T>
class matrix{
public:
	matrix(unsigned int rows=1, unsigned int columns=1, T fill=T() );
	matrix(const matrix<T>& original);
	~matrix();
	matrix<T> operator*(T value);
	void augment(const matrix<T>& left, const matrix<T>& right);
	void operator=(const matrix<T>& rhs);
	bool operator==(const matrix<T>& original);
	bool operator!=(const matrix<T>& original);
	matrix<T> operator+(const matrix<T>& original);
	T& operator()(unsigned int row,unsigned int column);
	matrix<T> operator()(unsigned int row);
	matrix<T> extractrow(unsigned int row);
	matrix<T> extractcolumn(unsigned int column);
	void pivot(unsigned int row1, unsigned int row2);
	void removerow(unsigned int row);
	void removecolumn(unsigned int column);
	T determinant(matrix<T>input);
	matrix<T> UpperTriangle(matrix<T>m);
	matrix<T> Adjoint(matrix<T>a);
	matrix<T> Inverse(matrix<T>a);
	matrix<T> Transpose(matrix<T>a);
	matrix<T> invert();
	unsigned int size();
	matrix<T> operator *(matrix<T>b);
	void insertion(std::ostream& os) const;
private:
	std::vector<T> data; //The vector stores all of the values
	int iDF; //A multiplying factor used for numerous functions
	unsigned int numberofcolumns; //The number of columns in the matrix
	unsigned int numberofrows; //The number of rows in the matrix
};

//<< Operator - calls the insertion function
//Written by Jonathan Dunder on October 8, 2005
template<class T>
std::ostream& operator<<(std::ostream& os, const matrix<T> mat){
	mat.insertion(os);
	return os;
}

//Size - Returns the number of rows
//Written by Jonathan Dunder on October 8, 2005
template<class T>
unsigned int matrix<T>::size(){
	return numberofrows;
}

//Constructor - Creates a matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T>::matrix(unsigned int rows, unsigned int columns, T fill):numberofrows(rows),numberofcolumns(columns){
	iDF=1;
	for(unsigned int i=0; i<(rows*columns); ++i){
		data.push_back(fill);
	}
}

//Copy Constructor - Creates a matrix from another
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T>::matrix(const matrix<T>& original){
	if(this!=&original){
		data=original.data;
		numberofcolumns=original.numberofcolumns;
		numberofrows=original.numberofrows;
	}
}

//Destructor - Destroys the matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T>::~matrix(){
	data.empty();
}

//() Operator - Returns a value from a matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
T& matrix<T>::operator()(unsigned int row,unsigned int column){
	if( numberofcolumns>=column && numberofrows>=column && data.size()>=((row*column)-1) ){
		if(row==1){
			return data[ (column)-1 ];
		}else{
			return data[ ( (row-1)*numberofcolumns+column )-1 ];
		}
	}else{
		std::cerr << "The specified location in the matrix does not exist\n\n";
		abort();
	}
}

//() Operator - Returns a row from a matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::operator()(unsigned int row){
	if( row<=numberofrows && row!=0 ){
		matrix<T>newmatrix(1,numberofcolumns,T());
		while(newmatrix.data.size()>0){
			newmatrix.data.pop_back();
		}
		if(row==1){
			for(unsigned int i=0; i<numberofcolumns; ++i){
				newmatrix.data.push_back(data[i]);
			}
		}else{
			for(unsigned int i=0; i<numberofcolumns; ++i){
				newmatrix.data.push_back(data[ (row-1)*numberofcolumns+i ]);
			}
		}
		return newmatrix;
	}else{
		std::cerr << "The specified row in the matrix does not exist\n\n";
		abort();
	}
}

//Extractrow - Returns a row from a matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::extractrow(unsigned int row){
	if( row<=numberofrows && row!=0 ){
		matrix<T>newmatrix(1,numberofcolumns,T());
		while(newmatrix.data.size()>0){
			newmatrix.data.pop_back();
		}
		if(row==1){
			for(unsigned int i=0; i<numberofcolumns; ++i){
				newmatrix.data.push_back(data[i]);
			}
		}else{
			for(unsigned int i=0; i<numberofcolumns; ++i){
				newmatrix.data.push_back(data[ (row-1)*numberofcolumns+i ]);
			}
		}
		return newmatrix;
	}else{
		std::cerr << "The specified row in the matrix does not exist\n\n";
		abort();
	}
}

//Extractcolumn - Returns a column from a matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::extractcolumn(unsigned int column){
	if( column<=numberofcolumns && column!=0 ){
		matrix<T>newmatrix(numberofrows,1,T());
		while(newmatrix.data.size()>0){
			newmatrix.data.pop_back();
		}
		if(column==1){
			for(unsigned int i=0; i<numberofrows; ++i){
				newmatrix.data.push_back( data[i*column] );
			}
		}else{
			for(unsigned int i=0; i<numberofrows; ++i){
				newmatrix.data.push_back( data[((numberofcolumns*i)-1)+column]);
			}
		}
		return newmatrix;
	}else{
		std::cerr << "The specified row in the matrix does not exist\n\n";
		abort();
	}
}

//== Operator - Checks if two matrices are equal
//Written by Jonathan Dunder on October 8, 2005
template<class T>
bool matrix<T>::operator==(const matrix<T>& original){
	bool returnval=0;
	if(original.data.size()==data.size() && original.numberofcolumns==numberofcolumns && original.numberofrows==numberofrows){
		returnval=1;
		for(unsigned int i=0; i<data.size(); ++i){
			if( data[i]!=original.data[i] ){
				returnval=0;
			}
		}
	}
	return returnval;
}

//= Operator - Sets one matrix equal to another
//Written by Jonathan Dunder on October 8, 2005
template<class T>
void matrix<T>::operator=(const matrix<T>& rhs){
	if(this!=&rhs){
		data=rhs.data;
		numberofcolumns=rhs.numberofcolumns;
		numberofrows=rhs.numberofrows;
	}
}

//!= Operator - Checks if two matrices are not equal
//Written by Jonathan Dunder on October 8, 2005
template<class T>
bool matrix<T>::operator!=(const matrix<T>& original){
	bool returnval=0;
	if(original.data.size()==data.size() && original.numberofcolumns==numberofcolumns && original.numberofrows==numberofrows){
		returnval=1;
		for(unsigned int i=0; i<data.size(); ++i){
			if( data[i]!=original.data[i] ){
				returnval=0;
			}
		}
	}
	return !returnval;
}

//* Operator - Multiplies a matrix by a scalar value
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::operator*(T value){
	matrix<T> returnmatrix=*this;
	for(unsigned int i=0; i<data.size(); ++i){
		returnmatrix.data[i]*=value;
	}
	return returnmatrix;
}

//Pivot - Swaps two rows
//Written by Jonathan Dunder on October 8, 2005
template<class T>
void matrix<T>::pivot(unsigned int row1, unsigned int row2){
	T temp;
	if(row1<=numberofrows && row2<=numberofrows && row1!=row2){
		for(unsigned int i=0; i<numberofcolumns; ++i){
			temp=data[ ((row1-1)*numberofcolumns)+i ];
			data[ ((row1-1)*numberofcolumns)+i ]=data[ ((row2-1)*numberofcolumns)+i ];
			data[ ((row2-1)*numberofcolumns)+i ]=temp;
		}
	}
}

//+ Operator - Allows addition of matrices
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::operator+(const matrix<T>& original){
	matrix<T>returnmatrix=*this;
	if( numberofrows==original.numberofrows && numberofcolumns==original.numberofcolumns ){
		for(unsigned int i=0; i<data.size(); ++i){
			returnmatrix.data[i]+=original.data[i];
		}
	}
	return returnmatrix;
}

//Augment - Combines two matrices
//Written by Jonathan Dunder on October 8, 2005
template<class T>
void matrix<T>::augment(const matrix<T>& left, const matrix<T>& right){
	if(left.numberofrows==right.numberofrows){
		numberofcolumns=left.numberofcolumns+right.numberofcolumns;
		numberofrows=left.numberofrows;
		while(data.size()>0){
			data.pop_back();
		}
		for(unsigned int i=0; i<numberofrows; ++i){
			for(unsigned int z=0; z<left.numberofcolumns; ++z){
				data.push_back(left.data[i*left.numberofcolumns+z]);
			}
			for(unsigned int z=0; z<right.numberofcolumns; ++z){
				data.push_back(right.data[i*right.numberofcolumns+z]);
			}
		}
	}
}

//Removerow - Removes a specified row from the matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
void matrix<T>::removerow(unsigned int row){
	if(row<=numberofrows){
		int currentrow=1;
		int currentcolumn=1;
		std::vector<T>data2;
		for(unsigned int i=0; i<data.size(); ++i){
			if(currentrow!=row){
				data2.push_back(data[i]);
			}
			currentcolumn++;
			if(currentcolumn>numberofcolumns){
				currentcolumn=1;
				currentrow++;
			}
		}
		data=data2;
		numberofrows--;
	}	
}

//Removecolumn - Removes a specified column from the matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
void matrix<T>::removecolumn(unsigned int column){
	if(column<=numberofcolumns){
		int currentcolumn=1;
		std::vector<T>data2;
		for(unsigned int i=0; i<data.size(); ++i){
			if(currentcolumn!=column){
				data2.push_back(data[i]);
			}
			currentcolumn++;
			if(currentcolumn>numberofcolumns)currentcolumn=1;
		}
		data=data2;
		numberofcolumns--;
	}
}

//Determinant - Finds the determinant of a matrix and returns it
//Written by Jonathan Dunder on October 8, 2005
template<class T>
T matrix<T>::determinant(matrix<T>input){
	int tms=input.numberofrows;
	T det=1;
	input=UpperTriangle(input);
	for(unsigned int i=1; i<=tms; i++){
		det=det*input(i,i);
	}
	det=det*iDF;
	return det;
}

//UpperTriangle - Returns a submatrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::UpperTriangle(matrix<T>m){
	T f1=0;
	T temp=0;
	int tms=m.numberofrows;
	int v=1;
	iDF=1;
	bool breaker=0;

	for(int col=1; col<=tms-1; col++){
		for(int row=col+1; row<=tms; row++){
			v=1;
			while(m(col,col)==0 && !breaker){
				if((col+v)>=tms){
					iDF=0;
					breaker=1;
				}else{
					for(int c=1; c<=tms; c++){
						temp=m(col,c);
						m(col,c)=m(col+v,c);
						m(col+v,c)=temp;
					}
					v++;
					iDF*=-1;
				}
			}

			if(m(col,col)!=0){
				f1=(-1)*m(row,col)/m(col,col);
				for(int i=col; i<=tms; i++){
					m(row,i)=f1*m(col,i)+m(row,i);
				}
			}
		}
	}
	return m;
}

//Adjoint - finds the adjoint of the matrix passed into it
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::Adjoint(matrix<T>a){
	int tms=a.numberofrows;
	matrix<T>m(tms,tms);
	int ii,jj,ia,ja;
	T det;
	for(unsigned int i=1; i<=tms; i++){
		for(unsigned int j=1; j<=tms; j++){
			ia=ja=1;
			matrix<T>ap(tms-1,tms-1);
			for(ii=1; ii<=tms; ii++){
				for(jj=1; jj<=tms; jj++){
					if((ii!=i)&&(jj!=j)){
						ap(ia,ja)=a(ii,jj);
						ja++;
					}
				}
				if((ii!=i)&&(jj!=j)){
					ia++;
				}
				ja=1;
			}
			det=determinant(ap);
			m(i,j)=pow(-1,i+j)*det;
		}
	}
	m=Transpose(m);
	return m;
}

//Inverse - finds the matrix
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::Inverse(matrix<T>a){
	int tms=a.numberofrows;
	matrix<T>m(tms,tms);
	matrix<T>mm=Adjoint(a);
	T det=determinant(a);
	T dd=0;
	dd=1/det;
	for(unsigned int i=1; i<=tms; i++){
		for(unsigned int j=1; j<=tms; j++){
			m(i,j)=dd*mm(i,j);
		}
	}
	return m;
}

//Transpose - Switches rows into columns
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::Transpose(matrix<T>a){
	matrix<T>m(a.numberofcolumns,a.numberofcolumns);
	for(unsigned int i=1; i<=a.numberofcolumns; i++){
		for(unsigned int j=1; j<=a.numberofcolumns; j++){
			m(j,i)=a(i,j);
		}
	}
	return m;
}

//invert - inverts the matrix passed into it
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::invert(){
	matrix<T>returnmatrix=Inverse(*this);
	return returnmatrix;
}

//* operator - supports cross-multiplication
//Written by Jonathan Dunder on October 8, 2005
template<class T>
matrix<T> matrix<T>::operator *(matrix<T>b){
	matrix<T>a=*this;
	assert(a.numberofcolumns==b.numberofrows);
	matrix<T>returnmatrix(a.numberofrows,b.numberofcolumns,0);
	for(unsigned int i=1; i<=a.numberofrows; ++i){
		for(unsigned int j=1; j<=b.numberofcolumns; ++j){
			T product=0;
			for(unsigned int z=1; z<=a.numberofcolumns; z++){
				product+=a(i,z)*b(z,j);
			}
			returnmatrix(i,j)=product;
		}
	}
	return returnmatrix;
}

//Insertion - Outputs the matrix values
//Written by Jonathan Dunder on October 8, 2005
template<class T>
void matrix<T>::insertion(std::ostream& os) const{
	matrix<T>temp=*this;
	for(unsigned int i=1; i<=numberofrows; ++i){
		for(unsigned int j=1; j<=numberofcolumns; ++j){
			std::cout << "Row " << i << ", Column " << j << "\=" << temp(i,j) << std::endl;
		}
	}
}

#endif