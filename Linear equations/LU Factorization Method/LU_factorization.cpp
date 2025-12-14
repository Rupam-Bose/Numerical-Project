#include<bits/stdc++.h>
#include<fstream>
using namespace std;


//-----A=LU : Decomposing the Conefficient Matrix into the product of Lower Triangular and Upper Triangular Matrices-----
void decompose(vector<vector<double>>&a, vector<vector<double>>&u, vector<vector<double>>&l, vector<double>&b, int n)
{
    u=a;

    //-----Finding the Upper Triangular Matrix using Gauss Elimination with Partial Pivoting-----
    for(int i=0; i<n-1; i++) { 
        int pivot = i;
        for(int j=i+1; j<n; j++) { 
            if(fabs(u[j][i]) > fabs(u[pivot][i])) pivot = j;
        } 

        if(pivot!=i) {
            swap(u[i], u[pivot]);
            swap(a[i], a[pivot]);
            swap(b[i], b[pivot]);
        }

        if(fabs(u[i][i]) < 1e-9) continue;

        double x = u[i][i];
        for(int k=i+1; k<n; k++) {
            double y=u[k][i];  

            for(int j=i; j<n; j++) { 
                u[k][j] -= u[i][j] * y/x;
            }
        } 
    }

    //------Computing Lower Triangular Matrix-----
    for(int i=0; i<n; i++) {
        l[i][i]=1; 

        for(int j=i+1; j<n; j++) {
            double sum = 0;

            for(int k=0; k<i; k++) {
                sum += l[j][k] * u[k][i];
            }

            if(fabs(u[i][i]) < 1e-9) l[j][i] = 0;  
            else l[j][i] = (a[j][i] - sum)/u[i][i];  
        
        }
    }
}


//-----LUx = b -> Ly = b : Forward Substitution Phase for y which is the solution of the Augmented Matrix [L|b]----
void forward(vector<vector<double>>&l, vector<double>&y, vector<double>&b, int n) 
{
    for(int i=0; i<n; i++) {
        double sum=0;
        for(int j=0; j<i; j++) {
            sum += y[j]*l[i][j]; 
        }
        y[i] = (b[i]-sum) / l[i][i]; 
    }
}


//-----Ux = y : Backward Substitution Phase for x which is the solution of the Augmented Matrix [U|y]-----
void backward(vector<vector<double>>&u, vector<double>&x, vector<double>&y, int n)
{
    for(int i=n-1; i>=0; i--) {
        double sum=0;
        for(int j=i+1; j<n; j++) {
            sum += x[j]*u[i][j];
        }
        x[i] = (y[i]-sum) / u[i][i];
    }
}


//------Ax = b : Printing the Linear System-------
void systemPrint(vector<vector<double>>&a, vector<double> &b, ofstream &fout)
{
    fout<<"The Linear System: \n";
    int n=a.size();
    for(int i=0; i<n; i++) {
        fout<<a[i][0]<<"*x"<<1<<" ";
        for(int j=1; j<n; j++) {
            if(a[i][j]>=0) fout<<" + ";
            fout<<a[i][j]<<"*x"<<j+1<<" ";
        }
        fout<<" = "<<b[i]<<'\n';
    }
}



void solve(ifstream &fin, ofstream &fout)
{
    //------Reading the number of linear equations of the linear system------
    int n; 
    fin>>n; 

    vector<vector<double>> a(n, vector<double>(n)); // a is the coefficient matrix
    vector<vector<double>> l(n, vector<double>(n, 0)); // l is the lower triangular matrix
    vector<vector<double>> u(n, vector<double>(n, 0)); // u is the upper triangular matrix
    vector<double>b(n); // b is the constants matrix 
    

    //-----Reading the Augmented Matrix into the Coefficient and Constant Matrices------
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            fin>>a[i][j];
            if(i==j) l[i][j]=1;
        }
        fin>>b[i];
    }

    systemPrint(a, b, fout);

    decompose(a,u,l,b,n);

    vector<double>y(n, 0); // y=Ux;
    forward(l,y,b,n);
    
    //-----Handling the solution types of the system-----
    if(u[n-1][n-1]==0)
    {
        if(y[n-1]) {
            fout<<"\nThe system has NO SOLUTION.\n";
        }
        else {
            fout<<"\nThe system has INFINITE SOLUTIONS.\n";
        }
    } 

    else {
        fout<<"\nThe system has UNIQUE SOLUTION. \n";

        vector<double>x(n, 0); 
        backward(u,x,y,n);

        //------Printing the solutions of the System------
        fout<<"\nSolutions: \n";
        for(int i=0; i<n; i++) {
            fout<<"x"<<i+1<<" = "<<x[i]<<'\n';
        }
    } 
}


int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout<<"Solution of Linear System using LU Decomposition Method\nMultiple Testcases\n\n";

    if(!fin) {
        fout<<"File not found.\n";
        return 0;
    }

    //------Mutiple Test Cases------
    int t;
    fin>>t;

    int cse=1;
    
    while(t--)
    {
        for(int i=0;i<45;i++) fout<<"=";
        fout<<'\n';
        fout<<"Test case: "<<cse<<"\n";
        solve(fin, fout);
        cse++; 
    }
}