#include<bits/stdc++.h>
#include<fstream>
using namespace std;


//-----Forward Elimination phase: Converting the Conefficient Matrix into Upper Triangular Matrix with pivoting------
void forwardd(vector<vector<double>>&cof)
{
    int n=cof.size();
    for(int i=0; i<n-1; i++) {
        int pivot = i;
        for(int j=i+1; j<n; j++) { 
            if(fabs(cof[j][i]) > fabs(cof[pivot][i])) pivot = j;
        } 
        
        if(pivot!=i) swap(cof[i], cof[pivot]);

        if(fabs(cof[i][i]) < 1e-9) continue;

        double x = cof[i][i];

        for(int k=i+1; k<n; k++) {
            double y=cof[k][i]; 
            for(int j=i; j<n+1; j++) { 
                cof[k][j]= cof[k][j] - cof[i][j]*y/x;
            }
        } 
    }
}


//------Back Substition Phase: Finding the solution vector------
void backwardd(vector<vector<double>>&cof, vector<double>&ans)
{
    int n=cof.size();
    for(int i=n-1; i>=0; i--) { 
        double sum=0;
        for(int j=i+1; j<n; j++) { 
            sum += ans[j]*cof[i][j];
        }  
        ans[i] = (cof[i][n]-sum) / cof[i][i]; 
    }
}


//------Printing the Linear System-------
void systemPrint(vector<vector<double>>&cof, ofstream &fout)
{
    fout<<"The Linear System: \n";
    int n=cof.size();
    for(int i=0; i<n; i++) {
        fout<<cof[i][0]<<"*x"<<1<<" ";
        for(int j=1; j<n; j++) {
            if(cof[i][j]>=0) fout<<" + ";
            fout<<cof[i][j]<<"*x"<<j+1<<" ";
        }
        fout<<" = "<<cof[i][n]<<'\n';
    }
}


void solve(ifstream &fin, ofstream &fout)
{
    //------Reading the number of linear equations of the linear system------
    int n; 
    fin>>n; 

    vector<vector<double>>cof(n, vector<double>(n+1));
    
    //------Reading the Conefficient Matrix------
    for(int i=0; i<n; i++) {
        for(int j=0; j<n+1; j++) {
            fin>>cof[i][j];
        }
    }

    systemPrint(cof, fout);

    forwardd(cof);

    //------Handling the solutions type of the system------
    if(cof[n-1][n-1]==0 && cof[n-1][n]) {
        fout<<"\nThe system has NO SOLUTION.\n";
    }

    else if(cof[n-1][n-1]==0 && cof[n-1][n]==0) {
        fout<<"\nThe system has INFINITE SOLUTIONS.\n";
    }

    else {
        fout<<"\nThe system has UNIQUE SOLUTION.\n";
        vector<double>ans(n, 0);
        backwardd(cof, ans);

        //------Printing the solutions of the System------
        fout<<"\nSolutions: \n";
        for(int i=0; i<n; i++) {
            fout<<"x"<<i+1<<" = "<<ans[i]<<'\n';
        }
    }
}


int main()
{
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    fout<<"Solution of Linear System using Gauss Elimination Method\nMultiple Testcases\n\n";

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