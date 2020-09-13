#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include<sstream>
#include<string>
#include <math.h>

using namespace std;

/**
 * The base Matrix class for the required basic matrix operations.
 */
class Matrix
{
public:
    Matrix(/* args */){};

    void printMatrix(vector<vector<int>> &AA){
        for(auto i = AA.begin(); i < AA.end();i++){
            for(auto j = i->begin(); j < i->end(); j++){
                cout<<*j<<" ";
            }
            cout<<endl;
        }
    };
};

/**
 * The base Quries class for storing input queries.
 */
class Quries{
public:
    /* Query vector */
    vector<string> query;

    Quries(char *filename){
        freopen(filename, "r", stdin);
        string str;
        while(getline(cin, str)){
            query.push_back(str);
        }
        cin.clear();
        fclose(stdin);
    };
};

/**
 * The base Matrix class for storing input attributes.
 */
class Relations{
public:
    /* attribute vector */
    vector<string> attr;

    Relations(char *filename){
        freopen(filename, "r", stdin);
        string str,key,value;
        getline(cin, str);
        while(getline(cin, str)){
            istringstream is(str);
            is>>key>>value;
            attr.push_back(value);
        }
        cin.clear();
        fclose(stdin);
    };
};

/**
 * Derived class for the use Matrix from the basic
 * Matrix class.
 */
class USE_Matrix : public Matrix
{
private:
    /**
    * Function to calculate the number of times an attribute is accessed by a query.
     */
    int count_attr(string query, string attr){
        int index = 0, count = 0;
        while( (index = query.find(attr,index)) < query.length() ){
            if(query[index-1] >= 'A'&&query[index-1] <= 'Z' 
            || query[index+attr.length()] >= 'A' && query[index+attr.length()] <= 'Z');
            else    count++;
		    index++;
        }
        return count>0;
    };
public:
    /* USE matrix */
    vector<vector<int>> use;

    USE_Matrix(Quries &q, Relations &r){
        use.resize(q.query.size());
        auto k = use.begin();
        for(auto i = q.query.begin(); i < q.query.end(); i++, k++){
            for(auto j = r.attr.begin(); j < r.attr.end(); j++){
                k->push_back(count_attr(*i, *j));
            }
        }
    };
};

/**
 * Derived class for the Access Frequence Matrix from the basic
 * Matrix class.
 */
class ACC_Matrix : public Matrix
{
public:
    /* Access Frequence matrix */
    vector<vector<int>> acc;

    ACC_Matrix(char *filename){
        freopen(filename, "r", stdin);
        string str,query,frequency;
        getline(cin, str);
        while(getline(cin, str)){
            istringstream is(str);
            vector<int> query_freq;
            is>>query;
            while(is>>frequency){
                query_freq.push_back(atoi(frequency.c_str()));
            }
            acc.push_back(query_freq);
        }
        cin.clear();
        fclose(stdin);
    };
};

/**
 * Derived class for the Attribute Affinity Matrix from the basic
 * Matrix class.
 */
class AA_Matrix : public Matrix
{
private:
    /* Attribute matrix */
    vector<vector<int>> A;
public:
    /* Attribute Affinity matrix */
    vector<vector<int>> AA;

    AA_Matrix(ACC_Matrix &acc, USE_Matrix &use){
        /*
        * Calculating the Attribute matrix A.
        * A[i][k] is the number of times Attribute A[i] is accessed by Query q[k]
        */
        A.resize(use.use[0].size());
        for(auto Ai = A.begin(); Ai < A.end(); Ai++){            
            for(auto k = acc.acc.begin(); k < acc.acc.end(); k++){
                int Aik = 0;
                for(auto j = k->begin(); j < k->end(); j++){
                    Aik += *j;
                }
                Ai->push_back(Aik * use.use[k - acc.acc.begin()][Ai - A.begin()]);
            }
        }
        /*
        * Calculating the Attribute Affinity Matrix.
        */
        AA.resize(A.size());
        for(auto AAi = AA.begin(); AAi < AA.end(); AAi++){
            int sum_i = 0;
            int i = AAi - AA.begin();
            for(auto Aik = A[i].begin(); Aik < A[i].end(); Aik++){
                sum_i += *Aik;
            }
            for(int j = 0; j < A.size(); j++){
                int sum_j = 0, sum_ixj=0;
                for(auto Ajk = A[j].begin(); Ajk < A[j].end();Ajk++){
                    int k = Ajk - A[j].begin();
                    sum_j += *Ajk;
                    sum_ixj += (A[i][k]* (*Ajk));
                }
                if(sum_ixj == 0)
                    AAi->push_back(0);
                else
                    AAi->push_back(ceil(1.0 * sum_ixj / sqrt(sum_i * sum_j)));
            }
        }
    };
};

int main(int argc, char **argv){
    if(argc == 4){
        auto attr = Relations(argv[1]);
        auto query = Quries(argv[2]);
        auto acc = ACC_Matrix(argv[3]);
        auto use = USE_Matrix(query, attr);
        auto aa = AA_Matrix(acc, use);
        aa.printMatrix(aa.AA);
    }
}