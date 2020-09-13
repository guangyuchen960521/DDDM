#include <iostream>
#include <fstream>
#include <cstdlib>
#include<sstream>
#include<string>

using namespace std;

int MATRIX_SIZE;


class Matrix {
protected:
 
    int **_ptr;


    int _length;

 
    int* _tempC;

public:

    Matrix(int n) {
        _length = n;
        _tempC = NULL;

        // Initialize the 2D matrix with zeros
        _ptr = (int**) calloc(n+1, sizeof(int*));
        for(int i = 0; i < n+1; i++) {
            _ptr[i] = (int*) calloc(n+1, sizeof(int));
        }

        // Number the colums in the top row. 0th index
        for(int i=1; i<= n; i++) {
            _ptr[0][i] = i;
        }
    }
    Matrix(){};
  
    class Proxy {
    private:
        int* _array;

        // Index of the first [] level
        int idx;
    public:
        Proxy(int* arr) : _array(arr) {}

        int& operator[] (int index) {
            idx = index;
            return _array[index];
        }

        void operator= (int rhs) {
            _array[idx] = rhs;
        }
    };


    Proxy operator[] (int index) {
        return Proxy(_ptr[index]);
    }


    int* getColumn(int colIndex) {
        free(_tempC);
        _tempC = (int*) calloc(_length+1, sizeof(int));
        for(int i=0; i<=_length; i++) {
            _tempC[i] = _ptr[i][colIndex];
        }

        return _tempC;
    }

 
    void printMatrix() {
        for(int i=0; i<=_length; i++) {
            if(i==0)
                cout<<"\t";
            else    cout<<"A"<<_ptr[0][i]<<"\t";
            for(int j=1; j<=_length; j++) {
                if(i==0){
                    cout<<"A";
                }
                cout<<_ptr[i][j]<<"\t";
            }
            cout<<endl;
        }
    }
};



class AttributeMatrix : public Matrix {
private:
    /**
     * Function to calculate the number of attribute.
     */
    int get_Matrix_size(char *filename){
        MATRIX_SIZE = 0;
        freopen(filename, "r", stdin);
        string str;
        getline(cin, str);
        istringstream is(str);
        int value;
        while(is>>value){
            MATRIX_SIZE++;
        }
        cin.clear();
        fclose(stdin);
        return MATRIX_SIZE;
    }
public:
    AttributeMatrix(char *filename) : Matrix(get_Matrix_size(filename)){}
    /**
     * Function to calculate the bond energy between two columns.
     */
    int calculateBond(int left, int right) {
        int sum = 0;
        for(int i = 1; i<=_length; i++) {
            sum = sum + (_ptr[i][left] * _ptr[i][right]);
        }

        return sum;
    }
};



class ClusteredMatrix : public Matrix {
private:
    // The index of the rightmost column filled in the clustered matrix.
    int _rightmostIndex;

    int _maxLeft;

    // contains the column number of the column that is being placed. (A2 in example)
    int _maxMid;


    int _maxRight;

public:
    ClusteredMatrix(int n) : Matrix(n) {
        _rightmostIndex = 2;
        _maxRight = 0;
        _maxMid = 0;
        _maxLeft = 0;

        for(int i=1; i<= n; i++) {
            _ptr[0][i] = 0;
        }
    }

 
    void copyToColumn(int colIndex, int* array) {
        for(int i=0; i<=_length; i++) {
            _ptr[i][colIndex] = array[i];
        }
    }


    int getRightmostIndex() {
        return _rightmostIndex;
    }

  
    void recordPlacement(int left, int mid, int right) {
        _maxLeft = left;
        _maxMid = mid;
        _maxRight = right;
    }

    /**
     * Function for placing a column in the clustered matrix from the affinity matrix
     */
    void placeColumnFrom(AttributeMatrix& AA) {

        // Handle placement in leftmost case
        if(_maxLeft == 0) {
            // Shift columns to the right
            for(int i = _rightmostIndex+1; i>1; i--) {
                copyToColumn(i, getColumn(i-1));
            }

            // Place the column.
            copyToColumn(1, AA.getColumn(_maxMid));

            // Increment the index of the rightmost filled column
            _rightmostIndex++;
            return;
        }

        // For other cases, find the index (start) after which to place the column.
        int start;
        for(start=1; start<=_length; start++) {
            if(_ptr[0][start] == _maxLeft)
                break;
        }

        // If start is the rightmost filled index, this is the case when placing in after the
        // rightmost column.
        if(start == _rightmostIndex) {
            _rightmostIndex++;
            copyToColumn(_rightmostIndex, AA.getColumn(_maxMid));
            return;
        }

        // Handle placement in between two columns.
        // Do shifting.
        for(int i = _rightmostIndex+1; i>start+1; i--) {
            copyToColumn(i, getColumn(i-1));
        }

        // Place the column from the affinity matrix.
        copyToColumn(start+1, AA.getColumn(_maxMid));

        // Increment the rightmost filled column index.
        _rightmostIndex++;
    }

    /**
     * Function to make the clustered matrix symmetrical in the final step.
     */
    void makeSymmetrical() {
        int **temp;

        temp = (int**) calloc(_length+1, sizeof(int*));
        for(int i = 0; i < _length+1; i++) {
            temp[i] = (int*) calloc(_length+1, sizeof(int));
        }


        for(int i=1; i<=_length; i++) {
            int row = _ptr[0][i];
            temp[0][i] = row;

            for(int j=1; j<=_length; j++) {
                temp[i][j] = _ptr[row][j];
            }
        }

        // free the previous pointer
        for(int i=0; i<=_length; i++) {
            free(_ptr[i]);
        }
        free(_ptr);

        // Point to new array
        _ptr = temp;
    }
};


/**
 * Function to calculate the bond between two columns
 */
int bond(int left, int right, AttributeMatrix& AA) {
    return AA.calculateBond(left, right);
}


/**
 * Function to calculate the contribution of a certain configuration of columns.
 */
int calculateContribution(int left, int middle, int right, AttributeMatrix& AA) {
    // if leftmost case
    if(left == 0) {
        return 2*bond(middle, right, AA);
    }

    // if rightmost case
    if(right == middle + 1) {
        return 2*bond(left, middle, AA);
    }

    // if placed in between two columns
    return 2*bond(left, middle, AA) + 2*bond(middle, right, AA) - 2*bond(left, right, AA);
}


/**
 * Function to perform the BEA algorithm.
 */
void doBea(AttributeMatrix& AA, ClusteredMatrix& CA) {

    // Copy the first and second columns from the
    // Attribute Affinity Matrix
    CA.copyToColumn(1, AA.getColumn(1));
    CA.copyToColumn(2, AA.getColumn(2));

    int index = 3;

    while(index <= MATRIX_SIZE) {
        // Declare and Initialize required variables
        int contrib = 0;
        int maxContribution = 0;

        // Calculate placement that gives highest contribution
        maxContribution = calculateContribution(CA[0][0], index, CA[0][1], AA);
        CA.recordPlacement(CA[0][0], index, CA[0][1]);
        for(int i = 1; i < index; i++) {
            contrib = calculateContribution(CA[0][i-1], index, CA[0][i], AA);

            if(contrib > maxContribution) {
                maxContribution = contrib;
                CA.recordPlacement(CA[0][i-1], index, CA[0][i]);
            }
        }

        contrib = calculateContribution(CA[0][index-1], index, index+1, AA);

        if(contrib > maxContribution) {
            maxContribution = contrib;
            CA.recordPlacement(CA[0][index-1], index, index+1);
        }

        // Place the column that has the highest contribution
        CA.placeColumnFrom(AA);
        index++;
    }
    CA.makeSymmetrical();
    CA.printMatrix();
}

int main(int argc, char** argv) {
    if(argc == 2){
        // The Attribute Affinity Matrix
        AttributeMatrix AA(argv[1]);

        // The Clustered Affinity Matrix
        ClusteredMatrix CA(MATRIX_SIZE);

        // Populate the Attribute Affinity Matrix from the input file.
        freopen(argv[1], "r", stdin);
        for (int i = 1; i <= MATRIX_SIZE; i++)
        {
            for (int j = 1; j <= MATRIX_SIZE; j++)
            {
                cin >> AA[i][j];
            }
            
        }
        cin.clear();
        fclose(stdin);

        // Bond Energy Algorithm
        doBea(AA, CA);
    }
    
}
