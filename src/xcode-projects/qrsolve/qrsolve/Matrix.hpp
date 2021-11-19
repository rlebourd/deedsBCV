//
//  Matrix.hpp
//  qrsolve
//
//  Created by Robert LeBourdais on 11/18/21.
//

#ifndef Matrix_hpp
#define Matrix_hpp

#include <stdio.h>
#include <array>

namespace lin_alg {

    template <size_t Rows, size_t Cols>
    class Matrix {
    public:
        static_assert(Rows > 0, "Rows must be greater than zero");
        static_assert(Cols > 0, "Cols must be greater than zero");

        Matrix(std::array<double, Rows*Cols> elems){
            for (int i = 0; i < Rows*Cols; i++){
                elements[i] = elems.at(i);
            }
        }

        Matrix(){
            for (int i = 0; i < Rows*Cols; i++){
                elements[i] = 0;
            }
        }

        static Matrix<Rows, Cols> identityMatrix(){
            static_assert(Rows == Cols, "Rows must equal columns");
            auto matrix = Matrix<Rows, Cols>{};
            for (int i = 0; i < Rows; i++){
                matrix.at(i, i) = 1;
            }
            return matrix;
        }
        
        static Matrix<Rows, Cols> zeroMatrix(){
            return Matrix<Rows, Cols>{};
        }
        
        template <size_t ColsRHS>
        Matrix<Rows, ColsRHS> times(const Matrix<Cols, ColsRHS> &rhs) const {
            Matrix<Rows, ColsRHS> result;
            
            // implement matrix multiplication
            for (int i = 0; i < Rows; i++){
                for (int j = 0; j < ColsRHS; j++){
                    double aij = 0;
                    for (int n = 0; n < Cols; n++){
                        const auto ain = at(i, n);
                        const auto bnj = rhs.at(n, j);
                        aij += ain * bnj;
                    }
                    result.at(i, j) = aij;
                }
            }
            
            return result;
        }
        
        Matrix<Rows, Cols> times(double scalar) const {
            Matrix<Rows, Cols> matrix;
            for (int i = 0; i < Rows*Cols; i++){
                matrix.elements[i] = scalar*elements[i];
            }
            return matrix;
        }
        
        bool operator==(const Matrix<Rows, Cols> &rhs) const {
            for (int i = 0; i < Rows*Cols; i++){
                if (elements[i] != rhs.elements[i]){
                    return false;
                }
            }
            return true;
        }
        
        double &at(size_t row, size_t col) {
            return elements[indexOfRowCol(row, col)];
        }

        double at(size_t row, size_t col) const {
            return elements[indexOfRowCol(row, col)];
        }
        
        double norm() const {
            double norm = 0;
            for (int i = 0; i < Rows*Cols; i++){
                norm += elements[i]*elements[i];
            }
            return norm;
        }

        Matrix<Rows, Cols> normalized() const {
            return this->times(norm());
        }

    private:
        size_t indexOfRowCol(size_t row, size_t col) const {
            return Cols*row + col;
        }
        
    private:
        std::array<double, Rows*Cols> elements{0};
    };

    using Vector3D  = Matrix<3, 1>;
    using Vector4D  = Matrix<4, 1>;
    using Matrix3x3 = Matrix<3, 3>;
    using Matrix4x4 = Matrix<4, 4>;

}

#endif /* Matrix_hpp */
