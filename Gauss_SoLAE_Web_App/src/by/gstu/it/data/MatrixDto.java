package by.gstu.it.data;

import java.util.Arrays;
import java.util.Objects;

public class MatrixDto {
    public Double[][] matrix;
    public int mpi;

    public MatrixDto(Double[][] matrix, int mpi) {
        this.matrix = matrix;
        this.mpi = mpi;
    }

    public MatrixDto() {}

    public Double[][] getMatrix() {
        return matrix;
    }

    public void setMatrix(Double[][] matrix) {
        this.matrix = matrix;
    }

    public int getMpi() {
        return mpi;
    }

    public void setMpi(int mpi) {
        this.mpi = mpi;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        MatrixDto matrixDto = (MatrixDto) o;
        return mpi == matrixDto.mpi &&
                Arrays.equals(matrix, matrixDto.matrix);
    }

    @Override
    public int hashCode() {
        int result = Objects.hash(mpi);
        for (Double[] arr : matrix)
            result = 31 * result + Arrays.hashCode(arr);
        return result;
    }

    @Override
    public String toString() {
        StringBuilder result = new StringBuilder("");
        result.append(matrix.length).append("\n");

        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                result.append(matrix[i][j]);
                if (j != matrix[i].length - 1)
                    result.append(" ");
            }
            if (i != matrix.length - 1)
                result.append("\n");
        }

        return result.toString();
    }
}
