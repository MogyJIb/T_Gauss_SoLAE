package by.gstu.it.domain;

public class ClassInjector {

    private static GaussSolver gaussSolver;

    public static GaussSolver provideGaussSolver() {
        if (gaussSolver == null) gaussSolver = new GaussSolver();
        return gaussSolver;
    }

}
