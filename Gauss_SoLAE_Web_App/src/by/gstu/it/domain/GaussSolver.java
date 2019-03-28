package by.gstu.it.domain;

import by.gstu.it.data.MatrixDto;
import java.util.Arrays;


public class GaussSolver {
    private static final String INPUT = "input.txt";
    private static final String OUTPUT = "output.txt";
    private static final String SOLVER = "Gauss_SoLAE_MPI.exe";

    private static final int PROCESS_COUNT = 5;
    private static int USE_MPI = 1;

    public Double[] solve(MatrixDto matrixDto) {
        USE_MPI = matrixDto.mpi;

        FileHelper.writeRes(INPUT, matrixDto.toString());
        ProcessHelper.runAndWait(buildCommands());
        String content = FileHelper.readRes(OUTPUT);
        return parseResult(content);
    }

    private Double[] parseResult(String content) {
        return Arrays.stream(content.split(" ")).map(Double::parseDouble).toArray(Double[]::new);
    }

    private String[] buildCommands() {
        String path = this.getClass().getClassLoader().getResource(SOLVER).getFile();
        if (path.startsWith("/") || path.startsWith("\\"))
            path = path.substring(1);
        path = path.substring(0, path.length() - SOLVER.length());

        return new String[]{ path + SOLVER, path };
        //return new String[]{ "mpiexec", "-n", String.valueOf(PROCESS_COUNT), path + SOLVER, path, String.valueOf(USE_MPI) };
    }

}
