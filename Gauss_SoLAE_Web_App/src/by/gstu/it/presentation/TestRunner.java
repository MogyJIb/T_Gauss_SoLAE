package by.gstu.it.presentation;

import by.gstu.it.data.MatrixDto;
import by.gstu.it.domain.ClassInjector;
import by.gstu.it.domain.FileHelper;
import by.gstu.it.domain.GaussSolver;
import by.gstu.it.domain.ProcessHelper;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;
import java.util.function.IntFunction;
import java.util.stream.Collectors;

public class TestRunner {

    public static void main(String[] args) {
        testSOlve();
    }

    private static void testSOlve() {
        GaussSolver solver = ClassInjector.provideGaussSolver();
        System.out.println(Arrays.toString(solver.solve(new MatrixDto(new Double[][]{
                {3.0, 2.0, 3.0, 1.0},
                { 4.0, 4.0, 3.0, 1.0},
                { 1.0, 4.0, 4.0, 2.0}
        }, 0))));
    }

    private static void makeTestFile() {
        String content = FileHelper.readRes("test1.txt");

        List<List<Double>> matrix= Arrays.stream(content.split("\n"))
                                            .map(line ->
                                                    Arrays.stream(line.trim().split(" "))
                                                            .map(Double::parseDouble)
                                                            .collect(Collectors.toList())
                                            )
                                            .collect(Collectors.toList());

        String contentB = FileHelper.readRes("testB.txt");
        List<Double> b = Arrays.stream(contentB.split("\n"))
                                    .map(Double::parseDouble)
                                    .collect(Collectors.toList());


        for (int i = 0; i < matrix.size(); i++) {
            matrix.get(i).add(b.get(i));
        }

        FileHelper.writeRes("test2.txt", new MatrixDto(matrix.stream().map(
                doubles -> doubles.stream().toArray(Double[]::new)
        ).toArray(Double[][]::new), 0).toString());




    }

}
