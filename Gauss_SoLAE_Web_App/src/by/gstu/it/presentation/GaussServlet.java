package by.gstu.it.presentation;

import by.gstu.it.data.MatrixDto;
import by.gstu.it.data.ResultDto;
import by.gstu.it.domain.ClassInjector;
import by.gstu.it.domain.GaussSolver;
import by.gstu.it.domain.WorkTimeCounter;
import com.google.gson.Gson;

import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.util.Arrays;

public class GaussServlet extends HttpServlet {

    private final GaussSolver solver = ClassInjector.provideGaussSolver();
    private final Gson gson = new Gson();

    public void doPost(HttpServletRequest request, HttpServletResponse response)
            throws IOException {

        MatrixDto matrixDto = readJson(request);

        final Double[][] res = new Double[1][];
        WorkTimeCounter timeCounter = new WorkTimeCounter(() -> res[0] = solver.solve(matrixDto));
        timeCounter.makeWork();

        String time = timeCounter.getWorkTimeFormatted();

        String gsonData = gson.toJson(new ResultDto(Arrays.toString(res[0]), time));

        response.setContentType("application/json");
        PrintWriter out = response.getWriter();
        out.print(gsonData);
        out.flush();
    }

    private MatrixDto readJson(HttpServletRequest request) {
        MatrixDto matrixDto = new MatrixDto( new Double[][]{{0.0, 0.0}}, 0);
        try {
            StringBuilder sb = new StringBuilder();
            String s;
            while ((s = request.getReader().readLine()) != null) sb.append(s);

            matrixDto = gson.fromJson(sb.toString(), MatrixDto.class);
        } catch (Exception ignored) {}

        return matrixDto;
    }

}