package by.gstu.it.domain;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;

public class ProcessHelper {

    public static Process run(String... commands) {
        ArrayList<String> cmd = new ArrayList<>();
        cmd.add("cmd");
        cmd.add("/c");
        cmd.addAll(Arrays.asList(commands));

        ProcessBuilder builder = new ProcessBuilder( cmd);
        builder.redirectErrorStream(true);

        Process process = null;

        try {
            process = builder.start();

            try(BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()))) {
                String line = null;
                do {
                    line = reader.readLine();
                    System.out.println(line);
                } while (line != null);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return process;
    }

    public static Process runAndWait(String... commands) {
        Process process = run(commands);
        try {
            process.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return process;
    }

}
