package by.gstu.it.domain;

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;

public class FileHelper {

    public static void writeRes(String resName, String data) {
        try(BufferedWriter writer = new BufferedWriter(new FileWriter(new File(toURI(resName))))) {
            writer.write(data);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static String readRes(String resName) {
        StringBuilder text = new StringBuilder("");

        try(BufferedReader reader = new BufferedReader(new FileReader(new File(toURI(resName))))) {

            String line = null;
            while (true) {
                line = reader.readLine();
                if (line == null) break;

                text.append(line).append("\r\n");
            }

            if (text.length() > 2)
                text.delete(text.length() - 2, text.length());

        } catch (Exception e) {
            e.printStackTrace();
        }

        return text.toString();
    }

    private static URI toURI(String resName) throws URISyntaxException {
        return FileHelper.class.getClassLoader().getResource(resName).toURI();
    }

}
