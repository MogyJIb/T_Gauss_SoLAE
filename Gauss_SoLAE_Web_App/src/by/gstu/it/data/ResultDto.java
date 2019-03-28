package by.gstu.it.data;

import java.io.Serializable;
import java.util.Objects;

public class ResultDto implements Serializable {
    public String solution;
    public String time;

    public ResultDto() {}

    public ResultDto(String solution, String time) {
        this.solution = solution;
        this.time = time;
    }

    public String getSolution() {
        return solution;
    }

    public void setSolution(String solution) {
        this.solution = solution;
    }

    public String getTime() {
        return time;
    }

    public void setTime(String time) {
        this.time = time;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ResultDto resultDto = (ResultDto) o;
        return Objects.equals(solution, resultDto.solution) &&
                Objects.equals(time, resultDto.time);
    }

    @Override
    public int hashCode() {
        return Objects.hash(solution, time);
    }

    @Override
    public String toString() {
        return "ResultDto{" +
                "solution=" + solution +
                ", time=" + time +
                '}';
    }
}
