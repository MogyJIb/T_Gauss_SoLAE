package by.gstu.it.domain;


import java.text.NumberFormat;
import java.util.Locale;

public class WorkTimeCounter {
    private Workable workableObj;
    private long workTime;

    public WorkTimeCounter(Workable work) {
        this.workableObj = work;
    }

    public void makeWork() {
        long before = System.nanoTime();
        workableObj.work();
        long after = System.nanoTime();
        this.workTime = after - before;
    }

    public long getWorkTime() {
        return this.workTime;
    }

    public String getWorkTimeFormatted() {
        return nanoSecondsToString(workTime);
    }

    @Override
    public String toString() {
        return nanoSecondsToString(workTime);
    }

    public static String nanoSecondsToString(long nanoSeconds) {
        StringBuilder result = new StringBuilder("");
        if (nanoSeconds == 0)
            return "0";

        long sec = nanoSeconds / 1000000000,
                ml_sec = (nanoSeconds / 1000000) % 1000,
                mk_sec = (nanoSeconds / 1000) % 1000,
                n_sec = (nanoSeconds) % 1000;


        result.append(String.format("%4ds %4dmls",
                sec, ml_sec));
        return result.toString();
    }

    public static Double nanoSecondsToSeconds(long nanoSeconds, int accuracy) {
        return roundDouble(nanoSeconds / 1000000.0, accuracy);
    }

    /**
     * Round double with accuracy
     * max fraction digits.
     *
     * @param value    double number
     * @param accuracy makeStatistic of number after comma
     * @return converted number
     */
    public static double roundDouble(double value, int accuracy) {
        NumberFormat numberFormat = NumberFormat.getInstance(Locale.ENGLISH);
        numberFormat.setMaximumFractionDigits(accuracy);
        numberFormat.setGroupingUsed(false);

        String newValue = numberFormat.format(value);
        return Double.parseDouble(newValue);
    }
}
