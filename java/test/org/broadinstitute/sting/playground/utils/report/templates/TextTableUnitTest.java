package org.broadinstitute.sting.playground.utils.report.templates;

import org.junit.Assert;
import org.junit.Test;

import java.util.ArrayList;

/**
 * test out the text table
 */
public class TextTableUnitTest {
    @Test
    public void testBasicSetCell() {
        TextTable table = new TextTable("name","description",new ArrayList<String>());
        table.setCell(1,5,"entry");
        int entriesSeen = 0;
        ArrayList<ArrayList<String>> rows = table.rows;
        for (int x = 0; x <= 1; x++)
            for (int y = 0; y <= 5; y++)
                if (x == 1 && y == 5) {
                    Assert.assertTrue(rows.get(x).get(y).equals("entry"));
                    entriesSeen++;
                } else
                    Assert.assertTrue(rows.get(x).get(y).equals(""));
        Assert.assertEquals("Incorrect number of entries seen",1,entriesSeen);
    }

    @Test
    public void testBasicSetTwoCells() {
        TextTable table = new TextTable("name","description",new ArrayList<String>());
        table.setCell(1,5,"entry");
        table.setCell(1,1,"entry");
        int entriesSeen = 0;
        ArrayList<ArrayList<String>> rows = table.rows;
        for (int x = 0; x <= 1; x++)
            for (int y = 0; y <= 5; y++)
                if ((x == 1 && y == 5) || (x == 1 && y == 1)) {
                    Assert.assertTrue(rows.get(x).get(y).equals("entry"));
                    entriesSeen++;
                }
                else
                    Assert.assertTrue(rows.get(x).get(y).equals(""));
        Assert.assertEquals("Incorrect number of entries seen",2,entriesSeen);
    }
}
