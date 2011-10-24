gsa.plot.venn <-
function(a, b, c=0, a_and_b, a_and_c=0, b_and_c=0,
                     col=c("#FF6342", "#63C6DE", "#ADDE63"),
                     pos=c(0.20, 0.20, 0.80, 0.82),
                     debug=0
                    ) {
    library(png);
    library(graphics);

    # Set up properties
    for (i in 1:length(col)) {
        rgbcol = col2rgb(col[i]);
        col[i] = sprintf("%02X%02X%02X", rgbcol[1], rgbcol[2], rgbcol[3]);
    }

    chco = paste(col[1], col[2], col[3], sep=",");
    chd = paste(a, b, c, a_and_b, a_and_c, b_and_c, sep=",");

    props = c(
        'cht=v',
        'chs=525x525',
        'chds=0,10000000000',
        paste('chco=', chco, sep=""),
        paste('chd=t:', chd, sep="")
    );
    proplist = paste(props[1], props[2], props[3], props[4], props[5], sep='&');

    # Get the venn diagram (as a temporary file)
    filename = tempfile("venn");
    cmd = paste("wget -O ", filename, " 'http://chart.apis.google.com/chart?", proplist, "' > /dev/null 2>&1", sep="");

    if (debug == 1) {
        print(cmd);
    }
    system(cmd);

    # Render the temp png file into a plotting frame
    a = readPNG(filename);
    
    plot(0, 0, type="n", xaxt="n", yaxt="n", bty="n", xlim=c(0, 1), ylim=c(0, 1), xlab="", ylab="");
    if (c == 0 || a >= b) {
        rasterImage(a, pos[1], pos[2], pos[3], pos[4]);
    } else {
        rasterImage(a, 0.37+pos[1], 0.37+pos[2], 0.37+pos[3], 0.37+pos[4], angle=180);
    }

    # Clean up!
    unlink(filename);
}

