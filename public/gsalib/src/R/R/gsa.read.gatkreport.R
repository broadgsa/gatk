# Load a table into the specified environment.  Make sure that each new table gets a unique name (this allows one to cat a bunch of tables with the same name together and load them into R without each table overwriting the last.
.gsa.assignGATKTableToEnvironment <- function(tableName, tableHeader, tableRows, tableEnv) {
    d = data.frame(tableRows, row.names=NULL, stringsAsFactors=FALSE);
    colnames(d) = tableHeader;
    
    for (i in 1:ncol(d)) {
        # use the general type.convert infrastructure of read.table to convert column data to R types 
        v = type.convert(d[,i])
        d[,i] = v;
    }
    
    usedNames = ls(envir=tableEnv, pattern=tableName);
    
    if (length(usedNames) > 0) {
        tableName = paste(tableName, ".", length(usedNames), sep="");
    }
    
    assign(tableName, d, envir=tableEnv);
}

# Read a fixed width line of text into a list.
.gsa.splitFixedWidth <- function(line, columnStarts) {
    splitStartStop <- function(x) {
        x = substring(x, starts, stops);
        x = gsub("^[[:space:]]+|[[:space:]]+$", "", x);
        x;
    }
    
    starts = c(1, columnStarts);
    stops = c(columnStarts - 1, nchar(line));
    
    sapply(line, splitStartStop)[,1];
}

# Old implementaton for v0.*
gsa.read.gatkreportv0 <- function(lines) {
    
    tableEnv = new.env();
    
    tableName = NA;
    tableHeader = c();
    tableRows = c();
    version = NA;
    
    for (line in lines) {
        if (length(grep("^##:GATKReport.v", line, ignore.case=TRUE)) > 0) {
            headerFields = unlist(strsplit(line, "[[:space:]]+"));
            
            if (!is.na(tableName)) {
                .gsa.assignGATKTableToEnvironment(tableName, tableHeader, tableRows, tableEnv);
            }
            
            tableName = headerFields[2];
            tableHeader = c();
            tableRows = c();
            
            # For differences in versions see
            #   $STING_HOME/public/java/src/org/broadinstitute/sting/gatk/report/GATKReportVersion.java
            if (length(grep("^##:GATKReport.v0.1[[:space:]]+", line, ignore.case=TRUE)) > 0) {
                version = "v0.1";
                
            } else if (length(grep("^##:GATKReport.v0.2[[:space:]]+", line, ignore.case=TRUE)) > 0) {
                version = "v0.2";
                columnStarts = c();
                
            }
            
        } else if (length(grep("^[[:space:]]*$", line)) > 0 | length(grep("^[[:space:]]*#", line)) > 0) {
            # do nothing
        } else if (!is.na(tableName)) {
            
            if (version == "v0.1") {
                row = unlist(strsplit(line, "[[:space:]]+"));
                
            } else if (version == "v0.2") {
                if (length(tableHeader) == 0) {
                    headerChars = unlist(strsplit(line, ""));
                    # Find the first position of non space characters, excluding the first character
                    columnStarts = intersect(grep("[[:space:]]", headerChars, invert=TRUE), grep("[[:space:]]", headerChars) + 1);
                }
                
                row = .gsa.splitFixedWidth(line, columnStarts);
            }
            
            if (length(tableHeader) == 0) {
                tableHeader = row;
            } else {
                tableRows = rbind(tableRows, row);
            }
        }
    }
    
    if (!is.na(tableName)) {
        .gsa.assignGATKTableToEnvironment(tableName, tableHeader, tableRows, tableEnv);
    }
    
    gatkreport = as.list(tableEnv, all.names=TRUE);
}

# Load all GATKReport v1 tables from file
gsa.read.gatkreportv1 <- function(lines) {
  #print("loading with optimized v1 reader")
  nLines = length(lines)
  tableEnv = new.env();
  
  tableName = NA;
  tableHeader = c();
  tableRows = NULL;
  version = "";
  rowCount = 0
  headerRowCount = -1;
  
  finishTable <- function() {
    if ( rowCount == 1 )
      # Workaround to avoid collapsing into an unstructured vector when 
      # there's only 1 row
      sub <- t(as.matrix(tableRows[1:rowCount,]))
    else
      sub <- tableRows[1:rowCount,]
    .gsa.assignGATKTableToEnvironment(tableName, tableHeader, sub, tableEnv);
  }
  
  for (line in lines) {
    
    if (length(grep("^#:GATKReport.v1", line, ignore.case=TRUE)) > 0) {
      version = "v1.0";
      headerRowCount = 0;
    }
    
    if ( (headerRowCount %% 2 == 1) && (version == "v1.0") ) {
      #print("Trying to start a table with line:");
      #print(line);
      
      #Get table header
      headerFields = unlist(strsplit(line, ":"));
      
      if (!is.na(tableName)) {
        finishTable()
      }
      
      tableName = headerFields[3];
      tableHeader = c();
      tableRows = NULL
      rowCount = 0
      
      columnStarts = c();
    }
    
    if (length(grep("^#:GATKTable", line, ignore.case=TRUE)) > 0) {
      headerRowCount = headerRowCount+1;
      #print("Header Row count is at:")
      #print(headerRowCount);
    } else if (!is.na(tableName)) {
      if ( version == "v1.0") {
        if (length(tableHeader) == 0) {
          headerChars = unlist(strsplit(line, ""));
          # Find the first position of non space characters, excluding the first character
          columnStarts = intersect(grep("[[:space:]]", headerChars, invert=TRUE), grep("[[:space:]]", headerChars) + 1);
          tableRows = matrix(nrow=nLines, ncol=length(columnStarts)+1);
        }
        
        row = .gsa.splitFixedWidth(line, columnStarts);
      }
      
      if (length(tableHeader) == 0) {
        tableHeader = row;
      } else if ( nchar(line) > 0 ) {
        rowCount = rowCount + 1
        tableRows[rowCount,] <- row
      }
    }
  }
  
  if (!is.na(tableName)) {
    finishTable()
  }
  
  gatkreport = as.list(tableEnv, all.names=TRUE);
}

# Load all GATKReport tables from a file
gsa.read.gatkreport <- function(filename) {
    con = file(filename, "r", blocking = TRUE);
    lines = readLines(con);
    close(con);
    
    # get first line
    line = lines[1];
    
    if (length(grep("^#:GATKReport.v1", line, ignore.case=TRUE)) > 0) {
        gsa.read.gatkreportv1(lines)
    }
    else if (length(grep("^##:GATKReport.v0", line, ignore.case=TRUE)) > 0) {
        gsa.read.gatkreportv0(lines)
    }
}
