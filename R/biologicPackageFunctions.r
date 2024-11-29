###############################################################################
## Make heatmap                                                              ##
###############################################################################
# Make heatmap function                                                       #
###############################################################################

make.hm <- function(
    m.df1,
    filename = "heatmap",
    k.number = 10,
    n.colors = 1000,
    hclust.method = "complete",
    dist.method = "euclidean",
    main = "",
    Colv = TRUE,
    showRowNames = FALSE,
    showColNames = TRUE,
    plotSeparationLines = FALSE
) {
    # library(RColorBrewer)
    # library(gplots)

    hclustfunc <- function(x) hclust(x, method = hclust.method)

    distfunc <- function(x) dist(x, method = dist.method)

    d <- distfunc(m.df1)

    fit <- hclustfunc(d)

    clusters <- cutree(fit, k = k.number)

    nofclust.height <- length(unique(as.vector(clusters)))

    a = max(abs(m.df1))

    breaks = seq((-1 * a), a, length.out = n.colors)


    ## Retired 20180228 ##
    #hmcols <- colorRampPalette(c("blue", "white", "red"))(n = (length(breaks) - 1))
    ## Start New:
    # blueWhiteRedVec <- rev(
    #    colorRampPalette(brewer.pal(9, "RdBu"))(3)
    # )

    blueWhiteRedVec <- c("#3060cf", "#fffbbc","#c4463a")

    hmcols <- grDevices::colorRampPalette(
        blueWhiteRedVec
    )(n = (length(breaks) - 1))

    selcol <- grDevices::colorRampPalette(brewer.pal(12, "Set3"))
    selcol2 <- grDevices::colorRampPalette(brewer.pal(9, "Set1"))
    clustcol.height = selcol2(nofclust.height)

    if (filename != ""){
        pdf(
            paste(
                filename,
                "pdf",
                sep = "."
            )
        )
    } else {
        pdf("temp.pdf")
    }

    if (showRowNames){
        labRowVec = row.names(m.df1)
    } else {
        labRowVec = rep("", nrow(m.df1))
    }

    if (showColNames){
        labColVec = colnames(m.df1)
    } else {
        labColVec = rep("", length(colnames(m.df1)))
    }

    if (plotSeparationLines) {
        colsep = c(0: ncol(m.df1))
        rowsep = c(0: nrow(m.df1))
    } else {
        colsep = c(0, ncol(m.df1))
        rowsep = c(0, nrow(m.df1))
    }

    hm = gplots::heatmap.2(
        m.df1,
        trace = "none",
        dendrogram = "both",
        density.info = "none",
        keysize = 1,
        key = TRUE,
        Colv = Colv,
        hclust = hclustfunc, distfun = distfunc, col = hmcols,
        symbreak = T,
        labRow = labRowVec,
        labCol = labColVec,
        RowSideColors = clustcol.height[clusters],
        margins = c(10, 10),
        cexCol = 1,
        cexRow = 0.5,
        srtCol = 45,
        srtRow = 0,
        main = main,
        breaks = breaks,
        sepcolor = "black",
        sepwidth = c(5e-04, 5e-05),
        colsep = colsep,
        rowsep = rowsep
    )

    if (filename != ""){
        dev.off()
    } else {
        dev.off()
        if (file.exists("temp.pdf")){
            unlink("temp.pdf")
        }
    }
    sorted = m.df1[
        match(
            rev(
                labels(hm$rowDendrogram)),
            rownames(m.df1)
        ),
    ]

    sorted = sorted[, hm$colInd]
    if (filename != ""){
        pdf(paste(filename, "colorkey.pdf", sep = "."))
    } else {
        pdf("temp.pdf")
    }
    plot.new()
    par(lend = 1)
    legend("topleft", legend = 1:nofclust.height, col = clustcol.height,
           lty = 1, lwd = 10)
    if (filename != ""){
        dev.off()
    } else {
        dev.off()
        if (file.exists("temp.pdf")){
            unlink("temp.pdf")
        }
    }
    df.res = list(sorted = sorted, clusters = clusters)
    return(df.res)
    #return(clusters)
}
## End make heatmap                                                          ##
###############################################################################



#' @export

createAndFormatExcelOutputFiles <- function(
    obj,
    metaCoreCountFilter = 1,
    customOutputCols = NULL,
    addedOutputCols = NULL
){
    ###############################################################################
    ## Create Excel output table                                                 ##
    if (length(customOutputCols) > 0){
        outCols <- customOutputCols
    } else {
        outCols <- c(
            obj@parameterList$geneIDcolumn,
            obj@parameterList$primaryAlignmentGeneID,
            "gene_description",
            "gene_type",
            names(obj@databaseTable)[grep("contrast_", names(obj@databaseTable))],
            names(obj@databaseTable)[grep("norm_counts_", names(obj@databaseTable))],
            names(obj@databaseTable)[grep("raw_counts_", names(obj@databaseTable))],
            "count_cut_off",
            "CoVar"
        )
    }

    if (length(addedOutputCols) > 0){
        outCols <- c(
            outCols,
            addedOutputCols
        )
    }

    outCols <- outCols[outCols %in% names(obj@databaseTable)]

    dfOutput <- unique(obj@databaseTable[,outCols])

    ## Rename columns ##
    names(dfOutput) <- gsub("norm_counts_", "", names(dfOutput))
    comparisons <- names(obj@dfDesign)[grep("comp_", names(obj@dfDesign))]

    for (i in 1:length(comparisons)){
        names(dfOutput) <- gsub(
            paste0("contrast_", i, "_"),
            "",
            names(dfOutput)
        )
    }


    outPutFN <- paste0(obj@parameterList$outputDir, obj@parameterList$project_id,".result.table.txt")

    write.table(
        dfOutput,
        outPutFN,
        row.names = FALSE,
        sep="\t"
    )

    ## Create Excel file ##
    #library(openxlsx)

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, paste0(obj@parameterList$project_id, "_full_DGE_result_list"))
    openxlsx::freezePane(wb, paste0(obj@parameterList$project_id, "_full_DGE_result_list") ,  firstActiveRow = 2)

    ## Filter is inactivated, as it does not appear to be compatible with the current version of Excel
    #addFilter(wb, 1, row = 1, cols = 1:ncol(dfOutput))

    ## Style headers ##
    hs1 <- openxlsx::createStyle(
        fontColour = "#ffffff",
        fgFill = "#000000",
        halign = "CENTER",
        textDecoration = "Bold"
    )

    openxlsx::writeData(wb, 1, dfOutput, startRow = 1, startCol = 1, headerStyle = hs1)

    openxlsx::saveWorkbook(
        wb,
        gsub(".txt", ".xlsx", outPutFN) ,
        overwrite = TRUE
    )

    ## Done creating Excel output table                                          ##
    ###############################################################################

    ###############################################################################
    ## Create metacore table                                                     ##

    ## Apply minimal filtering ##
    df.metacore <- obj@databaseTable[obj@databaseTable$count_cut_off >= metaCoreCountFilter,]

    ## select padj and logFC columns only ##
    sel.vec <- names(df.metacore)[grep("contrast_", names(df.metacore))]
    sel.vec <- sel.vec[-grep("stat", sel.vec)]
    sel.vec <- sel.vec[-grep("lg10p", sel.vec)]
    ## Remove contrast_x_ prefix > remove front 11 characters
    sel.vec <- append(obj@parameterList$geneIDcolumn, sel.vec)
    df.metacore <- unique(df.metacore[,sel.vec])

    ## Remove contrast_x_ tag from column labels ##
    comparisons <- names(obj@dfDesign)[grep("comp_", names(obj@dfDesign))]
    for (i in 1:length(comparisons)){
        names(df.metacore) <- gsub(
            paste0("contrast_", i, "_"),
            "",
            names(df.metacore)
        )
    }

    outPutFN <- paste0(obj@parameterList$outputDir, obj@parameterList$project_id,".metacore.input.file.txt")

    write.table(
        df.metacore,
        outPutFN,
        row.names = FALSE,
        sep="\t"
    )

    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, paste0(obj@parameterList$project_id, "_metacore_input_file"))

    ## Style headers ##
    hs1 <- openxlsx::createStyle(
        fgFill = "#4F81BD",
        halign = "CENTER",
        textDecoration = "Bold",
        border = "Bottom",
        fontColour = "white"
    )

    openxlsx::writeData(wb, 1, dfOutput, startRow = 1, startCol = 1, headerStyle = hs1)

    openxlsx::saveWorkbook(
        wb,
        gsub(".txt", ".xlsx", outPutFN) ,
        overwrite = TRUE
    )


    ## Done creating Metacore input file                                         ##
    ###############################################################################

    ###############################################################################
    ## Metacore analysis                                                         ##

    ## Open file in excel and save as metacore.input.file.xls (2003)
    #print("Perform enrichment analysis by subsetting data to logFC cut off: 1, padj 0.05")
    ## Select Workflows & Reports
    ## Select Enrichment analysis
    ## Set threshold: Threshold 1 p-value 0.05 direction both
    ## Run analysis
    ## If necessary, repeat with lower theshold
    ## If successful hit Get report button and safe as
    #print(paste0(project.code, ".metacore.results.enrichment.analysis.xls"))

    ## Next do transcription factor target analysis ##
    ## Select One-click Analysis > Transcription Factors
    ## Set FDR threshold to 0.05
    #print(paste0("Export result table as: ", project.code, ".metacore.TF.analysis.", "logFC_nonAligned_vs_aligned"))

    ## Save results as p111.metacore.result.

    ## For selected TF targets, export MC results an curate into project category ##
    ## Select transcription factor of interest (Object name)
    ## Limit selection: Direction: Outgoing Effect Activation Mechanism influence on expression and transcription
    ## regulation >> Aplly >> Build Network
    ## Select additional options
    ## Pre filters Interaction types transcription regulation
    ## Additional options Directions: Downstream
    ## Hit build network
    ## Select all >> File >> Export >> safe as [TFname.mc.targets.xls]


    ## End Module add metacore results                                           ##
    ###############################################################################
    print("Excel output files create and depoisted in project/outputs folder")
}

## End: (7c) Create Excel output tables                                      ##
###############################################################################


###############################################################################
## Kill connections                                                          ##



killDbConnections <- function () {
    #library(RMySQL)
    all_cons <- dbListConnections(RMySQL::MySQL())
    print(all_cons)
    for(con in all_cons)
        res <- DBI::dbDisconnect(con)
    #print(paste(length(all_cons), " connections killed."))
}
##                                                                           ##
###############################################################################


###############################################################################
# Function upload pca.table.to.db()                                           #
###############################################################################

#' @export
upload.pca.table.to.db <- function(
    df.pca,
    host = "www.biologic-db.org",
    port = 3306,
    prim.data.db = "vtl_data",
    password = "",
    db.user = "boeings",
    PCAdbTableName = "P79_VTL_ES_PCA"){


    ## Add little logic to compensate for missing port specification ##
    if (host == "clvd1-db-u-p-31" && port == 3306){
        port <- 6008
    }

    # Create color pallets
    if (length(grep("sample.group_colors", names(df.pca))) == 0){
        sample.group.vec <- names(df.pca)[grep("sample_group", names(df.pca))]
        sample.group.color.vec <- gsub("sample_group", "sample_group_colors", sample.group.vec)
        for (i in 1:length(sample.group.vec)){
            group.size <- length(unique(df.pca[,sample.group.vec[i]]))
            #library(RColorBrewer)
            #selcol <- colorRampPalette(brewer.pal(9,"YlOrBr"))
            #library(scales)
            group.cols <- scales::hue_pal()(group.size)
            assign(sample.group.color.vec[i], group.cols)
            assign(sample.group.vec[i], unique(df.pca[,sample.group.vec[i]]))
            df.col.pca <- data.frame(get(sample.group.vec[i]), get(sample.group.color.vec[i]))
            names(df.col.pca) <- c(sample.group.vec[i], sample.group.color.vec[i])
            df.pca <- merge(
                df.pca,
                df.col.pca,
                by.x = sample.group.vec[i],
                by.y = sample.group.vec[i],
                all = TRUE
            )
        }
    }
    # Adjust pca column names if neceessary
    if (length(grep("^pca", names(df.pca))) > 0){
        names(df.pca) <- gsub("^pca", "PC", names(df.pca))
    }

    if (length(grep("^PCA", names(df.pca))) > 0){
        gsub("^PCA", "PC", names(df.pca))
    }

    #library(RMySQL)
    #df.pca should contain three columns: sample_name, sample.group, sample.group_color, pca1, pca2, ..., pcaN
    df.pca[["row_names"]] <- 1:nrow(df.pca)

    dbDB <- DBI::dbConnect(drv = RMySQL::MySQL(), user = db.user, password = password, host = host, port = port)
    DBI::dbGetQuery(dbDB, paste("CREATE DATABASE IF NOT EXISTS ", prim.data.db, sep=""))
    dbDB <- DBI::dbConnect(drv = RMySQL::MySQL(), user = db.user, password = db.pwd, dbname=prim.data.db, host = host, port = port)
    DBI::dbGetQuery(dbDB, paste("DROP TABLE IF EXISTS ", PCAdbTableName, sep=""))
    dbWriteTable(dbDB, PCAdbTableName, df.pca, row.names= FALSE)

    dbDB <- DBI::dbConnect(drv = RMySQL::MySQL(), user = db.user, password = password, host = host, port = port, dbname=prim.data.db)
    DBI::dbGetQuery(dbDB, paste("ALTER TABLE `",PCAdbTableName,"` ADD UNIQUE(`row_names`)", sep=""))
    DBI::dbGetQuery(dbDB, paste("ALTER TABLE `",PCAdbTableName,"` ADD PRIMARY KEY(`row_names`)", sep=""))
    DBI::dbGetQuery(dbDB, paste0(
        "ALTER TABLE `",
        PCAdbTableName,
        "` CHANGE `sample_id` `sample_id` VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci"
    )
    )

    # Assign sample group and sample group color columns
    var.char.50.cols <- names(df.pca)[grep("sample.group", names(df.pca))]
    if (length(var.char.50.cols) > 0){
        for (i in 1:length(var.char.50.cols)){
            DBI::dbGetQuery(dbDB, paste0(
                "ALTER TABLE `",
                PCAdbTableName,
                "` CHANGE `",
                var.char.50.cols[i],
                "` `",
                var.char.50.cols[i],
                "` VARCHAR(100) CHARACTER SET latin1 COLLATE latin1_swedish_ci"
            )
            )
        }
    }

    dec.6.3.cols <- names(df.pca)[grep("^PC", names(df.pca))]
    for (i in 1:length(dec.6.3.cols)){
        DBI::dbGetQuery(dbDB, paste0("ALTER TABLE ",PCAdbTableName,"
                                CHANGE `",dec.6.3.cols[i],"` `",dec.6.3.cols[i],"` DECIMAL(6,3) NULL DEFAULT NULL"
        )
        )

    }

}

## End of function                                                           ##
###############################################################################

###############################################################################
# (19) upload.datatable.to.database()
###############################################################################

## Indexing of gene name column
# CREATE INDEX idx_mgi_symbol ON p268_rna_seq_table(mgi_symbol)
#' @export

upload.datatable.to.database <- function(
    host = "www.biologic-db.org",
    port = 3306,
    user = "db.user",
    password = "db.pwd",
    prim.data.db = "project.database",
    dbTableName = "rnaseqdbTableName",
    df.data = "df.data.to.upload",
    db.col.parameter.list = list(
        "VARCHAR(255) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("gene_description"),
        "VARCHAR(50) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ENSG", "ENSMUSG", "hgnc_symbol", "mgi_symbol", "uniprot", "entrezgene","display_ptm", "^sequence_window", "p_site_env","for_GSEA_gene_chip","associated_gene_name", "gene_type"),
        "VARCHAR(1) CHARACTER SET latin1 COLLATE latin1_swedish_ci" = c("ppos", "amino_acid", "charge","known_site"),
        "BIGINT(8) NULL DEFAULT NULL" = c("row_names"),
        "INT(8) NULL DEFAULT NULL" = c("row_id", "cluster_order","cluster_id", "count_cut_off", "^position$", "raw_counts"),
        "DECIMAL(6,3) NULL DEFAULT NULL" = c("norm_counts", "NES", "logFC", "lg2_avg", "intensity", "^int", "iBAQ","^localization_prob$"),
        "DECIMAL(6,5) NULL DEFAULT NULL" = c("padj", "pvalue","^pep$")
    ),
    increment = 5000,
    new.table = FALSE,
    first.row.name.index = 1,
    startOnlyWithConnectionCount1 = FALSE,
    cols2Index = NULL
){

    if (host == "clvd1-db-u-p-31" && port == 3306){
        port <- 6008
    }

    if (sum( nchar(names(df.data)) > 64) > 0){
        print("Table names clipped to 64 characters.")
        names(df.data) <- substr(names(df.data), 1, 64)
    }


    if (startOnlyWithConnectionCount1){

        ## helper function ##
        getConnectonCount <- function(
            user= "user",
            password = "password",
            dbname = "prim.data.db",
            host = "host"){
            dbDB <- DBI::dbConnect(
                drv = RMySQL::MySQL(),
                user = user,
                password = password,
                #dbname = prim.data.db,
                host = host
            )

            connectionCount <- as.numeric(
                DBI::dbGetQuery(
                    dbDB,
                    "SELECT COUNT(1) ConnectionCount, user FROM information_schema.processlist WHERE user <> 'system user' AND user = 'boeingS' GROUP BY user;"
                )$ConnectionCount
            )

            DBI::dbDisconnect(dbDB)

            return(connectionCount)
        }

        connectionCount <- getConnectonCount(
            user= user,
            password = password,
            dbname = prim.data.db,
            host = host,
            port = port
        )

        while (connectionCount > 2){
            print(paste(connectionCount, "connections open. Sleep 30 seconds and try again."))

            Sys.sleep(30)


            connectionCount <- getConnectonCount(
                user= user,
                password = password,
                dbname = prim.data.db,
                host = host,
                port = port
            )


        }

    }

    ## Uploading of data frame to database. Happens only if all columns are defined ##
    #library(RMySQL)
    ## Connect to MySQL to check existence of database ##
    dbDB <- DBI::dbConnect(
        drv = RMySQL::MySQL(),
        user = user,
        password = password,
        host = host,
        port = port,
        new.table = TRUE
    )

    ## Create the database if it does not exist already##
    res <- DBI::dbGetQuery(
        dbDB,
        paste(
            "CREATE DATABASE IF NOT EXISTS ",
            prim.data.db,
            sep = ""
        )
    )

    DBI::dbDisconnect(dbDB)

    ## Ensure that df.data has a row_names column ##
    df.data[["row_names"]] <- first.row.name.index:(first.row.name.index+nrow(df.data)-1)

    ## Check if all columns are assigned in db.col.parameter.list ##
    all.col.string.vec <- as.vector(do.call('c', db.col.parameter.list))

    ## Create a vector that contains all col names that contain at least in part the string in all.cols.vec

    ###############################################################################
    ## Function start                                                            ##
    get.all.col.names.with.these.strings <- function(all.col.string.vec){
        all.assigned.cols <- vector(mode="character", length=0)
        for (i in 1:length(all.col.string.vec)){
            pos <- grep(all.col.string.vec[i], names(df.data))
            if (length(pos) > 0){
                all.assigned.cols <- append(all.assigned.cols, names(df.data)[pos])
            }
        }
        return(all.assigned.cols)
    }

    ## End of function                                                           ##
    ###############################################################################
    all.assigned.cols <- get.all.col.names.with.these.strings(all.col.string.vec)

    ## Ensure that all database columns are assigned ##
    not.assigned <- names(df.data)[!(names(df.data) %in% all.assigned.cols)]

    if (length(not.assigned) == 0){
        print("All database columns are defined. Uploading to database...")
    } else {
        print(
            paste0(
                "The following database columns have not been defined: ",
                paste(not.assigned,
                      collapse = ', '
                ),
                ". Datatable not uploaded to database.")
        )
        stop(not.assigned)
    }


    ## Connect to database for dbtable upload  ##
    dbDB <- DBI::dbConnect(
        drv = RMySQL::MySQL(),
        user = user,
        password = password,
        dbname = prim.data.db,
        host = host
    )

    ## Remove all tables with the same name from db ##
    if (new.table){
        res <- DBI::dbGetQuery(
            dbDB,
            paste(
                "DROP TABLE IF EXISTS ",
                dbTableName,
                sep = ""
            )
        )
        DBI::dbDisconnect(dbDB)
    }

    ## Upload up to increment rows in one go ##
    iter <- nrow(df.data)%/%increment
    if (nrow(df.data)%%increment != 0){
        iter <- iter + 1
    }

    for (i in 1:iter){
        if (nrow(df.data) > increment){
            limit <- increment
        } else {
            limit <- nrow(df.data)
        }
        df.temp <- df.data[1:limit,]
        df.data <- df.data[(increment+1):nrow(df.data),]

        uploaded = FALSE
        while (!uploaded){
            tryCatch({
                killDbConnections()
                dbDB <- DBI::dbConnect(
                    drv = RMySQL::MySQL(),
                    user = user,
                    password = password,
                    dbname = prim.data.db,
                    host = host
                )

                ## Upload new dataframe to database ##
                res <- dbWriteTable(
                    dbDB,
                    dbTableName,
                    df.temp,
                    row.names = FALSE,
                    append = TRUE,
                    overwrite = FALSE
                )
                DBI::dbDisconnect(dbDB)
                uploaded = TRUE
                #dbDisconnect(dbDB)
            }, error=function(e){cat("Upload errror :",conditionMessage(e), "\n")})
        }

        print(paste0(i * increment, " rows uploaded to database..."))
        ## Connect to database for dbtable upload  ##
    }

    ####################################################
    ## Function alterDBtable
    alterDBtable <- function(
        cmd.string = "mysql command",
        user = "user",
        password = "password",
        dbname = "prim.data.db",
        host = "host",
        port = 3306
    ){
        if (host == "clvd1-db-u-p-31" && port == 3306){
            port <- 6008
        }
        dbDB <- DBI::dbConnect(
            drv = RMySQL::MySQL(),
            user = user,
            password = password,
            dbname = prim.data.db,
            host = host,
            port = port
        )

        tryCatch({
            DBI::dbGetQuery(
                dbDB,
                cmd.string
            )}, error=function(e) {paste0("Alter not executed. cmd.vector[", i, "]")})

        DBI::dbDisconnect(dbDB)


    }

    ## End of function ##
    ######################
    mysql.cmd = ""
    if (new.table){
        alterDBtable(
            cmd.string = paste(
                "ALTER TABLE `",
                dbTableName,
                "` ADD UNIQUE(`row_names`)",
                sep = ""
            ),
            user = user,
            password = password,
            dbname = dbname,
            host = host,
            port = port
        )

        ## Describe key columns in database table ##
        mysql.cmd <- paste(
            "ALTER TABLE `",
            dbTableName,
            "` ADD UNIQUE(`row_names`)",
            sep = ""
        )


        alterDBtable(
            cmd.string = paste(
                "ALTER TABLE `",
                dbTableName,
                "` ADD PRIMARY KEY(`row_names`)",
                sep = ""
            ),
            user = user,
            password = password,
            dbname = dbname,
            host = host,
            port = port
        )




        ###############################################################################
        ## Characterize and define secondary database columns                        ##
        for (i in 1:length(db.col.parameter.list)) {
            descriptor <- names(db.col.parameter.list[i])
            cols.in.class <-
                get.all.col.names.with.these.strings(db.col.parameter.list[[i]])

            if (length(cols.in.class) > 0) {
                print(
                    paste0(
                        "Assigned ",
                        paste0(cols.in.class, collapse = ', '),
                        " as ",
                        descriptor, "."
                    )
                )

                ## Assign column names to MySQL class ##
                alteration.string <-
                    paste0("ALTER TABLE ", dbTableName, " ")
                for (j in 1:length(cols.in.class)) {
                    alteration.string <- paste0(
                        alteration.string,
                        paste0(
                            "CHANGE `", cols.in.class[j], "` `", cols.in.class[j], "` ", descriptor, ", "
                        )
                    )
                }
                ## Remove last comma from string
                alteration.string <-
                    substr(alteration.string, 1, (nchar(alteration.string) - 2))

                ## Carry out alteration
                ## Connect to database for dbtable upload  ##
                ## Connection is repated to avoid loss of a short lived connection.
                alterDBtable(
                    cmd.string = alteration.string,
                    user = user,
                    password = password,
                    dbname = dbname,
                    host = host,
                    port = port
                )
            }
            #print(alteration.string)
            mysql.cmd <- append(mysql.cmd,
                                alteration.string)

            #DBI::dbGetQuery(dbDB,
            #           alteration.string)


        }

    }
    ## End characterize and define secondary database columns                  ##
    #############################################################################
    ## Add index based on row namems ##

    if (length(cols2Index) > 0){
        for (i in 1:length(cols2Index)){
            print("...indexing...")
            cmd.string <- paste0("CREATE INDEX idx_",cols2Index[i]," ON ",dbTableName," (",cols2Index[i],")")

            dbDB <- DBI::dbConnect(
                drv = RMySQL::MySQL(),
                user = user,
                password = password,
                dbname = prim.data.db,
                host = host,
                port = port
            )

            tryCatch({
                DBI::dbGetQuery(
                    dbDB,
                    cmd.string
                )}, error=function(e) {stop(paste0("Database table not uploaded. Problem adding index ",cols2Index[i],"."))})

            DBI::dbDisconnect(dbDB)

            print(paste0("Datatable ", dbTableName, " successfully uploaded and column(s) ",paste(cols2Index, collapse = " ")," indexed."))
        }
    }
    return(mysql.cmd)
}



#RENAME TABLE p131_rna_seq_table_part_1 TO p131_rna_seq_table
#INSERT INTO p131_rna_seq_table SELECT * FROM p131_rna_seq_table_part_2;
#INSERT INTO p131_rna_seq_table SELECT * FROM p131_rna_seq_table_part_3;
#INSERT INTO p131_rna_seq_table SELECT * FROM p131_rna_seq_table_part_4;
#INSERT INTO p131_rna_seq_table SELECT * FROM p131_rna_seq_table_part_5;

#DROP TABLE p131_rna_seq_table_part_2;
#DROP TABLE p131_rna_seq_table_part_3;
#DROP TABLE p131_rna_seq_table_part_4;
#DROP TABLE p131_rna_seq_table_part_5;

## End of function                                                           ##
###############################################################################

###################################
# (17) Load datatable from database#
###################################
#' @export

import.db.table.from.db <- function(dbtable = "interpro_categori",
                                    dbname = "reference_categories_db_new",
                                    user     = "boeings",
                                    password = "",
                                    host     = "www.biologic-db.org"
){
    ## Helper function ##
    oldw <- getOption("warn")
    options(warn = -1)



    ## End helper function ##


    #library(RMySQL)
    ## Create connection
    dbDB <- DBI::dbConnect(RMySQL::MySQL(), user = user, password = password, host = host, port = port, dbname=dbname)
    ## Get number of expected rows from query ##
    nrows.to.download <- DBI::dbGetQuery(dbDB, paste0("SELECT COUNT(*) FROM ",dbtable))



    DBI::dbDisconnect(dbDB)
    download = TRUE
    i=1
    while (download) {
        dbDB <- DBI::dbConnect(RMySQL::MySQL(), user = user, password = password, host = host, port = port, dbname=dbname)
        out <- tryCatch({
            df.data = DBI::dbGetQuery(dbDB, paste0("SELECT DISTINCT * FROM ", dbtable))
        },
        error=function() {
            message("Database error")
            # Choose a return value in case of error
        }
        )

        DBI::dbDisconnect(dbDB)
        if (nrow(df.data) == nrows.to.download){
            download = FALSE
            print(paste0(nrow(df.data), " of ", nrows.to.download, " rows downloaded."))
        } else {
            print(paste0("Expected: ", nrows.to.download, ". Downloaded: ", nrow(df.data), ". "))
            print(paste0("Download failed for the ", i, " time. Try again..."))
            i = i+1
        }
    }
    options(warn = oldw)
    return(df.data)
}

# End of function

###############################################################################
## (18) Retrieve gene category from db                                       ##
###############################################################################
#' @export

retrieve.gene.category.from.db <- function(
    cat_id         = "mysigdb_c5_MF__127",
    dbname         = "reference_categories_db_new",
    user           = "boeings",
    password       = "",
    host           = "www.biologic-db.org",
    port           = 3306,
    gene.symbol    = "mgi_symbol",
    print.cat.name = TRUE
){
    library(RMySQL)

    if (host == "clvd1-db-u-p-31" && port == 3306){
        port <- 6008
    }

    table <- unlist(
        strsplit(
            cat_id, "__"
        )
    )[1]

    ## Query category name ##
    drv = RMySQL::MySQL()

    sql.query = paste0(
        "SELECT cat_id, cat_name from ",
        table,
        " WHERE cat_id = '",
        cat_id, "'"
    )

    dbDB <- dbConnect(
        drv      = RMySQL::MySQL(),
        user     = user,
        password = password,
        host     = host,
        port = port,
        dbname   = dbname
    )

    cat.vec = dbGetQuery(
        dbDB,
        sql.query
    )

    DBI::dbDisconnect(dbDB)

    if (print.cat.name){
        print(
            paste(
                "Retrieved category: ",
                paste0(
                    cat.vec$cat_name,
                    collapse = ", "
                )
            )
        )

        print(
            paste(
                "Retrieved category ID: ",
                paste0(
                    cat.vec$cat_id,
                    collapse = ", "
                )
            )
        )

    }


    ## Query genes ##
    sql.query = paste0(
        "SELECT ",
        gene.symbol,
        " from ",
        table,
        " WHERE cat_id = '",
        cat_id, "'"
    )

    dbDB <- dbConnect(
        drv      = RMySQL::MySQL(),
        user     = user,
        password = password,
        host     = host,
        port = port,
        dbname   = dbname
    )

    cat.vec = dbGetQuery(
        dbDB,
        sql.query
    )[,gene.symbol]

    DBI::dbDisconnect(dbDB)

    cat.vec <- unlist(
        strsplit(
            cat.vec,
            ";"
        )
    )
    cat.vec <- cat.vec[!is.na(cat.vec)]
    cat.vec <- cat.vec[cat.vec != ""]
    cat.vec = as.vector(unique(cat.vec))
    return(cat.vec)
}

## End of function                                                           ##
###############################################################################


###############################################################################
## (2B) createSettingsFile()                                                 ##
#' @export

createSettingsFile <- function(
    obj = "object",
    df.data = 'database.table',
    defaultXcolName = NULL,
    defaultYcolName = NULL,
    timepointName = NULL,
    sample.order = "names(database.table)[grep('norm_counts', names(databse.table))]", #set to "" to go with default sorting
    heatmapSampleOrder = "lg2_avg vec",
    sample.names = "", # default is sample.order
    count.sample.colors = 'rainbow(length(sample.order))',
    ptm.colum = "display_ptm",
    count.table.headline = "PTM ratio H/L counts for all samples",
    count.table.sidelabel = "Counts",
    venn.slider.selector.strings = 'c("contrast_x_logFC", "constrast_x_padj")',
    plot.selection.strings = 'c(
    "_logFC",
    "_PCA_",
    "_lg10p"
)',
    webSiteDir = "/camp/stp/babs/working/boeings/Stefan/protocol_files/github/biologic/src/experiments",
    upper_heatmap_limit = 3,
    lower_heatmap_limit = -3,
    heamap.headline.text = "heamap.headline.text",
    project_id = "project_id",
    primDataTable = "p123_rna_seq_table",
    pcaDbTable = NULL,
    pointer = "Gene Symbol:"
){

    ###############################################################################
    ## Create timecourse string from dfDesign                                    ##

    createTSparams <- function(
        dfDesign = Obio@dfDesign,
        timepointName = "Timepoint"
    ) {

        tsOrder <- as.numeric(sort(unique(dfDesign$timepoint)))
        scriptVec <- as.vector(NULL, mode = "character")
        scriptVec <- c(
            scriptVec,
            "// New Begin: Timecourse",
            "'timecourse_chart' => array(",
            "    'timepoint_name' => 'Day',",
            "    'display_median' => 'calculate_median',",
            paste0("    'timepoint_array' => array(", paste(tsOrder, collapse = ","),"),"),
            "    'datasets' => array("
        )


        if (length(grep("dataseries_order", names(dfDesign))) > 0){
            if (length(grep("ts_color", names(dfDesign))) > 0){
                dfO <- unique(dfDesign[,c("dataseries", "dataseries_order","ts_color")])
                dfO <- dfO[order(dfO$dataseries_order, decreasing = F),]
                dataseriesVec <- as.vector(dfO$dataseries)
                dataseriesColVec <- as.vector(dfO$ts_color)
            } else {
                dfO <- unique(dfDesign[,c("dataseries", "dataseries_order")])
                dfO <- dfO[order(dfO$dataseries_order, decreasing = F),]
                dataseriesVec <- as.vector(dfO$dataseries)
                dataseriesColVec <- rainbow(length(dataseriesVec))
            }


        } else {
            if (length(grep("ts_color", names(dfDesign))) > 0){
                dfO <- unique(dfDesign[,c("dataseries", "ts_color")])
                dfO <- dfO[order(dfO$dataseries, decreasing = F),]
                dataseriesVec <- as.vector(dfO$dataseries)
                dataseriesColVec <- as.vector(dfO$ts_color)
            } else {
                dataseriesVec <- sort(unique(dfDesign$dataseries))
                dataseriesColVec <- rainbow(length(dataseriesVec))
            }
        }



        for (i in 1:length(dataseriesVec)){
            dfTemp <- unique(
                dfDesign[dfDesign$dataseries == dataseriesVec[i], c("sample.id", "dataseries", "sample.group", "timepoint")]
            )

            dfTemp <- dfTemp[order(dfTemp$timepoint, decreasing = F),]
            timepointVec <- unique(dfTemp$timepoint)
            sampleGroupVec <- unique(dfTemp$sample.group)

            scriptVec <- c(
                scriptVec,
                paste0("'",dataseriesVec[i],"' => array("),
                paste0("    'color' => '",dataseriesColVec[i],"',"),
                paste0("    'sample_group' => array(")
            )

            for (j in 1:length(sampleGroupVec)){
                dfTemp3 <- unique(dfDesign[dfDesign$sample.group %in% sampleGroupVec, c(timepointName, "sample.group")])
                dfTemp3 <- dfTemp3[order(dfTemp3[,timepointName], decreasing = F),]
                timepointVec <- as.numeric(dfTemp3[,timepointName])

                dfTemp2 <- unique(dfTemp[dfTemp$sample.group == sampleGroupVec[j],])
                scriptVec <- c(
                    scriptVec,
                    paste0("'",sampleGroupVec[j],"' => array("),
                    paste0("    'timepoint' =>  ",timepointVec[j],","),
                    paste0("    'sampleDbCols' =>  array("),

                    paste0(
                        sampleCols <- paste0("'norm_counts_", sort(dfTemp2$sample.id), "_TPM'"),
                        collapse = ","
                    ),

                    ")),"
                )
            }
            scriptVec[length(scriptVec)] <- gsub(")),", ")))),",scriptVec[length(scriptVec)])
        }

        scriptVec[length(scriptVec)] <- gsub(",", ")),",scriptVec[length(scriptVec)])

        scriptVec <- c(
            scriptVec,

            "// New End: Timecourse"
        )



        return(scriptVec)

    }

    ## Create timecourse string                                                  ##
    ###############################################################################


    if (sample.order[1] == "" | is.na(sample.order[1])){
        sample.order <- sort(names(database.table)[grep("norm_counts_", names(database.table))])
    }

    if (count.sample.colors[1] == "" | is.na(count.sample.colors[1])){
        count.sample.colors <- rainbow(length(sample.order))
    }

    if (sample.names[1] == "" | is.na(sample.names[1])){
        sample.names <- gsub("norm_counts_", "", sample.order)
        sample.names <- gsub("_", " ", sample.names)
    }

    settingsPhpVec <- c(
        "<?php",
        "",
        "return array(",
        "    'lab' => array(",
        paste0("        'name' => '",obj@parameterList$labname," DB'"),
        "    ),",
        "",
        "    /*",
        "    * Experiment settings",
        "    */",
        paste0("    'data_db_name' => '",obj@dbDetailList$primDataDB,"',"),
        "    'data_db' => array(",
        paste0("            'cat_table_name' => '",obj@parameterList$cat.ref.db.table,"'"),
        "    ),",
        "",
        paste0("    'rnaseq_db_table' => '",obj@parameterList$rnaseqdbTableName,"',"),
        paste0("    'primary_gene_symbol' => '",obj@parameterList$geneIDcolumn,"',"),
        paste0("    'ptm_display_column' => '",obj@parameterList$displayPTMcolumn,"',"),
        "",
        "    'heatmap' => array(",
        paste0("        'upper_limit' => ",upper_heatmap_limit,","),
        paste0("        'lower_limit' => ",lower_heatmap_limit,","),
        paste0("        'headline' => '",obj@parameterList$heamap.headline.text,"',"),
        paste0("        'pointer' => '",pointer,"'"),
        "    ),",
        ""
    )

    ## Add sample array ##
    settingsPhpVec <- c(
        settingsPhpVec,
        "    'samples' => array("
    )

    for (i in 1:length(sample.order)){
        settingsPhpVec <- c(
            settingsPhpVec,
            paste0("        '",sample.order[i],"' => array("),
            paste0("            'color' => '",sample.colors[i],"',"),
            paste0("            'name' => '",sample.names[i],"'"),
            "        )"
        )
        if (i < length(sample.order)){
            settingsPhpVec[length(settingsPhpVec)] <- paste0(
                settingsPhpVec[length(settingsPhpVec)], ","
            )
        }
    }
    settingsPhpVec <- c(
        settingsPhpVec,
        "    ), // End samples array"
    )

    ## Done adding samples ##

    ## Adding barchart parameters ##
    settingsPhpVec <- c(
        settingsPhpVec,
        "    // bar chart",
        "    'count_table' => array(",
        paste0("        'headline' => '", obj@parameterList$count.table.headline,"',"),
        paste0("        'sidelabel' => '", obj@parameterList$count.table.sidelabel,"'"),
        "    ),"
    )
    ## Done adding barchart parameters ##

    ## Adding timecourse parameters ##
    if (!is.null(timepointName)){
        tempVec <- createTSparams(
            dfDesign = Obio@dfDesign,
            timepointName = timepointName
        )

        settingsPhpVec <- c(
            settingsPhpVec,
            tempVec
        )
    }

    ## Done adding timecourse parameters ##


    ## Adding Venn section ##
    if (heatmapSampleOrder[1] == ""){
        heatmapSampleOrder <- names(df.data)[grep("lg2_avg", names(df.data))]
    }

    heatMapString <- paste(heatmapSampleOrder, collapse = "','")
    heatMapString <- paste0("'", heatMapString,"'")

    settingsPhpVec <- c(
        settingsPhpVec,
        "    // Venn Diagram Parameters",
        "    'venn' => array(",
        paste0("        'experiments' => array(", heatMapString,"),"),
        "",
        "    'table' => array(",
        "        'col_name_start' => 11,",
        "        'low_highlight' => -1,",
        "        'high_highlight' => 1",
        "    ),",
        "",
        "    'selection' => array("
    )

    vennCols <- as.vector(NULL, mode = "character")

    ## Make sure all venn cols are numeric ##
    df.data[,vennCols] <- apply(df.data[,vennCols], 2, as.numeric)

    for (i in 1:length(venn.slider.selector.strings)){
        vennCols <- c(
            vennCols,
            names(df.data)[grep(venn.slider.selector.strings[i], names(df.data))]
        )
    }

    for (i in 1:length(vennCols)){
        colMax <- ceiling(max(as.numeric(df.data[,vennCols[i]]), na.rm = TRUE))
        colMin <- floor(min(as.numeric(df.data[,vennCols[i]]), na.rm = TRUE))

        Vname <- vennCols[i]
        Vname <- substr(Vname ,11,100)
        Vname <- gsub("^_", "", Vname)
        Vname <- gsub("_", " ", Vname)

        if (is.numeric(colMax) & is.numeric(colMin)){
            settingsPhpVec <- c(
                settingsPhpVec,
                paste0("        '",vennCols[i],"' => array("),
                paste0("            'name' => '",Vname,"',"),
                paste0("            'slider_min' => ", colMin,","),
                paste0("            'slider_max' => ", colMax,","),
                paste0("            'default_low' => ", colMin,","),
                paste0("            'default_high' => ", colMax,""),
                "        )"
            )
        }

        if (i < length(vennCols)){
            settingsPhpVec[length(settingsPhpVec)] <- paste0(
                settingsPhpVec[length(settingsPhpVec)], ","
            )
        }


    }

    settingsPhpVec <- c(
        settingsPhpVec,
        "    )",  ## Done with venn array
        "    )," ## Done with venn array
    )
    ## Done adding Venn section

    ## Adding scatterplot ##
    scatterCols <- as.vector(NULL, mode = "character")
    for (i in 1:length(plot.selection.strings)){
        scatterCols <- c(
            scatterCols,
            names(df.data)[grep(plot.selection.strings[i], names(df.data))]
        )
    }

    if (!is.null(pcaDbTable)){
        settingsPhpVec <- c(
            settingsPhpVec,
            "    // Scatterplot Parameters",
            paste0("'pca' => '", pcaDbTable, "',")
        )
    }

    if (length(scatterCols) > 0){

        if (is.null(defaultXcolName)){
            defaultXcolName <- scatterCols[1]
        }

        if (is.null(defaultYcolName)){
            defaultXcolName <- scatterCols[2]
        }

        settingsPhpVec <- c(
            settingsPhpVec,
            "    // Scatterplot Parameters",
            "    'scatterplot' => array(",
            paste0("        'default-x' => '",defaultXcolName,"',"),
            paste0("        'default-y' => '",defaultYcolName,"',"),
            "        'selection' => array("

        )

        for (i in 1:length(scatterCols)){
            Sname <- scatterCols[i]
            Sname <- substr(Sname ,11,100)
            Sname <- gsub("^_", "", Sname)
            Sname <- gsub("_", " ", Sname)

            settingsPhpVec <- c(
                settingsPhpVec,
                paste0("            '",scatterCols[i],"' => array("),
                paste0("                'name' => '",Sname,"'"),
                "            )"
            )

            if (i < length(scatterCols)){
                settingsPhpVec[length(settingsPhpVec)] <- paste0(
                    settingsPhpVec[length(settingsPhpVec)], ","
                )
            }


        }

        settingsPhpVec <- c(
            settingsPhpVec,
            "        )", # close scatterplot selection array
            "    )", # close scatterplot  array
            "//End scatterplot" # close scatterplot  array
        )
    }

    ## Done adding scatterplot ##

    ## End of file ##
    settingsPhpVec <- c(
        settingsPhpVec,
        ");"
    )

    ###########################################################################
    ## Create settings.php file                                              ##
    setwd(webSiteDir)
    if (!dir.exists(project_id)){
        dir.create(project_id)
    }


    if (substr(webSiteDir, nchar(webSiteDir), nchar(webSiteDir)) != "/"){
        webSiteDir <- paste0(
            webSiteDir,
            "/"
        )
    }

    FN <- paste0(
        webSiteDir,
        project_id,
        "/settings.php"
    )

    sink(FN)
    for (i in 1:length(settingsPhpVec)){
        cat(settingsPhpVec[i]); cat("\n")
    }
    sink()

    ## Done creating settings.php file                                       ##
    ###########################################################################
}

## End: createSettingsFile()                                                 ##
###############################################################################


###############################################################################
## (1) Datatable.to.website.ptm                                              ##
###############################################################################

# df data input
# Requires a logFC mention in the contrast_X_
#' @export
datatable.to.website.ptm <- function (
    df.data,
    gene.id.column = "ENSMUSG",
    heatmap.genes = "", #Relevant genes has to be the same id class as in gene.id.column
    n.cluster.genes = 6000,
    count.data = FALSE,
    logFC.cut.off = 0, # Either 0 or 1. If 1, then df.data needs to contain
    # a logFC_cut_off column that is either 0 (exclude row in heatmap)
    # or 1 (include row in heatmap)

    selector4heatmap.cols = "logFC",
    heatmap.preprocessing = "lg2.row.avg", # possible: "lg2", "lg2.row.avg", "none"
    hm.cut.off = 4,
    n.hm.cluster = 10,
    count.cut.off.filter = 1
) {
    ###########################################################################
    ## Prepare data table                                                    ##

    # Remove all rows not featuring as an entry in the primary.gene.id.column
    df.data <- df.data[!is.na(df.data[, gene.id.column]), ]
    df.data                 <- unique(df.data)
    df.data[is.na(df.data)] <- ""

    # Enable filtering of low count rows

    if (length(grep("^count_cut_off$", names(df.data))) == 0 ){
        if (count.data){
            df.data[["count_cut_off"]] <- 0
            df.data[,"count_cut_off"]  <- rowSums(df.data[,grep("norm_counts", names(df.data))])
            df.data[,"count_cut_off"]  <- df.data[,"count_cut_off"]/length(grep("norm_counts", names(df.data)))
        } else {
            df.data[["count_cut_off"]] <- 5
        }
    }


    df.data <- df.data[df.data$count_cut_off > count.cut.off.filter,]

    df.data[["row_id"]]     <- paste(
        rep("R", nrow(df.data)),
        1:nrow(df.data),
        sep = ""
    )

    ## Calculate coeficient of variation based on norm_counts column for each row
    df.data["CoVar"] <- 0

    ## Ignore low-intesity rows ##
    df.data[df.data$count_cut_off > 1,"CoVar"] <- apply(
        df.data[df.data$count_cut_off > 1, grep("^norm_counts_", names(df.data))],
        1,
        function(x) sd(x)/mean(x)
    )

    df.data[is.na(df.data)] <- 0
    df.data[df.data$CoVar == Inf, "CoVar"] <- max(df.data[df.data$CoVar < Inf ,"CoVar"])


    # Order from highest to lowest CoVar
    df.data <- df.data[order(df.data$CoVar, decreasing = TRUE),]
    df.data[["CoVarOrder"]] <- 1:nrow(df.data)

    # Select columns for heatmaps and plot display

    df.lg2.row.avg.table <- df.data[, grep(selector4heatmap.cols, names(df.data))]

    row.names(df.lg2.row.avg.table) <- df.data[, "row_id"]

    ## Remove column handle from heatmap column ##
    if (length(grep("contrast_", names(df.lg2.row.avg.table))) > 0){
        names(df.lg2.row.avg.table) <- gsub("contrast_", "", names(df.lg2.row.avg.table))

        ## Remove contrast number ##
        names(df.lg2.row.avg.table)     <- substr(
            names(df.lg2.row.avg.table),
            2,
            100
        )
    } else if (length(grep("norm_counts_", names(df.lg2.row.avg.table))) > 0){
        names(df.lg2.row.avg.table) <- gsub(
            "norm_counts_",
            "",
            names(df.lg2.row.avg.table)
        )
    }

    ## Take care of double digit contrast numbers ##
    names(df.lg2.row.avg.table)  <- gsub(
        "^_",
        "",
        names(df.lg2.row.avg.table)
    )

    names(df.lg2.row.avg.table)     <- paste(
        "lg2_avg_",
        names(df.lg2.row.avg.table),
        sep = ""
    )

    # Ensure numericness
    df.lg2.row.avg.table[, grep("lg2_avg", names(df.lg2.row.avg.table))] <- apply(
        df.lg2.row.avg.table[, grep("lg2_avg", names(df.lg2.row.avg.table))],
        2,
        as.numeric
    )

    ## End df.lg2.row.avg.table creation for all rows                        ##

    ###########################################################################
    ## Create heatmap parameters and default selections                      ##

    # Data preprocessing accoring to selection #
    if (heatmap.preprocessing == "lg2"){
        for (i in 1:nrow(df.lg2.row.avg.table)) {
            df.lg2.row.avg.table[i, ] <- log2(df.lg2.row.avg.table[i,])
        }
    } else if (heatmap.preprocessing == "lg2.row.avg"){
        # Calculate row means
        row_means <- rep(0, nrow(df.lg2.row.avg.table))

        for (i in 1:nrow(df.lg2.row.avg.table)){
            temp.row <- df.lg2.row.avg.table[i, grep("lg2_avg", names(df.lg2.row.avg.table))]
            temp.row <- temp.row[temp.row != 0]
            if (length(temp.row) > 0){
                row_means[i] <- mean(temp.row)
            }
        }

        ## Retired 20160621 ## Start ##
        #row_means <- apply(
        #    df.lg2.row.avg.table[, grep("lg2_avg", names(df.lg2.row.avg.table))], 1, mean
        #)
        ## Retired 20160621 ## End ##

        # Avoid devison by 0
        row_means[row_means == 0] <- 0.001
        for (i in 1:nrow(df.lg2.row.avg.table)) {
            df.lg2.row.avg.table[i, ] <- log2(df.lg2.row.avg.table[i,]/row_means[i])

        }
    }

    # If 'none' or anything else is selected for heatmap processing, The values will be used as 'is' for
    # the heatmap display
    ## Set all Infs to 0 ##
    df.lg2.row.avg.table[df.lg2.row.avg.table == Inf ] <- 0
    df.lg2.row.avg.table[df.lg2.row.avg.table == -Inf ] <- 0


    # Limit top/bottom values of heatmap display
    df.lg2.row.avg.table[df.lg2.row.avg.table > hm.cut.off] <- hm.cut.off
    df.lg2.row.avg.table[df.lg2.row.avg.table < (-1) * hm.cut.off] <- (-1) * hm.cut.off

    df.lg2.row.avg.table[, "row_id"] <- row.names(df.lg2.row.avg.table)
    df.data = merge(df.data, df.lg2.row.avg.table, by.x <- "row_id", by.y = "row_id")

    df.data = na.omit(df.data)
    row.names(df.data) = make.names(df.data[, "row_id"])


    ## Make gene selection for heatmap                                       ##
    # If the selection is provided in the heatmap.genes vector
    # these genes are used
    # Create logFC_cut_off column if not present in dataset
    if (length(grep("logFC_cut_off", names(df.data))) == 0){
        df.data[["logFC_cut_off"]] <- 0
    }

    if (sum(df.data$logFC_cut_off) > 0){
        df.hm.sel <- df.data[df.data$logFC_cut_off == 1,]
    } else {
        df.hm.sel <- df.data
    }


    if (heatmap.genes[1] == "" | is.na(heatmap.genes[1])){
        ## Select gene subset for heatmap based on coefficient of variation
        heatmap.genes <- as.vector(
            unique(
                df.hm.sel[,gene.id.column]
            )
        )
    } else {
        # Make sure all listed gene IDs are present in the dataset
        heatmap.genes <- heatmap.genes[heatmap.genes %in% df.hm.sel[,gene.id.column]]
    }

    # Done selecting heatmap genes #

    # Limiting number of genes for heatmap display accoring to specifications #
    ## Query logFC limitation ##

    # Limit based on Coeficient of variation #
    if (length(heatmap.genes) > n.cluster.genes){
        dfSel <- df.hm.sel
        #row.names(dfSel) <- NULL
        dfSel <- unique(dfSel[,c("CoVar", "CoVarOrder", gene.id.column)])
        dfSel <- dfSel[order(dfSel$CoVarOrder),]
        heatmap.genes <- as.vector(dfSel[,gene.id.column])[1:n.cluster.genes]
    }

    # Create df.cluster #
    if (sum(df.data$logFC_cut_off) > 0){
        df.cluster <- df.data[
            df.data[,gene.id.column] %in% heatmap.genes &
                df.data$logFC_cut_off == 1,
            grep("lg2_avg", names(df.data))
        ]
    } else {
        df.cluster <- df.data[
            df.data[,gene.id.column] %in% heatmap.genes,
            grep("lg2_avg", names(df.data))
        ]
    }


    ## Done selecting genes for heatmap                                          ##
    ###############################################################################

    ###############################################################################
    # Make heatmap function                                                       #
    ###############################################################################

    ## Function definition moved to package ##

    ## Flatening ##
    df.cluster[df.cluster > hm.cut.off] = hm.cut.off
    df.cluster[df.cluster < (-1) * hm.cut.off] = (-1) * hm.cut.off
    m.cluster = data.matrix(df.cluster)
    m.cluster[is.na(m.cluster)] = 0
    colnames(m.cluster) <- gsub("lg2_avg_", "", colnames(m.cluster))
    #m.cluster[m.cluster == 0] = 0.01

    # Make col.sorted heatmap ##
    hm.res = make.hm(
        m.cluster,
        filename = "",
        k.number = n.hm.cluster,
        n.colors = 20,
        hclust.method = "ward.D2",
        dist.method = "euclidean",
        main = "",
        Colv = FALSE
    )

    ## Make col-clustered heatmap ##
    hm.res = make.hm(
        m.cluster,
        filename = "",
        k.number = n.hm.cluster,
        n.colors = 20,
        hclust.method = "ward.D2",
        dist.method = "euclidean",
        main = "",
        Colv = TRUE
    )

    ## Extract cluster order ##
    df.clust.order <- data.frame(hm.res$sorted)
    cluster.ordered.sample.vector <- names(df.clust.order)
    df.clust.order[["cluster_order"]] <- 1:nrow(df.clust.order)
    df.clust.order[, "row_id"] <- row.names(df.clust.order)
    df.clust.order <- df.clust.order[, c("row_id", "cluster_order")]

    ## Extract sample order ##
    df.sample.order <- hm.res$sorted
    sampleColClustOrder <- names(data.frame(df.sample.order))
    ## Re-order lg2_avg and norm_counts accordingly ##
    renameVec <- names(df.data)
    # Remove lg2_avg_entries #
    oldLog2AvgEntries <- names(df.data)[grep("lg2_avg", names(df.data))]
    newLog2AvgEntries <- paste0("lg2_avg_", sampleColClustOrder)

    if (sum(!(oldLog2AvgEntries %in% newLog2AvgEntries)) == 0){
        renameVec <- renameVec[!(renameVec %in% newLog2AvgEntries)]
        renameVec <- c(
            renameVec,
            newLog2AvgEntries
        )
    }

    # Remove norm_counts #
    oldNormCountsEntries <- names(df.data)[grep("norm_counts_", names(df.data))]
    newNormCountsEntries <- paste0("norm_counts_", sampleColClustOrder)

    if (sum(!(oldNormCountsEntries %in% newNormCountsEntries)) == 0){
        renameVec <- renameVec[!(renameVec %in% newNormCountsEntries)]
        renameVec <- c(
            renameVec,
            newNormCountsEntries
        )
    }

    ## Reorder columns in df.data ##
    if (sum(!(names(df.data) %in% renameVec)) == 0){
        df.data <- df.data[,renameVec]
    }

    ## Add to main data table ##
    remove <- as.vector(
        na.omit(
            match(
                df.clust.order[, "row_id"],
                df.data[, "row_id"]
            )
        )
    )

    id.vector = as.vector(df.data[-remove, "row_id"])
    df.rest = data.frame(id.vector, rep(0, length(id.vector)))
    names(df.rest) = names(df.clust.order)
    df.clust.order = rbind(df.clust.order, df.rest)
    df.data = merge(df.data, df.clust.order, by.x = "row_id",
                    by.y = "row_id")
    df.data = df.data[!is.na(df.data[, gene.id.column]), ]
    df.data = unique(df.data)

    # Add cluster id
    df.cluster.id <- data.frame(na.omit(hm.res$clusters))
    df.cluster.id[["row_id"]] <- row.names(df.cluster.id)
    names(df.cluster.id)<- c("cluster_id", "row_id")
    df.cluster.id <- df.cluster.id[grep("R", df.cluster.id$row_id),]
    # Adding all other ids
    row_id <- df.data[!(df.data$row_id %in% df.cluster.id$row_id), "row_id"]
    cluster_id <- rep(0, length(row_id))
    df.add <- rbind(df.cluster.id, data.frame(cluster_id, row_id))

    df.data <- merge(df.data, df.add, by.x = "row_id", by.y="row_id", all=TRUE)
    df.data[is.na(df.data)] = ""

    # # Add gene descripton
    # if (!exists("df.anno")){
    #     df.anno <- read.delim(
    #         gene.id.table,
    #         header = TRUE,
    #         sep = "\t",
    #         stringsAsFactors = FALSE
    #     )
    # }
    #
    # # Remove all entries from df.anno that are not present in df.data
    # df.anno <- df.anno[df.anno[,gene.id.column] %in% df.data[,gene.id.column],]
    # if (!add.uniprot.column){
    #     df.anno$uniprot = NULL
    #     df.anno <- unique(df.anno)
    # }
    #
    # df.anno <- unique(df.anno)
    #
    # df.data <- merge(
    #     df.data,
    #     df.anno,
    #     by.x = gene.id.column,
    #     by.y = gene.id.column,
    #     all=TRUE
    # )
    #
    # df.data[is.na(df.data)] = ""
    # df.data = unique(df.data)
    #
    # if (gene.id.column == "mgi_symbol" | gene.id.column == "ENSMUSG") {
    #     ENSG <- "ENSMUSG"
    # } else if (gene.id.column == "hgnc_symbol" | gene.id.column == "ENSG") {
    #     ENSG <- "ENSG"
    # }
    #
    # if (gene.id.column != "mgi_symbol" | gene.id.column != "hgnc_symbol"){
    #     df.data$gene_description <- paste0(
    #         df.data$gene_description,
    #         " (",
    #         df.data[,gene.id.column],
    #         ")"
    #     )
    # }

    ###########################################################################
    ## Deal with PTM datasets                                                ##
    if (length(grep("p_site_env", names(df.data))) > 0) {
        #Trim if the sequence window is to big
        length <- nchar(df.data$p_site_env)
        center <- ((length -1)/2)
        df.data$p_site_env <- ifelse(
            (length > 15),
            substr(df.data$p_site_env, center-6,center+8),
            df.data$p_site_env
        )

        one = tolower(substr(df.data$p_site_env, 1, 7))
        two = toupper(substr(df.data$p_site_env, 8, 8))
        three = tolower(substr(df.data$p_site_env, 9, 16))
        df.data$p_site_env = paste(one, two, three, sep = "")

        ################################################################################
        #Add ppos columns to datatable
        ################################################################################
        ppos.vec = c("ppos_minus_7","ppos_minus_6","ppos_minus_5","ppos_minus_4","ppos_minus_3","ppos_minus_2","ppos_minus_1","ppos",
                     "ppos_plus_1", "ppos_plus_2","ppos_plus_3","ppos_plus_4","ppos_plus_5","ppos_plus_6","ppos_plus_7")


        #In this dataset not all sequences are associated with an p_site_env
        #df.data[df.data$p_site_env == "", "p_site_env"] = substr(df.data[df.data$p_site_env == "", "sequence_window"], 9,23)

        for (i in 1:length(ppos.vec)){
            df.data[[ppos.vec[i]]] = sapply(
                df.data$p_site_env,
                function(x) substr(x, i,i))
        }

        # Done adding ppos columns
    }
    ## Done dealing with PTM datasets                                        ##
    ###########################################################################

    df.data[is.na(df.data)] = ""
    df.data = unique(df.data)
    df.data = df.data[!is.na(df.data[, gene.id.column]), ]
    df.data[["row_names"]] = 1:nrow(df.data)
    names(df.data) = gsub("[.]", "_", names(df.data))
    names(df.data) = gsub(" ", "_", names(df.data))

    return(df.data)
}
## End of function                                                           ##
###############################################################################


