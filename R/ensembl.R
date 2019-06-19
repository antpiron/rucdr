
#' @export
ensembl <- function(user="readonly", passwd="readonly",
                    db="homo_sapiens_core_95_38", host="localhost")
{
    ## print(DBI::dbListTables(storiesDb))
    ## DBI::dbDisconnect(storiesDb)

    env <- new.env()
    env$handle <- DBI::dbConnect(RMariaDB::MariaDB(), user=user,
                                 password=passwd, dbname=db, host=host)

    reg.finalizer(env,
                  function (e) {
                      logging("closing handle", .module="ensembl")
                      DBI::dbDisconnect(e$handle)
                  },
                  onexit=TRUE)
    
    return(structure(list(env=env),
                     class=c("ensembl", "genomicDB")))
}


#' @export
geneDescription <- function (ensembl, ...)
{
    UseMethod("geneDescription", ensembl)
}

#' @export
geneDescription.ensembl <- function(ensembl, genes)
{
    # You can fetch all results:
    res <- DBI::dbSendQuery(ensembl$env$handle,
                            "SELECT gene.gene_id, stable_id AS ensembl_id, value AS short_id, description, biotype, analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, display_xref_id, source, is_current, canonical_transcript_id, version, created_date, modified_date  FROM gene, gene_attrib WHERE gene.gene_id = gene_attrib.gene_id AND gene_attrib.attrib_type_id = 4 AND gene.stable_id  = ?")
    out <- t(sapply(genes, function(gene)
    {
        DBI::dbBind(res, gene)
        DBI::dbFetch(res)
    }))
    DBI::dbClearResult(res)

    return(out)
}

#' @export
transcriptDescription <- function (ensembl, ...)
{
    UseMethod("transcriptDescription", ensembl)
}

#' @export
transcriptDescription.ensembl <- function(ensembl, transcripts)
{
    # You can fetch all results:
    res <- DBI::dbSendQuery(ensembl$env$handle,
                            "SELECT *  FROM transcript, transcript_attrib WHERE transcript.transcript_id = transcript_attrib.transcript_id AND transcript_attrib.attrib_type_id = 4 AND transcript.stable_id  = ?")
    out <- t(sapply(transcripts, function(transcript)
    {
        DBI::dbBind(res, transcript)
        DBI::dbFetch(res)
    }))
    DBI::dbClearResult(res)

    return(out)
}

