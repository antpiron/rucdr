
#' @export
ensembl <- function(user="readonly", passwd="readonly",
                    db="homo_sapiens_core_95_38", host="localhost")
{
    handle <- DBI::dbConnect(RMariaDB::MariaDB(), user=user,
                             password=passwd, dbname=db, host=host)
    ## print(DBI::dbListTables(storiesDb))
    ## DBI::dbDisconnect(storiesDb)

    return(structure(list(handle=handle),
                     class=c("ensembl", "genomicDB")))
}


geneDescription <- function (ensembl, ...)
{
    UseMethod("geneDescription", ensembl)
}

geneDescription.ensembl <- function(ensembl, genes)
{
    # You can fetch all results:
    res <- DBI::dbSendQuery(ensembl$handle,
                            "SELECT gene.gene_id, stable_id AS ensembl_id, value AS short_id, description, biotype, analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, display_xref_id, source, is_current, canonical_transcript_id, version, created_date, modified_date  FROM gene, gene_attrib WHERE gene.gene_id = gene_attrib.gene_id AND gene_attrib.attrib_type_id = 4 AND gene.stable_id  = ?")
    out <- t(sapply(genes, function(gene)
    {
        DBI::dbBind(res, gene)
        DBI::dbFetch(res)
    }))
    DBI::dbClearResult(res)

    return(out)
}

