suppressPackageStartupMessages({library(DBI); library(RSQLite)})
con <- dbConnect(SQLite(), "labels/annotations.sqlite")
print(dbGetQuery(con, "SELECT batch_id, COUNT(*) AS n, MIN(batch_score) AS min_s, MAX(batch_score) AS max_s
                       FROM candidates WHERE batch_id IS NOT NULL GROUP BY batch_id"))
print(dbGetQuery(con, "SELECT batch_id, COUNT(*) AS n
                       FROM candidates WHERE batch_id IS NOT NULL AND queue_priority = 0 GROUP BY batch_id"))
print(dbGetQuery(con, "SELECT candidate_id, context, peak, gene, batch_id, ROUND(batch_score, 4) AS s
                       FROM candidates WHERE batch_id = 'review_v1_sync_repress'
                       ORDER BY batch_score DESC LIMIT 5"))
print(dbGetQuery(con, "SELECT candidate_id, context, peak, gene, batch_id, ROUND(batch_score, 4) AS s
                       FROM candidates WHERE batch_id = 'review_v1_low_confidence'
                       ORDER BY batch_score ASC LIMIT 5"))
dbDisconnect(con)
