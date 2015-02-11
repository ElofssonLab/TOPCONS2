
First some background:

Data is (often) trivially parallelizable, but still most bioinformatics pipelines are not developed to use more than rather limited set of computational units. A typical example is the search of a sequence against a database. Here a protein sequence of ~10^3 letters is searched against a database of ~10^7 sequence. Using dynamic programming this would require in the order of 10^13 calculations. A number within reach for modern computational resources. However, if you want to search all against all sequences you end up with 10^20 calculations, i.e not possible. However, already decades ago it was clear that heuristic methods such as BLAST (basically indexing to quickly identify possibly relevant hits) could dramatically speed up these searches without a significant loss in sensitivity. This worked well for many years and a typical blast search only took a few minutes. However, during the last years the database sizes has increased significantly faster than computational speed making the search time slower. As a BLAST search often is a crucial step in any bioinformatics pipeline and it is often not only linear with the number of input sequences this is starting to cause problems.


We are currently developing methods to solve this problem. In short we are developing methods that can do the full search in two steps. First a search is performed against a clustered small database of protein domains and secondly the hits found in this database are incorporated in a new database that is used for the final search. When developing a web-server (using the PDC cloud) we have used this strategy to decrease the response time from a few minutes to an average of 5 seconds, see http://topcons.net/. Similar strategies should be possible to apply to a number of related projects.
 


The best performance to predict membrane protein, as well as many other aspects of proteins, is to use a profile. The best profiles are obtained when searching a large database, such as uniprot. However, given the rapid increase in database sizes the search often can take several minutes given a single computer. This i not optimal for the experience of a web-server, where the user wants a rapid response. In the earlier topcons server this problem was circumvented by using a smaller database consisting only of membrane proteins. This did not significantly affect the predictions of topologies in membrane proteins. However, many non-membrane proteins had very few related proteins and therefore for some of these proteins incorrect membrane regions were predicted. In the new version of topcons we have switched to a two-step pipeline. First the query searches a domain-database and then all full-length sequences from this domain database are used to create a query-specific database that is used for the homology search. Because the domain database and the number of hits found found both are much smaller than a database containing all proteins this search is much faster and as almost all proteins have domain hits the resulting profile is virtually identical to the ones found when searching the entire database. Therefore, this two-step procedure combines both the speed in the earlier version of topcons  using a membrane specific database with the ability to separate membrane and non-membrane proteins obtained when using the entire uniprot"


Current implementation (and problems)

Correct me if I am wrong. Currently
+ We first use pfamscan to search Pfam for domain hits
+ We then extract all full length sequences and create a new database
+ We search this database using psiblast.

Limitations
+ For some proteins we have no Pfam hits, this makes the search for these slow.
+ We do currently use psiblast and could try jackhmmer.
+ We do not know how many homologs we miss using the Pfam search

Proposals
+ Implement everything using hmmer3 (better faster etc).
+ Benchmark better

Benchmarks
* Test on how many (and which) homologs are detected
* Test on the following predictors
 * Psipred
 * SCAMPI
 * NeturfP
 * Psicov
 * plmDVA
 * PconsC3b
 * ????
