import os
import pickle
import re
import sys, tempfile, httplib

from Bio import Entrez
from SAP.Exceptions import AnalysisTerminated
from SAP.ProgressBar import ProgressBar
import Taxonomy

def safe_generator(gen):
    """
    Wrapper generator for another generator
    to allow exception handling

    :param gen:
    :return:
    """
    while True:
        try:
            next = gen.next()
        except (StopIteration, httplib.IncompleteRead):
            raise
        except Exception as e:
          yield (None, e)
        else:
          yield (next, None)


def chunks(l, n):
    """
    Generator for returning n next elements from a list

    :param l:
    :param n:
    :return:
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def same_species(x, y):
    """
    Test if two species names represent the same species.
    (if the two first words are the same)

    :param x:
    :param y:
    :return:
    """
    a, b = x.split(), y.split()
    if a[0] == b[0] and a[1] == b[1]:
        return ' '.join(a[:2])
    else:
        return None


def retrieve_sequence_records(query_key, webenv, count, temp):
    """
    Retrieves and parses genbank XML records, writes them
    to a temp file, and creates a mapping from taxid to gi.
    :param query_key:
    :param webenv:
    :param count:
    :param temp:
    :return:
    """
    pg = ProgressBar(maxValue=count)
    pg(0)
    i = 1
    batch_size = 500
    not_downloaded = list()
    taxid2gi = dict()
    for start in range(0,count,batch_size):
        end = min(count, start+batch_size)
        fetch_handle = Entrez.efetch(db="nucleotide", retmode="xml",
            retstart=start, retmax=batch_size,
            webenv=webenv, query_key=query_key)
        records = Entrez.parse(fetch_handle, validate=False)
        for entry, exception in safe_generator(records):
            pg(i)
            i += 1

            if exception:
                not_downloaded.append((None, str(exception)))
                continue

            taxid = None
            gi = None
            seq = None
            try:
                for qual_dict in entry['GBSeq_feature-table'][0]['GBFeature_quals']:
                    if qual_dict['GBQualifier_value'].startswith('taxon:'):
                        taxid = qual_dict['GBQualifier_value'][6:]
                        break
                for x in entry['GBSeq_other-seqids']:
                    if x.startswith('gi|'):
                        gi = x[3:]
                seq = entry['GBSeq_sequence']
                if not (gi and taxid and seq):
                    not_downloaded.append((gi, 'error'))
                    continue
                taxid2gi.setdefault(taxid, []).append(gi)
                print >>temp, ">%s\n%s" % (gi, seq)
            except KeyError:
                not_downloaded.append((gi, 'xml error'))

    return taxid2gi, not_downloaded


def retrieve_taxonomies(taxid2gi):
    """
    Retrieves and parses XML taxonomy entries from genbank
    and to get a mapping from gi to taxonomy.
    """
    tax_not_retrieved = list()
    pg = ProgressBar(maxValue=len(taxid2gi))
    pg(0)
    gi2taxonomy = dict()
    batch_size = 500
    i = 1
    sep_regex = re.compile('[:;,]')
    for taxids in chunks(taxid2gi.keys(), batch_size):
        fetch_handle = Entrez.efetch(db="taxonomy", id=taxids, rettype="xml", retmax=batch_size)
        records = Entrez.parse(fetch_handle, validate=False)

        for entry, exception in safe_generator(records):
            pg(i)
            i += 1

            if exception:
                tax_not_retrieved.append((None, exception))
                continue

            taxid = entry['TaxId']

            taxonomy = Taxonomy.Taxonomy()

            # make sure we have a set of unique ranks:
            unique_ranks = dict()
            for x in entry['LineageEx']:
                rank, name = x['Rank'], x['ScientificName']
                if rank == 'no rank':
                    continue
                if rank not in unique_ranks:
                    unique_ranks[rank] = name
                elif rank == 'species' and same_species(name, unique_ranks[rank]):
                    unique_ranks[rank] = same_species(name, unique_ranks[rank])
                    # We do this because both species and subspecies (with three names) show up as 'species' in rank
                else:
                    tax_not_retrieved.append((taxid, "Error: duplicate ranks: %s %s %s" % (rank, name, unique_ranks[rank])))
                    continue

            # check for seperator characters in ranks and taxons
            if sep_regex.search(''.join(unique_ranks.keys() + unique_ranks.values())):
                tax_not_retrieved.append((taxid, "Interspersed seperator chars %s" % (str(unique_ranks))))
                continue

            # create taxonomy
            for rank, name in unique_ranks.items():
                taxonomy.add(Taxonomy.TaxonomyLevel(name, rank))
            taxonomy.organism = entry['ScientificName']

            # if not taxonomy.name('species'):
            #     species = ' '.join(taxonomy.organism.split()[:2])
            #     taxonomy.add(Taxonomy.TaxonomyLevel(species, 'species'))

            try:
                for gi in taxid2gi[taxid]:
                    gi2taxonomy[gi] = taxonomy
            except KeyError:
                tax_not_retrieved.append((taxid, "retrieved taxid not found in genbank record"))

    return gi2taxonomy, tax_not_retrieved


def write_database(gi2taxonomy, temp, output_file_name):
    """
    Writes the fasta database by mapping the genbank
    fasta to the taxonomy to the taxonomy information
    """
    included_gis = set()
    output_file = open(output_file_name, 'w')
    excluded = list()
    print >>sys.stderr, '\nWriting database'
    from Bio import SeqIO
    for record in SeqIO.parse(temp, "fasta") :
        try:
            taxonomy = gi2taxonomy[record.id]
            if not all((record.id, str(taxonomy), taxonomy.organism)):
                excluded.append((record.id, 'header fields missing'))
                continue
            if record.id in included_gis:
                excluded.append((record.id, 'dublicate gi'))
                continue
            included_gis.add(record.id)
            header = record.id + ' ; ' + str(taxonomy) + ' ; ' + taxonomy.organism
            seq = str(record.seq)
            print >>output_file, ">%s\n%s\n" % (header, seq)
        except KeyError:
            excluded.append((record.id, 'taxonomy not found for gi'))
    output_file.close()
    temp.close()

    return excluded

def compileDatabase(query, email, output_file_name):
    """
    Compiles a fasta database formatted with taxonomy information
    for use with SAP by retrieving sequence and taxonomy information
    defined by an entrez query

    :param query:
    :param email:
    :param output_file_name:
    :return:
    """

    print "Query:", query

    try:
        # Make the search
        Entrez.email = email
        search_handle = Entrez.esearch(db="nucleotide",term=query, usehistory="y", retmax=1000000)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        # Check if you got all results
        gi_list = search_results["IdList"]
        count = int(search_results["Count"])
        assert count == len(gi_list), (count, len(gi_list))

        # entrez history (session_cookie and query_key )
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        answer = raw_input('total nr of entries for download: %d. Proceed? yes/no: ' % (count))
        if answer != 'yes':
            sys.exit()

        print >>sys.stderr, "Downloading"
        temp = tempfile.TemporaryFile(mode='w+t')
        taxid2gi, not_downloaded = retrieve_sequence_records(query_key, webenv, count, temp)

        print >>sys.stderr, '\nRetrieving taxonomies'
        gi2taxonomy, tax_not_retrieved = retrieve_taxonomies(taxid2gi)

    except httplib.IncompleteRead:
        raise AnalysisTerminated(None, "Connection to NCBI broke - try again.")

    temp.seek(0) # rewind temp for reading
    excluded = write_database(gi2taxonomy, temp, output_file_name)

    with open(os.path.splitext(output_file_name)[0] + '.query', 'w') as f:
        print >>f, query

    print "Excluded %d of %d records:" % (len(not_downloaded) + len(tax_not_retrieved) + len(excluded), count)
    print "\tSequence download failed: %d" % (len(not_downloaded))
    print "\ttaxonomy retrieval failed: %d" % (len(tax_not_retrieved))
    print "\tTaxonomy not mapped to sequence entry: %d" % (len(excluded))
    excludedFileName = os.path.splitext(output_file_name)[0] + '.excluded'
    print "See", excludedFileName
    with open(excludedFileName, 'w') as f:
        for t in not_downloaded + tax_not_retrieved + excluded:
            print >>f, "\t".join(t)
