
import sys, tempfile, httplib

from Bio import Entrez
from SAP.ProgressBar import ProgressBar
import Taxonomy

def safe_generator(gen):
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
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def same_species(x, y):
    a, b = x.split(), y.split()
    if a[0] == b[0] and a[1] == b[1]:
        return ' '.join(a[:2])
    else:
        return None


def retrieve_sequence_records(query_key, webenv, count, temp):

    pg = ProgressBar(maxValue=count)
    pg(0)
    i = 1
    batch_size = 100
    nr_not_downloaded = 0
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
                print exception
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
                    print 'error', gi
                    nr_not_downloaded += 1
                    continue
                taxid2gi.setdefault(taxid, []).append(gi)
                print >>temp, ">%s\n%s" % (gi, seq)
            except KeyError:
                print 'xml error'
                nr_not_downloaded += 1

    return taxid2gi, nr_not_downloaded


def retrieve_taxonomies(taxid2gi):

    nr_tax_not_retrieved = 0
    pg = ProgressBar(maxValue=len(taxid2gi))
    pg(0)
    gi2taxonomy = dict()
    for taxids in chunks(taxid2gi.keys(), 100):
        fetch_handle = Entrez.efetch(db="taxonomy", id=taxids, rettype="xml", retmax=100)
        records = Entrez.parse(fetch_handle, validate=False)
        i = 1
        for entry, exception in safe_generator(records):
            pg(i)
            i += 1
            taxid = entry['TaxId']

            if exception:
                print exception
                continue

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
                    print "Error: duplicate ranks:", rank, name, unique_ranks[rank]
                    nr_tax_not_retrieved += 1
                    continue

            # create taxonomy
            for rank, name in unique_ranks.items():
                taxonomy.add(Taxonomy.TaxonomyLevel(name, rank))
            taxonomy.organism = entry['ScientificName']

            if not taxonomy.name('species'):
                species = ' '.join(taxonomy.organism.split()[:2])
                taxonomy.add(Taxonomy.TaxonomyLevel(species, 'species'))

            for gi in taxid2gi[taxid]:
                gi2taxonomy[gi] = taxonomy

    return gi2taxonomy, nr_tax_not_retrieved


def write_database(gi2taxonomy, temp, output_file_name):

    output_file = open(output_file_name, 'w')
    nr_excluded = 0
    print >>sys.stderr, '\nWriting database'
    from Bio import SeqIO
    for record in SeqIO.parse(temp, "fasta") :
        try:
            taxonomy = gi2taxonomy[record.id]
            fields = (record.id, str(taxonomy), taxonomy.organism)
            assert all(fields), fields
            header = record.id + ' ; ' + str(taxonomy) + ' ; ' + taxonomy.organism
            seq = str(record.seq)
            print >>output_file, ">%s\n%s\n" % (header, seq)
        except KeyError:
    #		print "Error", record.id
            nr_excluded += 1
    output_file.close()
    temp.close()

def compileDatabase(query, email, output_file_name):

    print "Query:", query

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
    taxid2gi, nr_not_downloaded = retrieve_sequence_records(query_key, webenv, count, temp)

    print >>sys.stderr, '\nRetrieving taxonomies'
    gi2taxonomy, nr_tax_not_retrieved = retrieve_taxonomies(taxid2gi)

    temp.seek(0) # rewind temp for reading
    write_database(gi2taxonomy, temp, output_file_name)

    print "Excluded %d of %d records" % (nr_not_downloaded + nr_tax_not_retrieved, count)
    print "(not downloaded: %d, no taxonomy %d)" % (nr_not_downloaded, nr_tax_not_retrieved)






