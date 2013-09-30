

from SAP.NW import align

#retval, one, two = align(2.3, 3.3)
#retval, one, two = align('ACCG','ACCT')

# retval, one, two = align("GGGATAAACGTTTCGCCGAT",
#                          "ATGTAACGCGATG")
# 
# print retval
# print one
# print two
# 
# retval, one, two = align("TTGGGATAAACGTTTCGCCGAT",
#                          "TTATGTAACGCGATG")
# 
# print retval
# print one
# print two
# 
# 
# retval, one, two = align("AGAGCACAGATTTTGGGATAAACGTTTCGCCGATAGAGCACATATT",
#                          "TTATGTAACGCGATG")
# 
# print retval
# print one
# print two

retval, one, two = align("TAAAATTGGGAAAATTTTAT",
                         "TAATTGGGAAAATTAT")

print retval
print one
print two
