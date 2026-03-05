import primer3

def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def reduce_to_desired_MT(seq, targetMT = 60, max_len = 36, min_len = 20):
    l = max_len  #length
    while primer3.calc_tm(seq[:l]) > targetMT:
        l -= 1
#         print (l,seq[:l], primer3.calcTm(seq[:l]) )
        if l == min_len:
            break
    l += 1
    return seq[:l], l,  primer3.calc_tm(seq[:l])

def annealing_region (seq, targetMT = 59):
    #reduce to desired melting temp
    p, l , tm= reduce_to_desired_MT(seq, targetMT = targetMT)
    #check for hairpin
    analysis = primer3.calc_hairpin(reduce_to_desired_MT(p)[0])
#     if analysis.tm > targetMT-15:
#         print ("warning, hairpin with Tm = %s" % analysis.tm)
    return p,l, tm, analysis.tm

def primers_for_CDS(CDS,
    #gibson overhangs (pEAQ)
    fwd_oh = 'TTTCGAACTTGGAGAAAGATTGTTAAGCTTCTGTATATTCTGCCCAAATTCGCGACCGGT'.lower(),
    rev_oh = 'CGAATAACAGTAAATTCAAACTAAAGAAAATTTAATGAAACCAGAGTTAAAGGCCTCGAG'.lower()):

    p,l,tm,hp = annealing_region(CDS, targetMT=59)
    fwd = fwd_oh[-(57-l):].lower() + p.upper()
    p,l,tm,hp = annealing_region(rev_comp(CDS), targetMT=59)
    rev = rev_oh[-(57-l):].lower() + p.upper()
    return fwd, rev

def main():
    CDS = input("CDS: ")
    fwd, rev = primers_for_CDS(CDS)
    print("Fwd: ", fwd)
    print("Rev: ", rev)

if __name__ == "__main__":
    main()


#Adapted from C. McClune (Sattely Lab)
