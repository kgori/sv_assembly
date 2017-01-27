import skbio
from skbio import DNA
from skbio.alignment import StripedSmithWaterman as SSW
from skbio.alignment import local_pairwise_align_ssw as ssw
import numpy as np
import os
from collections import defaultdict, Iterator
import pysam


class BreakpointInfo(object):
    __slots__ = ['lower_chromosome', 'lower_position', 'lower_ambiguity',
                 'upper_chromosome', 'upper_position', 'upper_ambiguity',
                 'overlap_sequence', 'untemplated_sequence']

    def __init__(self, **kwargs):
        for slot in self.__slots__:
            setattr(self, slot, kwargs.get(slot))

    def __str__(self):
        if self.overlap_sequence is not None:
            return '{}({})--{}--{}({})'.format(self.lower_ambiguity, _number_range(self.lower_ambiguity, self.lower_position),
                                               self.overlap_sequence,
                                               self.upper_position, _number_range(self.upper_position, self.upper_ambiguity))

        elif self.untemplated_sequence is not None:
            return '{}] {} [{}'.format(self.lower_position,
                                       self.untemplated_sequence,
                                       self.upper_position)

        else:
            return '{}][{}'.format(self.lower_position, self.upper_position)


class LocationInfo(object):
    __slots__ = ['lower_name', 'lower_pos', 'upper_name', 'upper_pos']

    def __init__(self, **kwargs):
        for slot in self.__slots__:
            setattr(self, slot, kwargs.get(slot))

    def chromosome_names(self):
        return self.lower_name, self.upper_name

    def is_close(self, read, limit=10000):
        """
        Return True if the read is near (within location +- limit)
        one of the locations in this LocationInfo
        """
        for slot in self.__slots__:
            if getattr(self, slot) is None:
                raise ValueError("Not fully populated")

        lower_match = (read.reference_name == self.lower_name and
                       self.lower_pos - limit < read.pos < self.lower_pos + limit)

        upper_match = (read.reference_name == self.upper_name and
                       self.upper_pos - limit < read.pos < self.upper_pos + limit)

        return lower_match, upper_match


class SVPlacer(object):
    """
    Used after SV calling, Tigra assembly and BWA alignment,
    SVPlacer predicts SV location to base precision, and 
    describes SV status, re. microhomology, clean-break, or
    non-template sequence inclusion.
    """

    def __init__(self, bamfile, contigfile, suppfile=None):
        """
        Required - bamfile: BAM file of contigs aligned to reference
                 - contigfile: FASTA file of complete contig sequences
        Optional - suppfile: BAM file of extracted supplementary reads
                   (used to generate list of potential breakpoints)
        """
        filecheck(bamfile)
        filecheck(contigfile)
        self.bamfile = bamfile
        self.contigfile = contigfile
        
        if suppfile is not None:
            filecheck(suppfile)
            with pysam.Samfile(suppfile, 'rb') as open_suppfile:
                supp_batches = _batchify(open_suppfile)
                self.potential_breaks = _get_breaks(supp_batches)
        else:
            self.potential_breaks = None

        reads = _extract_reads(bamfile)
        contigs = { contig.metadata['id']: contig for contig in _read_fasta(contigfile) }
        for k in reads:
            for contig_info in reads[k]:
                contig = contigs[contig_info.reads[0].qname]
                contig_info.contig = contig

        self.reads = reads

    def breakpoints(self, crossref=True, location_match=True, **kwargs):
        # Can only cross reference if there is something to check against
        if self.potential_breaks is None:
            crossref = False

        allbkpts = defaultdict(list)
        for k in self.reads:
            for info in self.reads[k]:
                bkpts = [info.breakpoint(*pair) for pair in info.compatible_reads()]
                if crossref:  # only accept if bkpt cross references with split reads
                    bkpts = [bk for bk in bkpts if _crossref_fn(bk, self.potential_breaks)]
                if location_match:  # only accept if bkpt matches suspected location
                    bkpts = [bk for bk in bkpts if _location_match_fn(info, pair, **kwargs)]

                allbkpts[k].extend(bkpts)
        return allbkpts
        

class ContigInfo(object):
    __slots__ = ['reads', 'location', 'contig']

    def __init__(self, **kwargs):
        for slot in self.__slots__:
            setattr(self, slot, kwargs.get(slot))
            if self.reads is None:
                self.reads = []

    def __str__(self):
        return _str_info(self)

    def add_read(self, read):
        """
        Add read to ContigInfo such that self.reads is in genome order
        """
        if self.reads is None:
            self.reads = [read]

        else:
            self.reads.append(read)
            self.reads.sort(key=lambda r: (r.reference_id, r.pos))


    def contig_index(self):
        if None not in (self.reads, self.contig):
            # return [_contig_index(read, self.contig) for read in self.reads]
            return [_cigar_to_index(read) for read in self.reads]
        else:
            raise ValueError("ContigInfo is not fully populated")

    def genome_index(self):
        if None not in (self.reads, self.contig):
            return [_genome_index(read) for read in self.reads]
        else:
            raise ValueError("ContigInfo is not fully populated")

    def full_genome_index(self):
        if None not in (self.reads, self.contig):
            res = []
            for read in self.reads:
                chrom = read.reference_id
                pos = _genome_index(read)
                res.append((chrom, pos[0], pos[1]))
            return res
        else:
            raise ValueError("ContigInfo is not fully populated")

    def reads_overlap(self, read1idx, read2idx):
        """
        Return True if reads at readidx1 and readidx2 overlap on the genome
        """
        fgi = self.full_genome_index()
        full_index_1 = fgi[read1idx]
        full_index_2 = fgi[read2idx]
        return full_index_1[0] == full_index_2[0] and _range_overlap(full_index_1[1:], full_index_2[1:])

    def compatible_reads(self):
        """
        Return indices of "compatible" read pairs, i.e. the read pairs that may indicate a breakpoint.
        These are adjacent on the contig.
        These pairs can then be used as inputs to the .breakpoint method.
        """
        import itertools, math
        pairs = []
        contig_order = _sort_indices(self.contig_index())
        for (read1idx, read2idx) in itertools.combinations(range(len(self.reads)), 2):
            # Compatible pairs are adjacent on the contig
            # overlap = self.reads_overlap(read1idx, read2idx) and nonoverlapping  # don't care about non-overlapping
            adjacent = int(math.fabs(contig_order[read1idx] - contig_order[read2idx])) == 1
            if adjacent:
                pairs.append((read1idx, read2idx))
        return pairs

    def breakpoint(self, read1idx=0, read2idx=1):
        for readidx in (read1idx, read2idx):
            assert 0 <= readidx < len(self.reads)

        if read1idx > read2idx:
            read1idx, read2idx = read2idx, read1idx
        
        # Make sure this ContigInfo is not missing any data
        if None not in (self.reads, self.contig):

            # Determining the end of the breakpoint means figuring out which
            # parts of the reference genome are joined together.
            # There are 8 possibilities:
            #     GENOME_LOWER       GENOME_UPPER
            # 1) ============== | | ==============
            #      a------->b   | |   b------->c  
            #  
            # 2) ============== | | ==============
            #      a------->b   | |   c<-------b  
            #  
            # 3) ============== | | ==============
            #      b<-------a   | |   b------->c  
            #  
            # 4) ============== | | ==============
            #      b<-------a   | |   c<-------b  
            #  
            # 5) ============== | | ==============
            #      b------->c   | |   a------->b  
            #  
            # 6) ============== | | ==============
            #      c<-------b   | |   a------->b  
            #  
            # 7) ============== | | ==============
            #      b------->c   | |   b<-------a  
            #  
            # 8) ============== | | ==============
            #      c<-------b   | |   b<-------a  
            #
            # N.B. Scenarios 5-8 are the reverse strand versions of 1-4.
            # We only want to report as if 1-4 is observed.
            
            contig_indices = self.contig_index()
            genome_indices = self.genome_index()
            lower_contig_idx = contig_indices[read1idx]
            upper_contig_idx = contig_indices[read2idx]
            lower_genome_idx = genome_indices[read1idx]
            upper_genome_idx = genome_indices[read2idx]
            lower_read = self.reads[read1idx]
            upper_read = self.reads[read2idx]
            lower_genome_name = lower_read.reference_name
            upper_genome_name = upper_read.reference_name
            
            if lower_contig_idx < upper_contig_idx:  # "lower/upper" w.r.t genome, index w.r.t. contig
                # This is a "1-4" scenario
                left_read, left_index = lower_read, lower_contig_idx
                right_read, right_index = upper_read, upper_contig_idx
                overlap, overlaplen, untemplated = _break_type(left_read, right_read, left_index, right_index,
                                                               self.contig, reverse=False)

                if (not left_read.is_reverse):
                    if (not right_read.is_reverse):
                        # Scenario 1
                        #    ============== | | ==============
                        #      a------->b   | |   b------->c  

                        breakpt = BreakpointInfo(lower_chromosome=lower_genome_name,
                                                 lower_position=lower_genome_idx[1],
                                                 lower_ambiguity=lower_genome_idx[1] - overlaplen,
                                                 upper_chromosome=upper_genome_name,
                                                 upper_position=upper_genome_idx[0],
                                                 upper_ambiguity=upper_genome_idx[0] + overlaplen,
                                                 overlap_sequence=overlap,
                                                 untemplated_sequence=untemplated)

                        # breakpt = (lower_genome_idx[1] - overlaplen, lower_genome_idx[1],
                        #            upper_genome_idx[0], upper_genome_idx[0] + overlaplen)

                    else:
                        # Scenario 2
                        #    ============== | | ==============
                        #      a------->b   | |   c<-------b  
                        breakpt = BreakpointInfo(lower_chromosome=lower_genome_name,
                                                 lower_position=lower_genome_idx[1],
                                                 lower_ambiguity=lower_genome_idx[1] - overlaplen,
                                                 upper_chromosome=upper_genome_name,
                                                 upper_position=upper_genome_idx[1],
                                                 upper_ambiguity=upper_genome_idx[1] - overlaplen,
                                                 overlap_sequence=overlap,
                                                 untemplated_sequence=untemplated)
                        
                        # breakpt = (lower_genome_idx[1] - overlaplen, lower_genome_idx[1],
                        #            upper_genome_idx[1], upper_genome_idx[1] - overlaplen)

                else:
                    if (not right_read.is_reverse):
                        # Scenario 3
                        #    ============== | | ==============
                        #      b<-------a   | |   b------->c  
                        breakpt = BreakpointInfo(lower_chromosome=lower_genome_name,
                                                 lower_position=lower_genome_idx[0],
                                                 lower_ambiguity=lower_genome_idx[0] + overlaplen,
                                                 upper_chromosome=upper_genome_name,
                                                 upper_position=upper_genome_idx[0],
                                                 upper_ambiguity=upper_genome_idx[0] + overlaplen,
                                                 overlap_sequence=overlap,
                                                 untemplated_sequence=untemplated)

                        # breakpt = (lower_genome_idx[0] + overlaplen, lower_genome_idx[0],
                        #            upper_genome_idx[0], upper_genome_idx[0] + overlaplen)

                    else:
                        # Scenario 4
                        #    ============== | | ==============
                        #      b<-------a   | |   c<-------b  
                        breakpt = BreakpointInfo(lower_chromosome=lower_genome_name,
                                                 lower_position=lower_genome_idx[0],
                                                 lower_ambiguity=lower_genome_idx[0] + overlaplen,
                                                 upper_chromosome=upper_genome_name,
                                                 upper_position=upper_genome_idx[1],
                                                 upper_ambiguity=upper_genome_idx[1] - overlaplen,
                                                 overlap_sequence=overlap,
                                                 untemplated_sequence=untemplated)

                        # breakpt = (lower_genome_idx[0] + overlaplen, lower_genome_idx[0],
                        #            upper_genome_idx[1], upper_genome_idx[1] - overlaplen)

            
            else:
                # "5-8" scenario
                left_read, left_index = upper_read, upper_contig_idx
                right_read, right_index = lower_read, lower_contig_idx
                overlap, overlaplen, untemplated = _break_type(left_read, right_read, left_index, right_index, 
                                                               self.contig, reverse=True)

                if (not left_read.is_reverse):
                    if (not right_read.is_reverse):
                        # Scenario 5
                        #    ============== | | ==============
                        #      b------->c   | |   a------->b
                        breakpt = BreakpointInfo(lower_chromosome=lower_genome_name,
                                                 lower_position=lower_genome_idx[0],
                                                 lower_ambiguity=lower_genome_idx[0] + overlaplen,
                                                 upper_chromosome=upper_genome_name,
                                                 upper_position=upper_genome_idx[1],
                                                 upper_ambiguity=upper_genome_idx[1] - overlaplen,
                                                 overlap_sequence=overlap,
                                                 untemplated_sequence=untemplated)

                        # breakpt = (upper_genome_idx[1] - overlaplen, upper_genome_idx[1],
                        #            lower_genome_idx[0], lower_genome_idx[0] + overlaplen)

                    else:
                        # Scenario 6
                        #    ============== | | ==============
                        #      c<-------b   | |   a------->b  
                        breakpt = BreakpointInfo(lower_chromosome=lower_genome_name,
                                                 lower_position=lower_genome_idx[1],
                                                 lower_ambiguity=lower_genome_idx[1] - overlaplen,
                                                 upper_chromosome=upper_genome_name,
                                                 upper_position=upper_genome_idx[1],
                                                 upper_ambiguity=upper_genome_idx[1] - overlaplen,
                                                 overlap_sequence=overlap,
                                                 untemplated_sequence=untemplated)

                        # breakpt = (upper_genome_idx[1] - overlaplen, upper_genome_idx[1],
                        #            lower_genome_idx[1], lower_genome_idx[1] - overlaplen)

                else:
                    if (not right_read.is_reverse):
                        # Scenario 7
                        #    ============== | | ==============
                        #      b------->c   | |   b<-------a  
                        breakpt = BreakpointInfo(lower_chromosome=lower_genome_name,
                                                 lower_position=lower_genome_idx[0],
                                                 lower_ambiguity=lower_genome_idx[0] + overlaplen,
                                                 upper_chromosome=upper_genome_name,
                                                 upper_position=upper_genome_idx[0],
                                                 upper_ambiguity=upper_genome_idx[0] + overlaplen,
                                                 overlap_sequence=overlap,
                                                 untemplated_sequence=untemplated)

                        # breakpt = (upper_genome_idx[0] + overlaplen, upper_genome_idx[0],
                        #            lower_genome_idx[0], lower_genome_idx[0] + overlaplen)

                    else:
                        # Scenario 8
                        #    ============== | | ==============
                        #      c<-------b   | |   b<-------a  
                        breakpt = BreakpointInfo(lower_chromosome=lower_genome_name,
                                                 lower_position=lower_genome_idx[1],
                                                 lower_ambiguity=lower_genome_idx[1] - overlaplen,
                                                 upper_chromosome=upper_genome_name,
                                                 upper_position=upper_genome_idx[0],
                                                 upper_ambiguity=upper_genome_idx[0] + overlaplen,
                                                 overlap_sequence=overlap,
                                                 untemplated_sequence=untemplated)

                        # breakpt = (upper_genome_idx[0] + overlaplen, upper_genome_idx[0],
                        #            lower_genome_idx[1], lower_genome_idx[1] - overlaplen)

                # Reverse complement the breakpoint, report as if scenario 1-4
                # breakpt = tuple(breakpt[::-1])

            return breakpt

        else:  # Guard against computing this with empty data
            raise ValueError("ContigInfo is not fully populated")


###################
# UTILITY FUNCTIONS
###################

def _number_range(a, b):
    def _len_common(stra, strb):
        c = 0
        for (a, b) in zip(stra, strb):
            if a == b:
                c += 1
            else:
                break
        return c
        #return len(''.join([x[0] for x in zip(stra, strb) \
        #             if reduce(lambda a, b:(a == b) and a or None,x)]))
    common = min(_len_common(str(a), str(b)), len(str(b)) - 2)
    return str(b)[common:]


def _str_info(info):
    def _get_read_string(read):
        read_strand = ('-' if read.is_reverse else '+')
        read_loc = read.reference_name
        read_lower = read.pos+1
        read_upper = read.aend+1
        return '{}{}:{}-{}'.format(read_strand, read_loc, read_lower, read_upper)

    return ' |...| '.join([_get_read_string(read) for read in info.reads])

    
def _read_fasta(filename):
    return [DNA(str(seq).upper(), seq.metadata) for seq in skbio.read(filename, format="fasta")]


def filecheck(filename):
    if not os.path.exists(filename):
        raise IOError('File {} not found'.format(filename))


def _parse_tigra_label(label):
    """
    Tigra labels are like this: Chr3_supercontig_000000381^2328502^Chr3_supercontig_000000381^2333647^INV^5145^++.Contig52
    This function extracts the first 4 '^' separated fields :- 
      Reference name (lower)
      Reference position (lower)
      Reference name (upper)
      Reference position (upper)
    """
    items = label.split('^')
    return ('^'.join(items[:4]), 
            LocationInfo(lower_name=items[0],
                         lower_pos=int(items[1]),
                         upper_name=items[2],
                         upper_pos=int(items[3])))


def _cigar_to_index(read):
    """
    Use the Cigar tuple to get indices of genome-aligned region w.r.t. contig
    Indices are SLICES, i.e. zero-based index of the first matching position
    and the first past-the-end position.
    Alignment length is rpos - lpos
    """
    MATCH = 0
    INS = 1
    DEL = 2
    REF_SKIP = 3
    SOFT_CLIP = 4
    HARD_CLIP = 5
    PAD = 6
    EQUAL = 7
    DIFF = 8

    length, begin, end = 0, 0, 0
    hasbegun = False
    
    for (state, size) in read.cigar:
        
        if state in (SOFT_CLIP, HARD_CLIP): # clipped end
            length += size
            if not hasbegun:
                begin += size
                end += size

        elif state == MATCH: # match
            hasbegun = True
            length += size
            end += size

        elif state in [INS, REF_SKIP]:
            length += size
            end += size

        elif state == DEL:
            pass

        elif state == REF_SKIP:
            length += size
            end += size

        elif state in [PAD, EQUAL, DIFF]:
            raise NotImplementedError('Didn\'t expect this state {}'.format(state))

        else:
            raise ValueError('Unrecognised state {} in {}'.format(state, read.cigar))

    if read.is_reverse:
        begin, end = length - end, length - begin

    return begin, end


def _ssw_to_index(read, contig):
    """
    Use Smith-Waterman alignment to get indices of genome-aligned region w.r.t. contig
    Indices are SLICES, i.e. zero-based index of the first matching position
    and the first past-the-end position.
    Alignment length is rpos - lpos
    """
    contigdb = SSW(contig)
    if read.is_reverse:
        aln = contigdb(_reverse_complement(read.query_alignment_sequence))
    else:
        aln = contigdb(read.query_alignment_sequence)
    return aln.query_begin, aln.query_end + 1


def _contig_index(read, contig):
    """
    Use exact matching to get indices of genome-aligned region w.r.t. contig
    Indices are SLICES, i.e. zero-based index of the first matching position
    and the first past-the-end position.
    Alignment length is rpos - lpos
    """
    seq = read.query_alignment_sequence.replace('-','')
    if read.is_reverse:
        rpos = len(contig) - contig.reverse_complement().index(read.query_alignment_sequence)
        lpos = rpos-len(seq)
    else:
        lpos = contig.index(read.query_alignment_sequence)
        rpos = lpos + len(seq)
    return lpos, rpos


def _reverse_complement(seq):
    import re
    replacements = dict(A='T', C='G', G='C', T='A',
                        a='t', c='g', g='c', t='a')
    pattern = '|'.join(sorted(re.escape(k) for k in replacements))
    return re.sub(pattern,
                  lambda m: replacements.get(m.group(0)),
                                             seq,
                                             flags=re.IGNORECASE)[::-1]


def _extract_reads(bamfile):
    sup = _get_supplementary_reads(bamfile)
    return _filter_reads(sup)


def _get_supplementary_reads(bamfile):
    filecheck(bamfile)
    bam = pysam.Samfile(bamfile, 'rb')
    names = set()
    for read in bam.fetch(until_eof=True):
        if read.is_supplementary:
            names.add(read.qname)

    reads = defaultdict(list)
    bam = pysam.Samfile(bamfile, 'rb')
    for read in bam.fetch(until_eof=True):
        if read.qname in names:
            reads[read.qname].append(read)

    return reads


def _filter_reads(matched_reads):
    filtered_reads = defaultdict(list)
    for name, readlist in matched_reads.items():

        contig_info = ContigInfo()
        id_, location = _parse_tigra_label(name)
        contig_info.location = location

        obs_matches = [False, False]
        for read in readlist:
            matches = location.is_close(read)
            if matches[0]:
                obs_matches[0] = True
            if matches[1]:
                obs_matches[1] = True

            if any(matches):
                contig_info.add_read(read)

        if (len(contig_info.reads) > 1) and (all(obs_matches)):
            filtered_reads[id_].append(contig_info)

    return filtered_reads


def _genome_index(read):
    """
    This is a 1-based INCLUSIVE index of the genome
    """
    return read.pos + 1, read.aend


def _break_type(left_read, right_read, left_index, right_index, contig, reverse=False):
    """
    Return any overlap or untemplated sequence at a breakpoint
    spanned by `left_read' and `right_read'
    """
    overlap = None
    overlaplen = 0
    untemplated = None

    if left_index[1] > right_index[0]:
        # Microhomology
        overlap = contig[right_index[0]:left_index[1]]
        if reverse:
            overlap = overlap.reverse_complement()
        overlaplen = len(overlap)

    elif left_index[1] < right_index[0]:
        # Untemplated sequence
        untemplated = contig[left_index[1]:right_index[0]]
        if reverse:
            untemplated = untemplated.reverse_complement()

    else:
        # Clean break
        pass

    return overlap, overlaplen, untemplated


def _sort_indices(container):
    """
    Return the indexes of the items in `container' AS IF it were
    sorted - not like numpy argsort
    """
    d = dict(zip(sorted(container), range(len(container))))
    return [d[k] for k in container]

def _range_overlap(x, y):
    """
    Return True if ranges x and y overlap. x and y must be well formed, i.e.
    x[0] <= x[1] && y[0] <= y[1]
    """
    if len(x) != 2 or len(y) != 2:
        raise ValueError("Ranges should have exactly two values, start and stop")

    if x[0] > x[1] or y[0] > y[1]:
        raise ValueError("Ranges must be in ascending order, i.e. x[0] <= x[1]")

    return x[0] <= y[1] and y[0] <= x[1]

def _checkindex(info):
    for i, read in enumerate(info.reads):
        bycigar = _cigar_to_index(read)
        byssw = _ssw_to_index(read, str(info.contig))
        bycontig = _contig_index(read, info.contig)
        print("Read {} - CIGAR {} SSW {} CONTIG {}".format(i, bycigar, byssw, bycontig))

def _batchify(collection):
    batches = []
    batch = []
    breaks = []

    if not isinstance(collection, Iterator):
        it = iter(collection)
    else:
        it = collection
    first_read = next(it)

    start, end = first_read.pos, first_read.aend
    clip = _clipside(first_read)
    batch.append(first_read)
    try:
        while True:
            next_read = next(it)
            s, e = next_read.pos, next_read.aend
            c = _clipside(next_read)
            if (s == start or e == end) and c == clip:
                batch.append(next_read)
            else:
                batches.append(batch)
                first_read = next_read
                batch = [first_read]
                start, end = s, e
                clip = c
    except StopIteration:
        pass
    return batches

def _batch_breakpoint(batch):
    starts, ends = zip(*[(read.pos+1, read.aend) for read in batch])
    clips = [_clipside(read) for read in batch]
    if len(set(starts)) == 1 and all(clip==(True, False) for clip in clips):
        return (True, batch[0].reference_name, starts[0])

    if len(set(ends)) == 1 and all(clip==(False, True) for clip in clips):
        return (True, batch[0].reference_name, ends[0])

    return (False, None, None)

def _clipside(read):
    """
    Which side is a read clipped on? Left, right, both or neither?
    """
    left = read.cigar[0][0] in (4, 5)
    right = read.cigar[-1][0] in (4, 5)
    return (left, right)

def _get_breaks(batches):
    breaks = []
    for batch in batches:
        success, chrom, breakpoint = _batch_breakpoint(batch)
        if success:
            breaks.append((chrom, breakpoint))
    return sorted(list(set(breaks)))

def _dist(bk1, bk2):
    if bk1[0] != bk2[0]:
        return np.inf
    dist = bk1[1] - bk2[1]
    if dist < 0:
        dist = -dist
    return dist

def _crossref_fn(bk, potential_breaks):
    import bisect
    lmatch = (bk.lower_chromosome, bk.lower_position) in potential_breaks
    umatch = (bk.upper_chromosome, bk.upper_position) in potential_breaks
    if lmatch and umatch:
        return True

    else:
        l = (bk.lower_chromosome, bk.lower_position)
        u = (bk.upper_chromosome, bk.upper_position)
        lpos = bisect.bisect_right(potential_breaks, l)
        upos = bisect.bisect_right(potential_breaks, u)
        if lpos >= len(potential_breaks):
            ldist = _dist(l, potential_breaks[lpos-1])
        else:
            ldist = min(_dist(l, potential_breaks[lpos-1]),
                        _dist(l, potential_breaks[lpos]))
        if upos >= len(potential_breaks):
            udist = _dist(l, potential_breaks[upos-1])
        else:
            udist = min(_dist(l, potential_breaks[upos-1]),
                        _dist(l, potential_breaks[upos]))
        if ldist < 5 and udist < 5:
            return True

    return False

def _location_match_fn(info, pair, **kwargs):
    a = info.location.is_close(info.reads[pair[0]], **kwargs)
    b = info.location.is_close(info.reads[pair[1]], **kwargs)
    return (a[0] or b[0]) and (a[1] or b[1])


if __name__ == '__main__':
    import sys
    import argparse

    epi='''
    There may be more than 1 inferred breakpoint for any given breakpoint prediction.
    This may be because there are multiple inputs, and the breakpoint looks slightly different
    each, or because there are multiple breakpoints close by each other, or because something
    has gone wrong. Use your judgement.
    '''
    parser = argparse.ArgumentParser(epilog=epi)
    parser.add_argument('bamfile', help='BAM file of contigs BWA mapped to reference')
    parser.add_argument('fastafile', help='FASTA file of complete contig sequences')
    parser.add_argument(
        '-l', '--limit',
        type=int,
        help='(Optional) limit of how far from the initial estimate a breakpoint can be inferred to be. Default = 5000',
        default=5000)

    args = parser.parse_args()

    sys.stderr.write("WARNING: This software is still under development and might not work %)\n")
    bam, fasta = args.bamfile, args.fastafile
    svs = SVPlacer(bam, fasta)
    breaks = svs.breakpoints(limit=args.limit)

    for (key, breaklist) in sorted(breaks.items()):
        print(key)
        for bk in set(str(x) for x in breaklist):
            print('\t{}'.format(bk))
    
