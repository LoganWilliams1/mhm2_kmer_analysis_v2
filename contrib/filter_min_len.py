#!/usr/bin/env python3

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

import gzip
import argparse

"""
This program takes an excessively large assembly (with > 2^31 scaffolds) and its sam[.gz] files
calculates the minimum_scaffold_length that <= max_num_scaffolds
and produces a filtered assembly and filtered sam.gz file using the new subset assembly
"""

include_unmapped = False

def maybe_gz_open(filename:str, mode:str):
    if filename.endswith(".gz"):
        return gzip.open(filename, mode + 't')
    else:
        return open(filename, mode)

def fold(fasta):
    return fasta

class Header():
    """
    SAM header
    """
    def __init__(self):
        self.hd = list() # list of @HD
        self.sq = dict() # dict of @SQ (name: dict().  Ordering is important and maintained)
        self.rg = list() # list of @RG
        self.pg = list() # list of @PG
        self.co = list() # list of @CO
        self.histogram = dict()

    def copy(self):
        """Return a new header"""
        header = Header()
        header.hd = self.hd.copy()
        header.sq = self.sq.copy()
        header.rg = self.rg.copy()
        header.pg = self.pg.copy()
        header.co = self.co.copy()
        return header

    def add_line(self, line):
        line = line.rstrip()
        if line.startswith("@HD"):
            self.hd.append(line)
        elif line.startswith("@SQ"):
            fields=line.split()
            sq = dict()
            for field in fields[1:]:
                logger.debug(f"field {field} in line {line}")
                k,v = field.split(":")
                sq[k] = v
            assert "SN" in sq
            assert "LN" in sq

            self.sq[sq["SN"]] = sq
            length = int(sq["LN"])
            if length not in self.histogram:
                self.histogram[length] = 1
            else:
                self.histogram[length] += 1
        elif line.startswith("@RG"):
            self.rg.append(line)
        elif line.startswith("@PG"):
            self.pg.append(line)
        elif line.startswith("@CO"):
            self.co.append(line)
        else:
            logger.error(f"Unknown header line: {line}")
            raise 

    def dump(self, fh):
        for line in self.hd:
            print(line, file=fh)
        for k,d in self.sq.items():
            line = f"@SQ"
            for k,v in d.items():
                line += f"\t{k}:{v}"
            print(line, file=fh)
        for line in self.rg:
            print(line, file=fh)
        for line in self.pg:
            print(line, file=fh)
        for line in self.co:
            print(line, file=fh)

    def load_header(self, sam_fh) -> str:
        """
        Loads a header from a sam file handle
        returns the first non-header line
        """
        logger.info(f"Reading header from sam file handle {sam_fh}")
        lines = 0
        for line in sam_fh:
            if not line[0] == "@":
                logger.info(f"Read {lines} header lines")
                return line
            self.add_line(line)
            lines += 1
        return None

    def filter_min_len(self, min_len:int):
        """
        Returns a new Header (with sequences sorted largest to smallest)
        returning only those scaffolds >= min_len
        """
        logger.info(f"Filtering {len(self.sq)} sequences with new header and min_len={min_len}")
        new_header = self.copy()
        new_sq = list()
        for name,sq in new_header.sq.items():
            if int(sq["LN"]) >= min_len:
                new_sq.append(sq)
        logger.info(f"Filtered {len(new_sq)} sequences with min_len={min_len} from {len(new_header.sq)}")
        new_sq.sort(key=lambda x: int(x["LN"]), reverse=True)
        logger.info(f"Sorted first is {new_sq[0]} last is {new_sq[-1]}")
        new_header.sq = dict()
        for sq in new_sq:
            new_header.sq[sq["SN"]] = sq
        return new_header

    def calculate_min_len(self, max_num_scaffolds:int = None) -> int:
        """Calculates the longest min_len that results in <= max_num_scaffolds"""
        if max_num_scaffolds is None:
            max_num_scaffolds = 2147483647
        logger.info(f"Finding the min_len that leaves at most {max_num_scaffolds}, starting with {len(self.sq)}")
        min_len = 0
        bins = sorted(self.histogram.keys())
        stop = 2 ** 31 -1
        if len(bins):
            stop = bins[-1]
        logger.info(f"bins={self.histogram}")
        num = len(self.sq)
        while num > max_num_scaffolds:
            if min_len > stop:
                break
            if min_len in self.histogram:
                num -= self.histogram[min_len]
            min_len += 1
        logger.info(f"Found that min_len={min_len} yields {num} scaffolds out of {len(self.sq)}")
        return min_len

    def filter_sam(self, sam_file_name:str, new_file_name:str):
        """
        Writes a new sam file replacing the header and filtering out lines not covered by the new header
        """
        
        global include_unmapped
        is_good = True
        out_fh = maybe_gz_open(new_file_name, 'w')
        pass_records = 0
        records = 0
        logger.info(f"Filtering SAM {sam_file_name} writing to {new_file_name}")
        with maybe_gz_open(sam_file_name, 'r') as in_fh:
            old_header = Header()
            line = old_header.load_header(in_fh)
            for name,sq in self.sq.items():
                if name not in old_header.sq:
                    logger.warning(f"No entry of {name} from {sam_file_name} is in this header")
                    is_good = False
            if not is_good:
                logger.warning(f"Please ensure {sam_file_name} is from the same assembly")
                raise
            self.dump(out_fh)
            while line is not None and len(line) > 0:
                fields = line.split()
                #logger.debug(f"processing {fields}")
                if len(fields) < 3:
                    logger.warning(f"what is {line} {fields}")
                    continue
                if fields[2] is None or fields[2] == "*":
                    if include_unmapped:
                        print(line, file=out_fh)
                        pass_records += 1
                elif fields[2] in self.sq:
                    print(line, file=out_fh)
                    pass_records += 1
                records += 1
                line = in_fh.readline().rstrip()
        logger.info(f"Filtered {pass_records} of {records}")
        out_fh.close()

    def filter_assembly(self, assem_file_name:str, new_assem_file_name:str):
        """
        Prints the filtered assembly in the same order as this header
        Skipping any records not in this header
        """
        logger.info(f"Filtering assembly {assem_file_name} saving to {new_assem_file_name}")
        printable = dict()
        out_fh = maybe_gz_open(new_assem_file_name, 'w')
        printed = 0
        with maybe_gz_open(assem_file_name, 'r') as in_fh:
            # print in the correct order
            last_n = ""
            name = ""
            fasta = ""
            print_it = False
            line_num = 0
            for tgt,sq in self.sq.items():
                logger.debug(f"Looking for {tgt} at {line_num}")
                if tgt in printable:
                    logger.debug(f"Found {tgt} before {line_num}")
                    print(f"{printable[tgt]}", file=out_fh)
                    del printable[tgt]
                    printed += 1
                    continue

                # read and remember file until tgt is found
                for l in in_fh:
                    l = l.rstrip()
                    line_num += 1
                    if l[0] == ">":
                        fields = l.split()
                        n = fields[0][1:] # name minus '>' prefix and anything after the first whitespace
                        if tgt == last_n:
                            assert print_it
                            print(f"{name}\n{fold(fasta)}", file=out_fh)
                            logger.debug(f"Hit {tgt} at {line_num}")
                            printed += 1
                            name = l
                            last_n = n
                            fasta = ""
                            print_it = n in self.sq
                            break
                        elif print_it:
                            printable[last_n] = f"{name}\n{fold(fasta)}"
                            logger.debug(f"Recorded {last_n} at {line_num}")
                        name = l
                        last_n = n
                        fasta = ""

                        if n in self.sq:
                            print_it = True
                            name = l
                            fasta = ""
                            logger.debug(f"Keeping {self.sq[n]}")
                        else:
                            print_it = False
                            name = None
                            fasta = ""
                            logger.debug(f"Skipping {n}")
                    else:
                        fasta += l
            if print_it:
                printable[last_n] = f"{name}\n{fasta}"

        if printed != len(self.sq):
            logger.warning(f"Should not get here!! name={name} last_n={last_n} printed={printed} seqs={len(self.sq)}")
                
        out_fh.close()
            
if __name__ == "__main__":
    argparser = argparse.ArgumentParser(add_help=True, prog = "filter_min_len.py", description = """
    Filters SAM[.gz] and assembly.fa[.gz] files where all references have a minimum length
    Optionally determines the optimal min_len to conform to the BAM standard (<=2^31-1 scaffolds)
    """, epilog="")
    argparser.add_argument("--assembly", default=None, type=str, help="Optional. The assembly file to filter")
    argparser.add_argument("--min-len", default=None, type=int, help="If not set (default), then determine the optimal minimum scaffold length")
    argparser.add_argument("--max-scaffolds", default=None, type=int, help="If not set (default), then use the maximum allowd by the BAM sepcification")
    argparser.add_argument("--gzip-output", action='store_true', help="If set then the filtered files will be compressed")
    argparser.add_argument("sam", nargs='+', type=str, help="one or more sam files to process serially")
    options, unknown_options = argparser.parse_known_args()
    logger.info(f"options {options} unknown {unknown_options}")

    if len(options.sam) == 0:
        logger.warning(f"You must specify at least one SAM file")
        raise
    
    initial_header = Header()
    sam_fh = maybe_gz_open(options.sam[0], 'r')
    initial_header.load_header(sam_fh)
    sam_fh.close()

    min_len = options.min_len
    if min_len is None:
        min_len = initial_header.calculate_min_len(options.max_scaffolds)
    else:
        min_len = int(min_len)

    suffix = ".gz" if options.gzip_output else ""
    new_header = initial_header.filter_min_len(min_len)
    if options.assembly:
        new_header.filter_assembly(options.assembly, f"{options.assembly}-min{min_len}.fa{suffix}")

    for sam in options.sam:
        new_header.filter_sam(sam, f"{sam}-min{min_len}.sam{suffix}")
    


