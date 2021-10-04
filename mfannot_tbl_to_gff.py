#!/usr/bin/env python3
'''
- Multi-exon rnl is missing usually, parsing rnl gene might be erroneous.
- Feature table (tbl) format: https://www.ncbi.nlm.nih.gov/genbank/feature_table
'''

import sys

fin = sys.argv[1]

class MFGene:
    def __init__(self):
        self.chr = ''
        self.start = 0
        self.end = 0
        self.gene = ''
        self.strand = '+'
        self.type = 'mRNA'
        self.note = ''
        self.product = ''
        self.protein_id = ''
        self.anticodon = ''
        self.cdsCoord = []
        self.cdsNote = ''
        self.intronCoord = []
        self.intronNote = []

    def __str__(self):
        provider = 'mfannot'
        out = []
        if self.type == 'mRNA':
            # gene
            info = f'ID={self.protein_id};gene={self.gene}'
            if self.note:
                info += f';note={self.note}'
            out.append(f'{self.chr}\t{provider}\tgene\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{info}')
            # mRNA
            info = f'ID={self.protein_id}:mRNA;gene={self.gene};product={self.product}'
            # for orf, should be only one CDS
            if self.gene.startswith('orf') and self.cdsNote:
                info += f';note={self.cdsNote}'
            info += f';Parent={self.protein_id}'
            out.append(f'{self.chr}\t{provider}\tmRNA\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{info}')
            # CDS
            outcds = []
            for n, c in enumerate(self.cdsCoord, 1):
                info = f'ID={self.protein_id}:mRNA:CDS:{n};Parent={self.protein_id}:mRNA'
                outcds.append(f'{self.chr}\t{provider}\tCDS\t{c[0]}\t{c[1]}\t.\t{self.strand}\t.\t{info}')
            # exon
            outexon = []
            for n, c in enumerate(self.cdsCoord, 1):
                info = f'ID={self.protein_id}:mRNA:exon:{n};Parent={self.protein_id}:mRNA'
                outexon.append(f'{self.chr}\t{provider}\texon\t{c[0]}\t{c[1]}\t.\t{self.strand}\t.\t{info}')
            # intron
            outintron = []
            if len(self.intronCoord) > 0:
                for n, (c, note) in enumerate(zip(self.intronCoord, self.intronNote), 1):
                    info = f'ID={self.protein_id}:mRNA:intron:{n};note={note};Parent={self.protein_id}:mRNA'
                    outintron.append(f'{self.chr}\t{provider}\tintron\t{c[0]}\t{c[1]}\t.\t{self.strand}\t.\t{info}')
            # flip if complementary strand
            if self.strand == '-':
                out = out + outcds[::-1] + outexon[::-1] + outintron[::-1]
            else:
                out = out + outcds + outexon + outintron
        elif self.type == 'tRNA':
            # gene
            info = f'ID={self.protein_id};gene={self.gene}'
            if self.note:
                info += f';note={self.note}'
            out.append(f'{self.chr}\t{provider}\tgene\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{info}')
            # tRNA
            ac = self.gene[5:8]
            pos = self.anticodon[4:-5]
            aa = self.anticodon[-1].upper()
            info = f'ID={self.protein_id}:tRNA;gene={self.gene};product={self.product};anticodon={ac};pos={pos};aa={aa};Parent={self.protein_id}'
            out.append(f'{self.chr}\t{provider}\ttRNA\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{info}')
            # exon
            info = f'ID={self.protein_id}:tRNA:exon:1;Parent={self.protein_id}:tRNA'
            out.append(f'{self.chr}\t{provider}\texon\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{info}')
        elif self.type == 'rRNA':
            # gene
            info = f'ID={self.protein_id};gene={self.gene}'
            if self.note:
                info += f';note={self.note}'
            out.append(f'{self.chr}\t{provider}\tgene\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{info}')
            # rRNA
            info = f'ID={self.protein_id}:rRNA;gene={self.gene};product={self.product};Parent={self.protein_id}'
            out.append(f'{self.chr}\t{provider}\trRNA\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{info}')
            # exon
            for n, c in enumerate(self.cdsCoord, 1):
                info = f'ID={self.protein_id}:rRNA:exon:{n};Parent={self.protein_id}:rRNA'
                out.append(f'{self.chr}\t{provider}\texon\t{c[0]}\t{c[1]}\t.\t{self.strand}\t.\t{info}')
            # intron
            if len(self.intronCoord) > 0:
                for n, (c, note) in enumerate(zip(self.intronCoord, self.intronNote), 1):
                    info = f'ID={self.protein_id}:rRNA:intron:{n};note={note};Parent={self.protein_id}:rRNA'
                    out.append(f'{self.chr}\t{provider}\tintron\t{c[0]}\t{c[1]}\t.\t{self.strand}\t.\t{info}')
        return '\n'.join(out)

with open(fin) as f:
    gene = None
    genes = []
    for line in f:
        if line.startswith('>'):
            CHROM = line.split()[1]
            continue
        line = line.rstrip().split('\t')
        if len(line) == 3:
            if line[2] == 'gene':
                feat = 'gene'
                if gene:
                    genes.append(gene)
                gene = MFGene()
                gene.chr = CHROM
                p = [int(line[0]), int(line[1])]
                if p[0] > p[1]:
                    gene.strand = '-'
                gene.start, gene.end = sorted(p)
                # next line is gene name
                line = next(f).rstrip().split('\t')
                gene.gene = line[4]
            elif line[2] == 'CDS':
                feat = 'CDS'
                gene.cdsCoord.append(sorted([int(line[0]), int(line[1])]))
            elif line[2] == 'exon':
                feat = 'exon'
                # same to cds, just skip
                next(f)
            elif line[2] == 'intron':
                feat = 'intron'
                gene.intronCoord.append(sorted([int(line[0]), int(line[1])]))
            elif line[2] == 'tRNA':
                feat = 'tRNA'
                gene.type = 'tRNA'
            elif line[2] == 'rRNA':
                feat = 'rRNA'
                gene.type = 'rRNA'
                gene.cdsCoord.append(sorted([int(line[0]), int(line[1])]))
        elif len(line) == 5:
            if feat == 'gene':
                if line[3] == 'note':
                    # eg. copy 1
                    gene.note = line[4]
            elif feat == 'CDS':
                if line[3] == 'product':
                    gene.product = line[4]
                elif line[3] == 'protein_id':
                    gene.protein_id = line[4][7:]
                elif line[3] == 'note':
                    # only in orf
                    gene.cdsNote = line[4]
            elif feat == 'intron':
                if line[3] == 'note':
                    gene.intronNote.append(line[4])
            elif feat == 'tRNA':
                if line[3] == 'product':
                    gene.product = line[4]
                elif line[3] == 'protein_id':
                    gene.protein_id = line[4][7:]
                elif line[3] == 'anticodon':
                    gene.anticodon = line[4].strip('()')
            elif feat == 'rRNA':
                if line[3] == 'product':
                    gene.product = line[4]
                elif line[3] == 'protein_id':
                    gene.protein_id = line[4][7:]
        elif len(line) == 2:
            if feat == 'CDS':
                gene.cdsCoord.append(sorted([int(line[0]), int(line[1])]))
            else:
                # maybe multi-exon rnl
                raise NotImplementedError(line)

genes.sort(key=lambda x: x.start)
for g in genes:
    # print(g.__dict__)
    print(g)
