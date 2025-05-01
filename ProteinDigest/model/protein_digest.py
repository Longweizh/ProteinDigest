import re


class Peptide:
    def __init__(self, name, seq, miss, l, r):
        self.name = name
        self.seq = seq
        self.miss = miss
        self.l = l
        self.r = r

    def length(self):
        return len(self.seq)

    def mw(self):
        return mol_weight(self.seq)


def peptide_slice(slice, miss):
    results = {}
    for i in range(1, miss + 2):
        result_list = ["".join(slice[j:j + i]) for j in range(len(slice) - i + 1)]
        if result_list:
            results[i - 1] = set(result_list)
    return results


def protein_digest(seq, enzyme, l_min, l_max, mw_min, mw_max, miss):
    seq = seq.upper()
    seq = seq.replace(' ','')
    pattern_data = {"Trypsin": re.compile(r'(?<=[KR])(?=[^P])'),
                    "Trypsin (C-term to K/R, even before P)": re.compile(r'(?<=[KR])(?=.)'),
                    "Lys C": re.compile(r'(?<=[K])(?=.)'),
                    "Lys N": re.compile(r'(?<=.)(?=K)'),
                    "CNBr": re.compile(r'(?<=[M])(?=.)'),
                    "Arg C": re.compile(r'(?<=[R])(?=[^P])'),
                    "Asp N": re.compile(r'(?<=.)(?=D)'),
                    "Glu C (bicarbonate)": re.compile(r'(?<=[E])(?=[^PE])'),
                    "Glu C (phosphate)": re.compile(r'(?<=[DE])(?=[^PE])'),
                    "Microwave-assisted formic acid hydrolysis (C-term to D)": re.compile(r'(?<=[D])(?=.)'),
                    "Pepsin (pH 1.3)": re.compile(r'(?<=[FL])(?=.)'),
                    "Pepsin (pH > 2)": re.compile(r'(?<=[FLWYAEQ])(?=.)'),
                    "Proteinase K": re.compile(r'(?<=[AFYWLIV])(?=.)'),
                    "Thermolysin": re.compile(r'(?<=[DE])(?=[AFILMV])')}
    pattern = pattern_data[enzyme]
    slice = pattern.split(seq)
    peptide_list = peptide_slice(slice, miss)
    candidates = []
    results = []
    for miss_number in peptide_list.keys():
        for peptide in peptide_list[miss_number]:
            pos_start = [i.start() for i in re.finditer(peptide,seq)]
            s = set()
            for i in pos_start:
                if i == 0:
                    s.add('Start')
                else:
                    s.add(seq[i-1])
            pos_end = [i.end() for i in re.finditer(peptide, seq)]
            e = set()
            for i in pos_end:
                try:
                    e.add(seq[i])
                except IndexError:
                    e.add('End')
            a = Peptide(name=peptide, seq=peptide, miss=miss_number, l=s, r=e)
            candidates.append(a)

    for candidate in candidates:
        if int(l_min) <= candidate.length() <= int(l_max) and int(mw_min) <= candidate.mw() <= int(mw_max):
            results.append(candidate)

    results = sorted(results, key=lambda x:x.mw(), reverse=True)
    return results


def mol_weight(seq):
    amino_acid_masses = {
        'A': 89.1,
        'R': 174.2,
        'N': 132.1,
        'D': 133.1,
        'C': 121.2,
        'E': 147.1,
        'Q': 146.2,
        'G': 75.1,
        'H': 155.2,
        'I': 131.2,
        'L': 131.2,
        'K': 146.2,
        'M': 149.2,
        'F': 165.2,
        'P': 115.1,
        'S': 105.1,
        'T': 119.1,
        'U': 204.2,
        'W': 204.2,
        'Y': 181.2,
        'V': 117.1,
    }
    mw_seq = round(sum(amino_acid_masses[aa] for aa in seq) - 18.0 * (len(seq) - 1),1)
    return mw_seq


a = 'MGKKQNRKTGNSKTQSASPPPKERSSSPATEQSWMENDFDELREEGFRRSNYSELREDIQTKGKEVENFEKNLEECITRITNTEKCLKELMELKTKARELREECRSLRSRCDQLEERVSAMEDEMNEMKREGKFREKRIKRNEQSLQEIWDYVKRPNLRLIGVPESDVENGTKLENTLQDIIQENFPNLARQANVQIQEIQRTPQRYSSRRATPRHIIVRFTKVEMKEKMLRAAREKGRVTLKGKPIRLTADLSAETLQARREWGPIFNILKEKNFQPRISYPAKLSFISEGEIKYFIDKQMLRDFVTTRPALKELLKEALNMERNNRYQPLQNHAKM'
slice = protein_digest(a, "Trypsin", 1, 100000, 0, 100000000, 5)
for i in slice:
    print(i.seq, i.l, i.r)