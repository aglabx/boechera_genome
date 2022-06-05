#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 01.07.2022
#@author: Aleksey Komissarov
#@contact: ad3002@gmail.com

from collections import defaultdict, Counter

from PyExp import sc_iter_filepath_folder
from trseeker.seqio.fasta_file import sc_iter_fasta_brute
from trseeker.models.trf_model import TRModel
from trseeker.tools.assembly_tools import get_n50


def aggregate_repeats_from_assemblies(trf_folder):
    assembly2repeats = defaultdict(list)
    repeat2locus = defaultdict(list)

    for file_path in sc_iter_filepath_folder(trf_folder, mask="genomic.trf"):
        name = file_path.split("/")[-1].split(".")[0]
        print(file_path)
        with open(file_path) as fh:
            for line in fh:
                d = line.strip().split("\t")
                tr = TRModel()
                tr.set_with_list(d)
                tr.project = name
                rep = tr.trf_consensus
                assembly2repeats[name].append(tr)
                repeat2locus[rep].append(tr)

    return assembly2repeats, repeat2locus


def get_contig_sizes(fasta_folder):
    assembly2contig_sizes = defaultdict(list)
    for file_path in sc_iter_filepath_folder(fasta_folder):
        print(file_path)
        name = file_path.split("/")[-1].split(".")[0]
        for header, seq in sc_iter_fasta_brute(file_path):
            assembly2contig_sizes[name].append(len(seq))
    return assembly2contig_sizes

def get_short_assembly_names(fasta_folder):

    assembly2name = {}
    for file_path in sc_iter_filepath_folder(fasta_folder):
        print(file_path)
        name = file_path.split("/")[-1].split(".")[0]
        for header, seq in sc_iter_fasta_brute(file_path):
            family, species = header.split()[1:3]
            assembly2name[name] = (family, species, f"{family} {species}")
            break
    return assembly2name

def print_assembly_trf_stats(assembly2contig_sizes, assembly2repeats):

    print("\t".join(map(str, ("Assembly", "Length (Mb)", "Contigs", "N50 (kb)", "L50", "Min contig", "Max contig (kb)", "All TRs", "Micro", "Large TRs (1kb)", "Large TRs (3kb)"))))
    for name in assembly2contig_sizes:
        n50, l50, min_contig, max_contig = get_n50(assembly2contig_sizes[name])
        size = round(sum(assembly2contig_sizes[name])/1000000, 1)
        n50 = round(n50/1000, 1)
        trs_all = len(assembly2repeats[name])
        trs_micro = len([x for x in assembly2repeats[name] if x.trf_period < 5 and x.trf_array_length < 500])
        trs_1kb = len([x for x in assembly2repeats[name] if x.trf_array_length > 1000])
        trs_3kb = len([x for x in assembly2repeats[name] if x.trf_array_length > 3000])
        print("\t".join(map(str, (name, size, len(assembly2contig_sizes[name]), n50, l50, min_contig, max_contig, trs_all, trs_micro, trs_1kb, trs_3kb))))


if __name__ == '__main__':

    fasta_folder = "/media/eternus1/nfs/projects/boechera/assembly"
    trf_folder = "/media/eternus1/nfs/projects/boechera/trf"

    assembly2repeats, repeat2locus = aggregate_repeats_from_assemblies(trf_folder)
    assembly2contig_sizes = get_contig_sizes(fasta_folder)
    
    assembly2name = get_short_assembly_names(fasta_folder)

    print_assembly_trf_stats(assembly2contig_sizes, assembly2repeats)

    print(len(repeat2locus), len(assembly2contig_sizes))

    '''Мы обнаружили 6 616 662 различных тандемных повторов в 75 собранных геномах. 
    Из них 5 577 378 повторов являются найденными только один раз и с длиной поля тандемного повтора короче 1000 bp. После фильтрации такие повторов в анализе осталось 1 039 285 повторов. Нашей основной целью был поиск больших тандемных повторов характерных для таких областей генома как центромерный и перицентрометный районы хромосом. Для этого мы отобрали только те повторы для хотя одного локуса которого собрано хотя бы одное поле длинней 3000 bp и количеством копий мономера тандемного повтора более 10. Таких повторов оказалось 54069 штук без учета количества копий и 29919 штук с учетом количества копий. Эти повторы мы отсортировали по длине поля тандемного повтора и проанализировали самые большие.'''

    ok = 0
    removed_single_short = 1
    for repeat in repeat2locus:
        N = len(repeat2locus[repeat])
        arrays_bp = sum([x.trf_array_length for x in repeat2locus[repeat]])
        if N == 1 and arrays_bp < 1000:
            removed_single_short += 1
            continue
        ok += 1

    print(f"OK {ok}, removed_single_short {removed_single_short}")


    large_trs = {}
    has_large_n = 0
    for repeat in repeat2locus:
        N = len(repeat2locus[repeat])
        has_large = sum([1 for x in repeat2locus[repeat] if x.trf_array_length > 3000 and x.trf_n_copy > 10])
        if not has_large:
            continue
        has_large_n += 1
        large_trs[repeat] = repeat2locus[repeat][::]

    print(f"has_large_n {has_large_n}")

    species_specific = 0
    common_trs = defaultdict(list)
    for repeat in large_trs:
        assemblies = len(set([x.project for x in repeat2locus[repeat]]))
        if assemblies == 1:
            species_specific += 1
        else:
            common_trs[repeat] = repeat2locus[repeat]
    #     print(repeat, len(repeat), max([x.trf_array_length for x in repeat2locus[repeat]]), assemblies)
        
    print(f"species_specific {species_specific} common {has_large_n-species_specific}")

    table = []
    taxons = [(assembly,family, species, both) for assembly, (family, species, both) in assembly2name.items()]
    taxons.sort(key=lambda x: x[-1])
    for repeat in common_trs:
        number_of_projects = len(set([x.project for x in repeat2locus[repeat]]))
        max_array_length = max([x.trf_array_length for x in repeat2locus[repeat]])
        projects = Counter([x.project for x in repeat2locus[repeat]])
        if number_of_projects > 10:
            print(repeat.upper(), number_of_projects, max_array_length, projects)
            row = [repeat.upper(), number_of_projects, max_array_length] + ["" for x in range(len(assembly2name))]
            for i,(assembly,family, species, both) in enumerate(taxons):
                if assembly in projects:
                    row[i+3] = projects[assembly]
            table.append(row)

    header = [""]
    sizes1 = [""]
    sizes2 = [""]
    table.sort(key=lambda x: -x[1])
    for j, (rep, number_of_projects, max_array_length, *items) in enumerate(table):
        header.append(rep)
        sizes1.append(number_of_projects)
        sizes2.append(max_array_length)
    trans_table = [header, sizes1, sizes2]

    for i,(v,a,b,c) in enumerate(taxons):
        row = [c]
        for _row in table:
            row.append(_row[3+i])
        trans_table.append(row)

    for row in trans_table:
        print("\t".join(map(str, row)))

    '''Из 29919 больших тандемных повтров 28509 повторов были видоспецифичные. И только 1410 повторов встречались более чем в одной сборке.'''

    period2n = defaultdict(int)
    for repeat in large_trs:
        period2n[len(repeat)] += 1

    keys = list(period2n.keys())
    keys.sort(key=lambda x: period2n[x], reverse=True)

    for i,key in enumerate(keys):
        print(f"{i+1}\t{key}\t{period2n[key]}")