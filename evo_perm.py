import random
import copy
import time
import multiprocessing

# spectrum = []  # fragments of DNA sequence
# cached_distances = {}  # cached overlapping value for dna fragments


def compute_back_track_table(s):
    t = [0] * len(s)
    cnd, pos, t[0], t[1] = 0, 2, -1, 0
    while pos < len(s):
        if s[pos - 1] == s[cnd]:
            t[pos] = cnd + 1
            pos += 1
            cnd += 1
        elif cnd > 0:
            cnd = t[cnd]
        else:
            t[pos] = 0
            pos += 1
    return t


class GeneticAlgorithm:

    def __init__(self, spectrum, population_size, iterations_without_improvement, mutation_chance):
        self.spectrum = spectrum
        self.cached = {}
        for s in self.spectrum:
            self.cached[s] = {}
        self.n = len(self.spectrum)
        self.l = len(self.spectrum[0])
        self.seq_limit = self.n + self.l - 1
        self.size = population_size
        self.iterations_to_wait = iterations_without_improvement
        self.mutation_chance = mutation_chance

        self.population = []
        self.best_solution = []
        self.best_value = 0
        self.best_gene = None

    def compute(self):  # main function
        self.population = [Gene() for _ in range(self.size)]  # generate population
        self.best_gene = copy.deepcopy(self.population[0])  # assume first one is best

        no_progress = 0  # iterations without progress
        while self.iterations_to_wait - no_progress > 0:  # run until there is no improvement

            # sort population basing on length of their sequence
            self.population.sort(key=lambda x: x.get_solution(self.seq_limit)[0], reverse=True)

            no_progress += 1  # assume there is no progress
            if self.is_improved():  # check if progress was made
                best_value = self.best_gene.get_solution(self.seq_limit)
                print(no_progress, best_value)
                if best_value[0] == self.n:  # finish if the solution uses entire spectrum
                    break
                no_progress = 0  # reset progress count

            self.create_new_generation()  # perform crossovers and add random new genes
            self.perform_mutations()  # perform mutations basing on a given chance percentage
        return self.best_gene  # return the best gene ever!

    def is_improved(self):
        # check whether present best gene uses more spectrum elements than the best gene ever
        if self.population[0].get_solution(self.seq_limit)[0] > self.best_gene.get_solution(self.seq_limit)[0]:
            self.best_gene = copy.deepcopy(self.population[0])
            return True
        return False

    def create_new_generation(self):
        crossovers = 0  # count the number of crossovers performed
        while crossovers < self.size // 4:
            first_parent = second_parent = random.randint(0, self.size - crossovers - 1)
            while first_parent == second_parent:
                second_parent = random.randint(0, self.size - crossovers - 1)
            g1, g2 = self.population[first_parent], self.population[second_parent]
            p = multiprocessing.Process(target=Gene.crossover,
                                        args=(g1, g2, self.population, self.size // 2 + crossovers))
            p.start()
            p.run()
            # Gene.crossover(g1, g2, self.population, self.size // 2 + crossovers)
            crossovers += 1
        for i in range(self.size // 2 + crossovers, self.size):
            self.population[i] = Gene()
            # p.join()

    def perform_mutations(self):
        for gene in self.population:
            if random.random() < self.mutation_chance:
                gene.mutate()

    def overlap(self, string1, string2):
        c = self.cached[string1].get(string2, None)
        if c is not None:
            return c

        t = compute_back_track_table(string2)
        m, i = 0, 0
        while m + i < len(string1):
            if string2[i] == string1[m + i]:
                i += 1
            else:
                m += i - t[i]
                if i > 0:
                    i = t[i]
        self.cached[string1][string2] = i
        return i


class Gene:
    @staticmethod
    def crossover(g1, g2, arr, ind):
        neigh_list = {}
        length = len(g1.solution)
        for i, base in enumerate(g1.solution):
            neigh_list[base] = [g1.solution[i - 1], g1.solution[(i + 1) % length]]
        for i, base in enumerate(g2.solution):
            neigh_list[base].append(g2.solution[i - 1])
            neigh_list[base].append(g2.solution[(i + 1) % length])

        neigh_chosen = [g1.solution[0], g2.solution[0]][random.randint(0, 1)]
        child = [neigh_chosen]

        while len(child) < length:
            min_key = 5
            min_neigh_list = []
            for k in neigh_list:
                if neigh_chosen in neigh_list[k]:
                    neigh_list[k].remove(neigh_chosen)
                if k in neigh_list[neigh_chosen]:
                    if len(neigh_list[k]) < min_key:
                        min_key = len(neigh_list[k])
                        min_neigh_list = [k]
                    elif len(neigh_list[k]) == min_key:
                        min_neigh_list.append(k)
            del neigh_list[neigh_chosen]
            if len(min_neigh_list) > 0:
                neigh_chosen = max(min_neigh_list, key=lambda x: overlapping(neigh_chosen, x))
            else:
                neigh_chosen = max(neigh_list.keys(), key=lambda x: overlapping(neigh_chosen, x))
            child.append(neigh_chosen)
        arr[ind] = Gene(child)

    def __init__(self, permutation=None):
        self.solution = permutation
        if self.solution is None:
            self.solution = list(spectrum)
            random.shuffle(self.solution)
        self.valid = False
        self.s_valid = False
        self.best_value = 0
        self.best_solution = None

    def solution_value(self):
        if not self.valid:
            self.best_value = 0
            for codon in range(len(self.solution) - 1):
                self.best_value += overlapping(self.solution[codon], self.solution[codon + 1])
            self.valid = True
        return self.best_value

    def get_solution(self, limit):
        if self.s_valid:
            return self.best_solution
        sol = self.solution[0]
        fits = [overlapping(self.solution[s], self.solution[s + 1]) for s in range(len(self.solution) - 1)]
        max_slice = (0, 0)
        for slice_s in range(len(fits)):
            for slice_f in range(slice_s + max_slice[1] - max_slice[0], len(fits)):
                total = sum(fits[slice_s:slice_f])
                if (slice_f - slice_s + 1) * l - total >= limit:
                    break
                if slice_f - slice_s > max_slice[1] - max_slice[0]:
                    max_slice = (slice_s, slice_f)
        for s in range(max_slice[0], max_slice[1] + 1):
            common = overlapping(self.solution[s], self.solution[s + 1])
            sol += self.solution[s + 1][common:]
        self.best_solution = max_slice[1] - max_slice[0] + 2, max_slice, sol[:limit], len(sol[:limit])
        self.s_valid = True
        return max_slice[1] - max_slice[0] + 2, max_slice, sol[:limit], len(sol[:limit])

    def mutate(self):
        solution = self.get_solution(n + l - 1)
        seq_s, seq_f = solution[1]
        minimal = (-1, 0, -1)  # index / with value / at index
        for i in range(seq_s, min(seq_f, len(self.solution) - 2)):
            if overlapping(self.solution[i], self.solution[i + 1]) == l - 1:
                continue
            for j in range(0, seq_s):
                total_overlap = overlapping(self.solution[i], self.solution[j]) + overlapping(self.solution[j],
                                                                                              self.solution[i + 2])
                if minimal[1] < total_overlap:
                    minimal = (j, total_overlap, i + 1)
            for j in range(seq_f, len(self.solution)):
                total_overlap = overlapping(self.solution[i], self.solution[j]) + overlapping(self.solution[j],
                                                                                              self.solution[i + 2])
                if minimal[1] < total_overlap:
                    minimal = (j, total_overlap, i + 1)
            if minimal[1] >= (l - 1) * 2:
                break
        if minimal[0] >= 0:
            rnd1 = minimal[0]
            rnd2 = minimal[2]
            # print(rnd1, rnd2, minimal[1])
            # print(seq_s, seq_f)
            # print(self.solution[rnd1], 'at', self.solution[rnd2-1], self.solution[rnd2], self.solution[rnd2+1])
            temp = self.solution[rnd1]
            self.solution[rnd1] = self.solution[rnd2]
            self.solution[rnd2] = temp
            self.valid = False
            # rnd1 = random.randint(0, len(self.solution) - 1)
            # rnd2 = random.randint(0, len(self.solution) - 1)
            # temp = self.solution[rnd1]
            # self.solution[rnd1] = self.solution[rnd2]
            # self.solution[rnd2] = temp
            # self.valid = False
            # self.s_valid = False


if __name__ == '__main__':
    dict_test = {'end/': 12}  # ,
    # 'neg/': 22,
    # 'pos/': 12,
    # 'rep/': 5}
    start_path = 'inst/'

    for k, v in dict_test.items():
        for i in range(1, v + 1):
            file_name = str.format('{0}{1}{2}.bin', start_path, k, i)
            spectrum = []
            # cached_distances = {}
            file = open(file_name)
            line = file.readline().replace('\n', '')
            while len(line) > 1:
                spectrum.append(line)
                # cached_distances[line] = {}
                line = file.readline().replace('\n', '')
            __n = len(spectrum)
            __l = len(spectrum[0])
            ppl = 100
            print(str.format('file: {0}, n = {1}, l = {2}', file_name, __n, __l))
            print('Population: ' + str(ppl))
            start_stamp = time.time()
            result_gene = GeneticAlgorithm(spectrum, ppl, __n ** 2, 0.001).compute()
            print(str.format('{0} elements out of {1} in {2} seconds', result_gene.get_solution(__n + __l - 1)[0], __n,
                             time.time() - start_stamp))

            # RAW_DATA = 'AGTTCGCCCCAGTAATGTTGCCAATAAGGACCACCAAATCCGCATGTTACAGGACTTCTTATAAATTCTTTTTTCGTGGGGAGCAGCGGATCTTAATGGATGGCGCCAGCTGGTATGGAAGCTAATAGCGCCGGTGAGAGGGTAATCAGCCGTCTCCACCAACACAACGCTATCGGGTCATATTATAAGATTCCGCAATGGGACTACTTATAGGTTGCCTTAACGATATCCGCAACTTGCGATGTGCCTGCTATGCTTAAATACATACCTCGCCCAGTAGCTTTCCAATATGGGAACATCAATTGTACATCGGGCCGGGATAATCATGTCGTCACGGAACTTACTGTAAGAGTAATAATTCAAAAGAGATGTCGGTTTGCTAGTTCACGTAAAGGTGCCTCGCGCCACCTCTAAGTAAGTGAGCCGTCGAGACATTATCCCTGATTTTCTCACTACTATTAGTACTCACGGCGCAATACCACCACAGCCTTGTCTCGCCAGAATGCCGGTCAGCATATGGAAGAGCTCAAGGCAGGTCAATTCGCACTGTGAGGGTCACATGGGCGTTTGGCACTACCGACACGAACCTCAGTTAGCGTACATCCTACCAGAGGTCTGTGGCCCCGTGGTCAAAAGTGCGGGTTTCGTATTTGCTGCTCGTCTGTACTTTCAGAATCTTGACCTGCACGGCAAAGAGACGCTTTTTATGGAGCTCGACATGGCAACAACGCGACGGATCTACGTCACAACGAGAATAGTGTAAACGAAGCTGCTGACGGCGGAAGCGACATAGGGATCTGTGAGTTGTTATTCGCGAAAAACATCCGTCCCCGTGGGGGATAGTCACTGACGCGGTTTTGTAGAAGCCTAGGGGAACAGGTTAGTTTGACTAGCTTAAGAATGTAAATTCTGGGATTATACTGTAGTAATCACTAATTAACGGTGAGGGTTTTAAGACGGATCTTTGCAAATTCAAGCGAGGTGATTTCAACAAATTTTG'
            # l = 10
            # n = 50
            # for i in range(0, n):
            #     spectrum.append(RAW_DATA[i:i + 10])
            #     cached_distances[spectrum[i]] = {}
            #
            # print(RAW_DATA[:len(spectrum)+9])
            # print((len(spectrum)-1)*(len(spectrum[0])-1))
            #
            # size = int(2*(n**0.5))
            # print('Populacja: ' + str(size))
            # result_gene = GeneticAlgorithm(100, 1000, 0.001).compute()
            # print(result_gene.solution_value(), result_gene.get_solution(len(spectrum)+9), len(spectrum))
