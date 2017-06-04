import random
import copy
import time

spectrum = []  # fragments of DNA sequence
cached = {}  # cached overlapping value for dna fragments
n, l = 0, 0  # power of spectrum and length of an element

def overlap(string1, string2):
    c = cached[string1].get(string2, None)  # check if overlap was already found
    if c is not None:
        return c
    # otherwise perform algorithm
    t = compute_back_track_table(string2)
    m, i = 0, 0
    while m + i < len(string1):
        if string2[i] == string1[m + i]:
            i += 1
        else:
            m += i - t[i]
            if i > 0:
                i = t[i]
    cached[string1][string2] = i  # remember the overlap value
    return i


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
    def __init__(self, population_size, iterations_without_improvement, mutation_chance):
        self.size = population_size
        self.iterations_to_wait = iterations_without_improvement
        self.mutation_chance = mutation_chance
        self.population = []
        self.best_gene = None

    def compute(self, verbose=False):  # main function
        self.population = [Gene() for _ in range(self.size)]  # generate population
        self.best_gene = copy.deepcopy(self.population[0])  # assume first one is the best
        no_progress = 0  # iterations without progress
        TIMESTAMP = time.time()
        time_result = 0
        while self.iterations_to_wait - no_progress > 0:  # run until there is no improvement
            # sort population basing on length of their sequence
            self.population.sort(key=lambda x: x.get_value(), reverse=True)
            no_progress += 1  # assume there is no progress
            if self.is_improved():  # check if progress was made
                time_result = time.time() - TIMESTAMP
                best_value = self.best_gene.get_solution()
                if verbose:
                    print(no_progress, best_value)
                no_progress = 0  # reset progress count
            self.create_new_generation()  # perform crossovers and add random new genes
            self.perform_mutations()  # perform mutations basing on a given chance percentage
        return self.best_gene, time_result  # return the best gene ever!

    def is_improved(self):
        # check whether present best gene uses more spectrum elements than the best gene ever
        if self.population[0].get_solution()[0] > self.best_gene.get_solution()[0]:
            self.best_gene = copy.deepcopy(self.population[0])
            return True
        return False

    def create_new_generation(self):
        crossovers = 0  # count the number of crossovers performed
        while crossovers < self.size // 4:  # create a fourth of an original size genes in crossovers
            first_parent = second_parent = random.randint(0, self.size - 1)  # get random parents
            while first_parent == second_parent:  # ensure that parents are not duplicate elements
                second_parent = random.randint(0, self.size - 1)
            g1, g2 = self.population[first_parent], self.population[second_parent]
            Gene.crossover(g1, g2, self.population, self.size//2 + crossovers)  # perform crossover
            crossovers += 1
        for i in range(self.size//2 + crossovers, self.size):  # the weakest fourth of genes replace with random ones
            self.population[i] = Gene()

    def perform_mutations(self):
        for gene in self.population:
            if random.random() < self.mutation_chance:
                gene.mutate()


class Gene:
    @staticmethod
    def crossover(g1, g2, arr, ind):  # edge recombination operator
        neigh_list = {}  # adjacency list
        length = len(g1.permutation)  # expected length of a child
        for i, base in enumerate(g1.permutation):  # create nodes
            neigh_list[base] = {g1.permutation[i - 1], g1.permutation[(i + 1) % length]}
        for i, base in enumerate(g2.permutation):  # add neighbours to each node
            neigh_list[base].add(g2.permutation[i - 1])
            neigh_list[base].add(g2.permutation[(i + 1) % length])

        # a starting point of a child is a starting point of one of the parents
        neigh_chosen = [g1.permutation[0], g2.permutation[0]][random.randint(0, 1)]
        child = [neigh_chosen]

        while len(child) < length:  # run until child has desired length
            min_length = 5  # each list has lower length than 5
            min_neigh_list = []
            for k in neigh_list:  # for every node
                if neigh_chosen in neigh_list[k]:  # remove a chosen fragment from the node
                    neigh_list[k].remove(neigh_chosen)
                if k in neigh_list[neigh_chosen]:  # if a node is a neighbour of previously chosen
                    if len(neigh_list[k]) < min_length:  # remember nodes with the fewest neighbours
                        min_length = len(neigh_list[k])
                        min_neigh_list = [k]
                    elif len(neigh_list[k]) == min_length:
                        min_neigh_list.append(k)
            del neigh_list[neigh_chosen]  # delete list of the chosen node
            if len(min_neigh_list) > 0:  # if the chosen node has any neighbours
                # get the best match out of neighbours as next
                max_overlap = overlap(neigh_chosen, max(min_neigh_list, key=lambda x: overlap(neigh_chosen, x)))
                possibilities = list(filter(lambda x: overlap(neigh_chosen, x) == max_overlap, min_neigh_list))
                neigh_chosen = possibilities[random.randint(0, len(possibilities) - 1)]
            else:
                # get the best match out of every node as next
                max_overlap = overlap(neigh_chosen, max(neigh_list, key=lambda x: overlap(neigh_chosen, x)))
                possibilities = list(filter(lambda x: overlap(neigh_chosen, x) == max_overlap, neigh_list))
                neigh_chosen = possibilities[random.randint(0, len(possibilities) - 1)]

            child.append(neigh_chosen)  # add the node to the solution
        arr[ind] = Gene(child)

    def __init__(self, permutation=None):
        if permutation is None:  # if the permutation is not provided, shuffle the entire spectrum
            self.permutation = list(spectrum)
            random.shuffle(self.permutation)
        else:
            self.permutation = list(permutation)
        self.value_is_valid = False
        self.value = 0
        self.solution_is_valid = False
        self.solution = None

    def get_value(self):  # sum the total of all overlaps
        if not self.value_is_valid:
            self.value = 0
            for codon in range(len(self.permutation) - 1):
                self.value += overlap(self.permutation[codon], self.permutation[codon + 1])
            self.value_is_valid = True
        return self.value

    def get_solution(self):  # get the slice with the best sum of overlaps
        if self.solution_is_valid:
            return self.solution
        limit = n + l - 1  # maximum length of a sequence
        sequence = self.permutation[0]
        # count the overlap between fragments
        fits = [overlap(self.permutation[s], self.permutation[s + 1]) for s in range(len(self.permutation) - 1)]
        slice = (0, 0)
        for initial in range(len(fits)):  # for every initial position
            for final in range(initial + slice[1] - slice[0], len(fits)):  # for every (not worse) final position
                total = sum(fits[initial:final])
                if (final - initial + 1) * l - total >= limit:  # check if the constraint is met
                    break
                if final - initial > slice[1] - slice[0]:  # if the slice covers more than previous
                    slice = (initial, final)
        for s in range(slice[0], slice[1] + 1):  # assemble the fragments covered by the slice
            common = overlap(self.permutation[s], self.permutation[s + 1])
            sequence += self.permutation[s + 1][common:]
        self.solution = slice[1] - slice[0] + 2, slice, sequence[:limit]
        self.solution_is_valid = True
        return slice[1] - slice[0] + 2, slice, sequence[:limit]

    def mutate(self):
        solution = self.get_solution()
        seq_s, seq_f = solution[1]  # get the slice with the best result
        minimal = (-1, 0, -1, 0)  # replace index / with overlap / with index
        for i in range(seq_s, min(seq_f, len(self.permutation) - 2)):  # for every element in the solution slice
            iteration_overlap = overlap(self.permutation[i], self.permutation[i + 1]) + \
                                overlap(self.permutation[i + 1], self.permutation[i + 2])
            if iteration_overlap == 2 * (l - 1):  # if the match is perfect
                continue
            for j in range(0, seq_s):  # for every element before slice
                # check if the fragment out of the slice fits better then some one from the slice
                total_overlap = overlap(self.permutation[i], self.permutation[j]) + overlap(self.permutation[j],
                                                                                            self.permutation[i + 2])
                if minimal[1] < total_overlap and total_overlap > iteration_overlap:
                    minimal = (j, total_overlap, i + 1, iteration_overlap)
            for j in range(seq_f, len(self.permutation)):  # for every element after slice do the same
                total_overlap = overlap(self.permutation[i], self.permutation[j]) + overlap(self.permutation[j],
                                                                                            self.permutation[i + 2])
                if minimal[1] < total_overlap and total_overlap > iteration_overlap:
                    minimal = (j, total_overlap, i + 1, iteration_overlap)
            if minimal[1] >= (l - 1) * 2:  # if a perfect replacement was found
                break
        if minimal[0] >= 0:  # if any fragment actually fits better, swap them
            rnd1 = minimal[0]
            rnd2 = minimal[2]
            temp = self.permutation[rnd1]
            self.permutation[rnd1] = self.permutation[rnd2]
            self.permutation[rnd2] = temp
            self.value_is_valid = False
            self.solution_is_valid = False


if __name__ == '__main__':
    dict_test = {'end/': 12,
                 'neg/': 22,
                 'pos/': 12,
                 'rep/': 5}
    start_path = 'inst/'
    result_file = open('results.txt', 'a')

    for k, v in dict_test.items():
        result_file.write(k[:-1] + '\n')
        for i in range(1, v + 1):
            file_name = str.format('{0}{1}{2}.bin', start_path, k, i)
            cached = {}
            spectrum = []
            file = open(file_name)  # read an instance from file
            line = file.readline().replace('\n', '')
            while len(line) > 1:
                spectrum.append(line)
                cached[line] = {}  # prepare cached overlap structure
                line = file.readline().replace('\n', '')
            n = len(spectrum)
            l = len(spectrum[0])
            ppl = n
            iters = (1000-n)//4
            print(str.format('file: {0}, n = {1}, l = {2}', file_name, n, l))
            print('Population: ' + str(ppl) + ' for ' + str(iters))
            result_gene, time_r = GeneticAlgorithm(ppl, iters, 0.25).compute(verbose=True)
            print(str.format('{0} elements out of {1} in {2} seconds',
                             result_gene.get_solution()[0],
                             n,
                             time_r))
            result_file.write(str.format('{0};{1};{2:.3f}\n',
                                         n,
                                         result_gene.get_solution()[0],
                                         time_r).replace('.', ','))
    result_file.close()
