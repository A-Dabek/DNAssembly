import random
import time
import numpy as np


def overlap(string1, string2):
    t = compute_back_track_table(string2)
    m, i = 0, 0
    while m + i < len(string1):
        if string2[i] == string1[m + i]:
            i += 1
        else:
            m += i - t[i]
            if i > 0:
                i = t[i]
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
        self.best_value = 0
        self.order = np.array(range(n))

    def compute(self, verbose=False):  # main function
        self.population = np.zeros(shape=(self.size, n+1), dtype=np.int16)  # generate population
        for i in range(self.size):
            self.population[i, 0:n] = np.random.permutation(n)
        self.best_gene = self.population[0]  # assume first one is the best
        no_progress = 0  # iterations without progress
        time_result = 0
        TIMESTART = time.time()
        while self.iterations_to_wait - no_progress > 0:  # run until there is no improvement
            TIMESTAMP = time.time()
            for i in range(self.size):
                total = 0
                for j in range(n-1):
                    total += cached[self.population[i, j], self.population[i, j+1]]
                self.population[i, n] = total
            # sort population basing on length of their sequence
            self.order = self.population[:, n].argsort()
            no_progress += 1  # assume there is no progress
            if self.is_improved():  # check if progress was made
                time_result = time.time() - TIMESTART
                best_value = self.assemble_fragments(self.best_gene)
                if verbose:
                    print(no_progress, best_value)
                no_progress = 0  # reset progress count
            self.create_new_generation()  # perform crossovers and add random new genes
            self.perform_mutations()  # perform mutations basing on a given chance percentage
            print(time.time() - TIMESTAMP)
        return self.assemble_fragments(self.best_gene), time_result  # return the best gene ever!

    def is_improved(self):
        # check whether present best gene uses more spectrum elements than the best gene ever
        if self.assemble_fragments(self.population[self.order[-1]])[0] > self.assemble_fragments(self.best_gene)[0]:
            self.best_gene = self.population[self.order[-1]]
            return True
        return False

    def create_new_generation(self):
        crossovers = 0  # count the number of crossovers performed
        while crossovers < self.size // 4:  # create a fourth of an original size genes in crossovers
            first_parent = second_parent = random.randint(0, self.size - 1)  # get random parents
            while first_parent == second_parent:  # ensure that parents are not duplicate elements
                second_parent = random.randint(0, self.size - 1)
            self.crossover(first_parent, second_parent, self.order[crossovers])  # perform crossover
            crossovers += 1
        for i in range(self.size // 2, crossovers, -1):  # the weakest fourth of genes replace with random ones
            ind = self.order[i]
            self.population[ind, 0:n] = np.random.permutation(n)
            self.population[ind, n] = 0

    def perform_mutations(self):
        for gene in range(self.size):
            if random.random() < self.mutation_chance:
                self.mutate(gene)

    def crossover(self, g1, g2, at):  # edge recombination operator
        # neigh_list = np.full(shape=(n, 5), fill_value=-1, dtype=np.int16)  # adjacency list
        neigh_list = {}
        # neigh_list[:, 4] = 4
        child = np.full(shape=(n, ), fill_value=-1, dtype=np.int16)
        child_size = 0

        # for i in range(n):  # add neighbours to each node
        #     neigh_list[self.population[g1, i], 0] = self.population[g1, (i-1) % n]
        #     neigh_list[self.population[g1, i], 1] = self.population[g1, (i+1) % n]
        #     neigh_list[self.population[g2, i], 2] = self.population[g2, (i-1) % n]
        #     neigh_list[self.population[g2, i], 3] = self.population[g2, (i+1) % n]
        for i in range(n):
            neigh_list[self.population[g1, i]] = {self.population[g1, (i-1) % n], self.population[g1, (i+1) % n]}
        for i in range(n):
            neigh_list[self.population[g2, i]].add(self.population[g2, (i-1) % n])
            neigh_list[self.population[g2, i]].add(self.population[g2, (i+1) % n])

        # a starting point of a child is a starting point of one of the parents
        neigh_chosen = self.population[g1, 0] if random.random() > 0.5 else self.population[g2, 0]
        child[0] = neigh_chosen
        child_size += 1

        while child_size < n:  # run until child has desired length
            # min_neigh_list = -1
            for k in neigh_list:
                if neigh_chosen in neigh_list[k]:
                    neigh_list[k].remove(neigh_chosen)
            min_neigh_list = list(neigh_list[neigh_chosen])
            del neigh_list[neigh_chosen]  # delete list of the chosen node
            if len(min_neigh_list) > 0:  # if the chosen node has any neighbours
                # get the best match out of neighbours as next
                max_overlap = cached[neigh_chosen, max(min_neigh_list, key=lambda x: cached[neigh_chosen, x])]
                possibilities = list(filter(lambda x: cached[neigh_chosen, x] == max_overlap, min_neigh_list))
                neigh_chosen = possibilities[random.randint(0, len(possibilities) - 1)]
            else:
                # get the best match out of every node as next
                max_overlap = cached[neigh_chosen, max(neigh_list, key=lambda x: cached[neigh_chosen, x])]
                possibilities = list(filter(lambda x: cached[neigh_chosen, x] == max_overlap, neigh_list))
                neigh_chosen = possibilities[random.randint(0, len(possibilities) - 1)]

            child[child_size] = neigh_chosen  # add the node to the solution
            child_size += 1
        self.population[at, 0:n] = child
        self.population[at, n] = 0

    @staticmethod
    def assemble_fragments(sequence):  # get the slice with the best sum of overlaps
        limit = n + l - 1  # maximum length of a sequence
        # count the overlap between fragments
        fits = np.zeros(shape=(len(sequence)-2, ), dtype=np.int16)
        for i in range(len(fits) - 1):
            fits[i] = cached[sequence[i], sequence[(i+1) % n]]
        r_slice = (0, 0)
        for initial in range(len(fits)):  # for every initial position
            for final in range(initial + r_slice[1] - r_slice[0], len(fits)):  # for every (not worse) final position
                total = sum(fits[initial:final])
                if (final - initial + 1) * l - total >= limit:  # check if the constraint is met
                    break
                if final - initial > r_slice[1] - r_slice[0]:  # if the slice covers more than previous
                    r_slice = (initial, final)
        assembled = str(spectrum[r_slice[0]])[3:-1]
        for s in range(r_slice[0], r_slice[1] + 1):  # assemble the fragments covered by the slice
            common = cached[sequence[s], sequence[s+1]]
            assembled += str(spectrum[sequence[s+1]][common:])[3:-1]
        return r_slice[1] - r_slice[0] + 1, r_slice, assembled[:limit]


if __name__ == '__main__':
    # dict_test = {'end/': 12,
    #              'neg/': 22,
    #              'pos/': 12,
    #              'rep/': 5}
    dict_test = {'end/': 12}
    start_path = 'inst/'
    result_file = open('results.txt', 'a')

    for k, v in dict_test.items():
        result_file.write(k[:-1] + '\n')
        for i in range(1, v + 1):
            file_name = str.format('{0}{1}{2}.bin', start_path, k, i)
            file = open(file_name)  # read an instance from file
            spectrum = np.array(file.readlines(), dtype='|S10')
            n = len(spectrum)
            l = 10
            cached = np.zeros(shape=(n, n), dtype=np.int16)
            for row in range(len(cached)):
                for cell in range(len(cached[row])):
                    cached[row, cell] = overlap(spectrum[row], spectrum[cell])
            ppl = n
            iters = (1000 - n) // 8
            print(str.format('file: {0}, n = {1}, l = {2}', file_name, n, l))
            print('Population: ' + str(ppl) + ' for ' + str(iters))
            result_gene, time_r = GeneticAlgorithm(ppl, iters, 0.25).compute(verbose=True)
            print(str.format('{0} elements out of {1} in {2} seconds',
                             result_gene[0],
                             n,
                             time_r))
            result_file.write(str.format('{0};{1};{2:.3f}\n',
                                         n,
                                         result_gene[0],
                                         time_r).replace('.', ','))
    result_file.close()
